#!/usr/bin/env bash
# =============================================================================
# Genome-wide 50 bp tracks: replicate-POOLED (all available depth), CPM, then
# spike-in scaled.  10 pooled tracks + 4 delta tracks.
#
#   ENL x0.80   H3K9pr x0.62   H3K27ac x1.00   H3K4me3 x1.00 (assumed)  IgG x1.00
#
# HOW THE POOLING IS DONE, AND WHY IT LOOKS INDIRECT
# --------------------------------------------------
# bamCoverage takes one BAM, so the obvious route is `samtools merge` then
# bamCoverage. That would materialise a ~7 GB merged BAM for ENL-shBCAT1 alone,
# on a volume with ~23 GiB free, on top of the filtered inputs.
#
# Instead each replicate is written as RAW coverage pre-scaled by
#       s_i = (1e6 / N_pool) * spikein
# and the replicates are then summed with `bigwigCompare --operation add`.
# That is algebraically identical to pooling the reads:
#       sum_i  c_i(bin) * 1e6/N_pool * f  ==  (sum_i c_i(bin)) / N_pool * 1e6 * f
#                                        ==  pooled CPM * spike-in factor
# ...with no merged BAM ever existing. N_pool is the summed fragment total of the
# replicates on chr1-22,X at MAPQ30 -- the same number used everywhere else in
# this analysis, read from the promoter/category count pass.
#
# The per-replicate intermediates are in pooling units (already divided by the
# POOL depth, not their own), so they are NOT meaningful on their own and are
# deleted once the pooled track exists. KEEP_PARTS=1 retains them.
#
# RUN DETACHED:
#   cd ${BCAT1_REPO:-.}/CUTRUN
#   nohup ./scripts/run_pooled_spikein_tracks.sh > logs/pooled_$(date +%F_%H%M).log 2>&1 &
# =============================================================================
set -euo pipefail

ROOT="${ROOT:-$HOME/Desktop/BCAT1/BCAT1/CUTRUN}"
BAMDIR="${BAMDIR:-$ROOT/hg38}"
QC="${QC:-$ROOT/analysis_QC_normalization/multimark_2026-07}"
OUT="${OUT:-$ROOT/multimark_tracks/pooled}"
BLACKLIST_GZ="${BLACKLIST_GZ:-$ROOT/hg38-blacklist.v2.bed.gz}"

BINSIZE="${BINSIZE:-50}"
KEEP_PARTS="${KEEP_PARTS:-0}"
# Pseudocount for the delta tracks. Genome-wide mean coverage here is ~0.07 CPM,
# so the phase-1 value of 1 would make every low-coverage bin almost pure
# pseudocount. 0.05 is chosen at the scale of the data: big enough that empty
# regions do not explode, small enough that real signal is not compressed.
PSEUDO="${PSEUDO:-0.05}"

ENVBIN="${ENVBIN:-$HOME/miniconda3/envs/k4width/bin}"
[[ -x "$ENVBIN/samtools" ]] || { echo "[fatal] no samtools at $ENVBIN" >&2; exit 1; }
export PATH="$ENVBIN:$PATH"

if PCORES=$(sysctl -n hw.perflevel0.logicalcpu 2>/dev/null); then :; else
  PCORES=$(sysctl -n hw.logicalcpu 2>/dev/null || nproc)
fi
JOBS="${JOBS:-2}"
THREADS="${THREADS:-$(( PCORES / JOBS ))}"; (( THREADS < 2 )) && THREADS=2
echo "[cfg] jobs=$JOBS threads/job=$THREADS bin=$BINSIZE pseudocount=$PSEUDO"
echo "[cfg] free: $(df -h "$ROOT" | awk 'NR==2{print $4}')"

mkdir -p "$OUT"/{parts,tmp} "$ROOT/logs"
BLACKLIST="$OUT/hg38-blacklist.v2.bed"
[[ -s "$BLACKLIST" ]] || gunzip -c "$BLACKLIST_GZ" > "$BLACKLIST"

GB2UCSC="CM000663.2:chr1 CM000664.2:chr2 CM000665.2:chr3 CM000666.2:chr4 \
CM000667.2:chr5 CM000668.2:chr6 CM000669.2:chr7 CM000670.2:chr8 CM000671.2:chr9 \
CM000672.2:chr10 CM000673.2:chr11 CM000674.2:chr12 CM000675.2:chr13 \
CM000676.2:chr14 CM000677.2:chr15 CM000678.2:chr16 CM000679.2:chr17 \
CM000680.2:chr18 CM000681.2:chr19 CM000682.2:chr20 CM000683.2:chr21 \
CM000684.2:chr22 CM000685.2:chrX"

spike_of () {   # $1 mark  $2 condition
  [[ $2 == shBCAT1 ]] || { echo 1.0; return; }
  case "$1" in
    ENL)     echo 0.80 ;;
    H3K9pr)  echo 0.62 ;;
    *)       echo 1.00 ;;    # H3K27ac measured invariant; H3K4me3 assumed
  esac
}

frag_total () { awk -F'\t' '$1=="total"{print $2}' "$QC/promcat/$1.category.tsv"; }

MANIFEST="$QC/sample_manifest.tsv"
[[ -s "$MANIFEST" ]] || { echo "[fatal] missing $MANIFEST" >&2; exit 1; }

# ---- pool depths -------------------------------------------------------------
POOLSPEC="$OUT/pool_depths.tsv"
if [[ ! -s "$POOLSPEC" ]]; then
  : > "$POOLSPEC.partial"
  for mark in H3K4me3 H3K27ac H3K9pr ENL IgG; do
    for cond in shNT shBCAT1; do
      labs=$(awk -F'\t' -v m="$mark" -v c="$cond" 'NR>1 && $2==m && $3==c{print $1}' "$MANIFEST")
      [[ -n "$labs" ]] || { echo "[fatal] no libraries for $mark $cond" >&2; exit 1; }
      n=0
      for l in $labs; do
        t=$(frag_total "$l")
        [[ -n "$t" ]] || { echo "[fatal] no fragment total for $l" >&2; exit 1; }
        n=$(( n + t ))
      done
      printf '%s\t%s\t%s\t%d\t%s\n' "$mark" "$cond" "$(echo $labs | tr ' ' ',')" \
             "$n" "$(spike_of "$mark" "$cond")" >> "$POOLSPEC.partial"
    done
  done
  mv "$POOLSPEC.partial" "$POOLSPEC"
fi
echo "[pools] mark cond replicates pooled_fragments spikein"
column -t "$POOLSPEC"

# ---- per-replicate scaled coverage -------------------------------------------
part_one () {          # $1 label $2 naming $3 src $4 scalefactor
  local label=$1 naming=$2 src=$3 sf=$4
  local bw="$OUT/parts/${label}.scaled.bw"
  [[ -s "$bw" ]] && { echo "[skip] part $label"; return; }
  local bam="$BAMDIR/${src}"
  local tmp="$OUT/tmp/${label}.main.bam"

  local keep map regions
  if [[ $naming == gb ]]; then
    keep=$(for p in $GB2UCSC; do printf '%s ' "${p%%:*}"; done); map="$GB2UCSC"
  else
    keep=$(printf 'chr%s ' {1..22} X); map=""
  fi
  regions=$(samtools view -H "$bam" \
    | awk -v K="$keep" 'BEGIN{n=split(K,a," "); for(i=1;i<=n;i++) if(a[i]!="") k[a[i]]=1}
          $1=="@SQ"{ for(i=2;i<=NF;i++) if($i ~ /^SN:/){ s=substr($i,4); if(s in k) print s } }' \
    | tr '\n' ' ')

  if [[ ! -s "$tmp.bai" ]]; then
    echo "[filter] $label"
    # shellcheck disable=SC2086
    samtools view -h -q 30 "$bam" $regions \
      | awk -v MAP="$map" -v REG="$regions" 'BEGIN{FS=OFS="\t"
             n=split(MAP,kv," "); for(i=1;i<=n;i++){split(kv[i],p,":"); m[p[1]]=p[2]}
             n=split(REG,r," ");  for(i=1;i<=n;i++) if(r[i]!="") keep[r[i]]=1 }
           /^@SQ/ { sn=""; for(i=2;i<=NF;i++) if($i ~ /^SN:/){ sn=substr($i,4);
                      if(sn in m) $i="SN:" m[sn] }
                    if(sn=="" || !(sn in keep)) next
                    print; next }
           /^@/   { print; next }
           { $3 = ($3 in m) ? m[$3] : $3
             if($7 == "=")          { }
             else if($7 in m)       { $7 = m[$7] }
             else if(!($7 in keep)) { $7 = "*"; $8 = 0 }
             print }' \
      | samtools view -b -@ "$THREADS" -o "$tmp.partial" -
    mv "$tmp.partial" "$tmp"
    samtools index -@ "$THREADS" "$tmp"
  fi

  echo "[coverage] $label  scaleFactor=$sf"
  bamCoverage -b "$tmp" -o "$bw.partial" \
    --binSize "$BINSIZE" --normalizeUsing None --scaleFactor "$sf" \
    --extendReads --samFlagInclude 66 --samFlagExclude 1804 \
    --blackListFileName "$BLACKLIST" -p "$THREADS"
  mv "$bw.partial" "$bw"
  rm -f "$tmp" "$tmp.bai"
  echo "[clean] removed $tmp"
}
export -f part_one
export OUT BAMDIR THREADS BINSIZE BLACKLIST GB2UCSC

JOBSPEC="$OUT/part_jobs.tsv"; : > "$JOBSPEC"
while IFS=$'\t' read -r mark cond labs npool spike; do
  for l in $(echo "$labs" | tr ',' ' '); do
    src=$(awk -F'\t' -v l="$l" 'NR>1 && $1==l{print $6}' "$MANIFEST")
    nm=$(awk -F'\t'  -v l="$l" 'NR>1 && $1==l{print $5}' "$MANIFEST")
    sf=$(awk -v n="$npool" -v f="$spike" 'BEGIN{printf "%.10f", 1e6/n*f}')
    printf '%s\t%s\t%s\t%s\n' "$l" "$nm" "$src" "$sf" >> "$JOBSPEC"
  done
done < "$POOLSPEC"
echo "[plan] $(wc -l < "$JOBSPEC" | tr -d ' ') per-replicate coverage jobs"
column -t "$JOBSPEC"
[[ "${PLAN_ONLY:-0}" == 1 ]] && { echo "[plan-only] stopping"; exit 0; }

# biggest BAMs first, while the disk is empty
sort -t$'\t' -k4,4g "$JOBSPEC" \
  | xargs -P "$JOBS" -n4 bash -c 'part_one "$0" "$1" "$2" "$3"'

# ---- sum replicates into the pooled track ------------------------------------
while IFS=$'\t' read -r mark cond labs npool spike; do
  out="$OUT/${cond}-${mark}.pooled.spikein.bw"
  [[ -s "$out" ]] && { echo "[skip] pool $cond-$mark"; continue; }
  set -- $(echo "$labs" | tr ',' ' ')
  if [[ $# -eq 1 ]]; then
    cp "$OUT/parts/$1.scaled.bw" "$out"          # single library: already scaled
    echo "[pool] $cond-$mark  <- $1 (single library)"
  else
    echo "[pool] $cond-$mark  <- $1 + $2"
    bigwigCompare -b1 "$OUT/parts/$1.scaled.bw" -b2 "$OUT/parts/$2.scaled.bw" \
      --operation add --binSize "$BINSIZE" \
      --blackListFileName "$BLACKLIST" -p "$THREADS" -o "$out.partial"
    mv "$out.partial" "$out"
  fi
done < "$POOLSPEC"

# ---- delta tracks ------------------------------------------------------------
for mark in H3K4me3 H3K27ac H3K9pr ENL; do
  out="$OUT/delta_${mark}.pooled.spikein.log2.bw"
  [[ -s "$out" ]] && { echo "[skip] delta $mark"; continue; }
  a="$OUT/shBCAT1-${mark}.pooled.spikein.bw"; b="$OUT/shNT-${mark}.pooled.spikein.bw"
  [[ -s "$a" && -s "$b" ]] || { echo "[warn] $mark: pooled track missing, no delta"; continue; }
  echo "[delta] $mark"
  bigwigCompare -b1 "$a" -b2 "$b" --operation log2 --pseudocount "$PSEUDO" \
    --binSize "$BINSIZE" --skipZeroOverZero \
    --blackListFileName "$BLACKLIST" -p "$THREADS" -o "$out.partial"
  mv "$out.partial" "$out"
done

if [[ "$KEEP_PARTS" != "1" ]]; then rm -rf "$OUT/parts"; echo "[clean] removed per-replicate parts"; fi
rmdir "$OUT/tmp" 2>/dev/null || true

echo
echo "[done] pooled spike-in tracks in $OUT"
ls -la "$OUT"/*.bw
df -h "$ROOT" | awk 'NR==2{print "[disk] free: "$4}'
