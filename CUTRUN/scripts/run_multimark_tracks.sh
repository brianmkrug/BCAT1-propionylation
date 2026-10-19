#!/usr/bin/env bash
# =============================================================================
# 50 bp CPM tracks for H3K27ac, ENL and H3K9pr — same treatment as the
# H3K4me3/Pol2 tracks built by scripts/run_phase1_tracks.sh, so that promoter
# deltas are comparable across marks.
#
# Pipeline per sample (identical to phase 1 stages 1/3/4/5, fused for disk):
#   chr-filter (+ GenBank->UCSC rename for H3K27ac) + MAPQ30  ->  depth-match
#   ->  bamCoverage CPM 50 bp, blacklisted  ->  (i) IP minus matched IgG
#                                              (ii) shBCAT1 vs shNT log2
#
# Differences from run_phase1_tracks.sh, and why:
#  * stages 1 and 3 are FUSED into one `samtools view` (chr-filter, header
#    rebuild, rename and subsample in a single pass). Phase 1 materialised a
#    24 GB `01_main/` intermediate; this machine has ~14 GiB free, so the
#    un-downsampled copy is never written. Nothing else changes — `-s` selects
#    on a hash of the read NAME, so mates stay together either way.
#  * fragment totals for depth matching are read from the 10 kb count pass
#    (analysis_QC_normalization/multimark_2026-07/counts/) instead of a fresh
#    `samtools view -c`. Same filters (-f 66 -F 1804 -q 30, chr1-22 + chrX),
#    so the same number — just already paid for.
#  * `-q 30` is applied here. Every other library in the project was MAPQ30-
#    filtered at alignment; the two H3K27ac r1 libraries were filtered at
#    `-bq 1` only, and would otherwise carry multi-mapper signal the rest of
#    the panel does not have.
#  * ds BAMs are deleted once their bigWig exists (KEEP_DS=1 to retain). They are
#    exactly reproducible from the fragment totals + SEED below.
#
# RUN DETACHED, NOT INSIDE CLAUDE CODE:
#   cd ${BCAT1_REPO:-.}/CUTRUN
#   nohup ./scripts/run_multimark_tracks.sh > logs/multimark_$(date +%F_%H%M).log 2>&1 &
# =============================================================================
set -euo pipefail

ROOT="${ROOT:-$HOME/Desktop/BCAT1/BCAT1/CUTRUN}"
BAMDIR="${BAMDIR:-$ROOT/hg38}"
QCDIR="${QCDIR:-$ROOT/analysis_QC_normalization/multimark_2026-07}"
OUT="${OUT:-$ROOT/multimark_tracks}"
PHASE1_BW="${PHASE1_BW:-$ROOT/H3K4me3/phase1_tracks/03_bw}"   # existing H3K4me3 + IgG CPM tracks
BLACKLIST_GZ="${BLACKLIST_GZ:-$ROOT/hg38-blacklist.v2.bed.gz}"

BINSIZE="${BINSIZE:-50}"
SEED="${SEED:-42}"                # matches run_phase1_tracks.sh — do not change
KEEP_DS="${KEEP_DS:-0}"
TIER_TOL="${TIER_TOL:-1.05}"      # build a replicate-blocked (tierB) track only
                                  # when the pair minimum exceeds the group
                                  # minimum by more than this factor

ENVBIN="${ENVBIN:-$HOME/miniconda3/envs/k4width/bin}"
[[ -x "$ENVBIN/samtools" ]] || { echo "[fatal] no samtools at $ENVBIN" >&2; exit 1; }
export PATH="$ENVBIN:$PATH"

if PCORES=$(sysctl -n hw.perflevel0.logicalcpu 2>/dev/null); then :; else
  PCORES=$(sysctl -n hw.logicalcpu 2>/dev/null || nproc)
fi
JOBS="${JOBS:-3}"
THREADS="${THREADS:-$(( PCORES / JOBS ))}"; (( THREADS < 2 )) && THREADS=2
echo "[cfg] P-cores=$PCORES jobs=$JOBS threads/job=$THREADS bin=$BINSIZE seed=$SEED"
echo "[cfg] free space: $(df -h "$ROOT" | awk 'NR==2{print $4}')"

mkdir -p "$OUT"/{02_ds,03_bw,04_cmp,qc} "$ROOT/logs"
BLACKLIST="$OUT/hg38-blacklist.v2.bed"
[[ -s "$BLACKLIST" ]] || gunzip -c "$BLACKLIST_GZ" > "$BLACKLIST"

# ----------------------------------------------------------------- manifest --
# label|mark|condition|rep|naming|source basename (no .bam)
SAMPLES=(
  "shNT-H3K27ac-r1|H3K27ac|shNT|r1|gb|shNT-H3K27ac-r1_nodup_sorted_unique"
  "shNT-H3K27ac-r2|H3K27ac|shNT|r2|gb|shNT-H3K27ac-r2_nodup_sorted_unique"
  "shBCAT1-H3K27ac-r1|H3K27ac|shBCAT1|r1|gb|shBCAT1-H3K27ac-r1_nodup_sorted_unique"
  "shBCAT1-H3K27ac-r2|H3K27ac|shBCAT1|r2|gb|shBCAT1-H3K27ac-r2_nodup_sorted_unique"
  "shNT-ENL-r1|ENL|shNT|r1|ucsc|DIPG13-shNT-ENL-r1_nodup_sorted_unique"
  "shNT-ENL-r2|ENL|shNT|r2|ucsc|DIPG13-shNT-ENL-r2_nodup_sorted_unique"
  "shBCAT1-ENL-r1|ENL|shBCAT1|r1|ucsc|DIPG13-shBCAT1-ENL-r1_nodup_sorted_unique"
  "shBCAT1-ENL-r2|ENL|shBCAT1|r2|ucsc|DIPG13-shBCAT1-ENL-r2_nodup_sorted_unique"
  "shNT-H3K9pr-r1|H3K9pr|shNT|r1|ucsc|DIPG13-shNT-H3K9pr-r1_nodup_sorted_unique"
  "shNT-H3K9pr-r2|H3K9pr|shNT|r2|ucsc|DIPG13-shNT-H3K9pr-r2_nodup_sorted_unique"
  "shBCAT1-H3K9pr-r1|H3K9pr|shBCAT1|r1|ucsc|DIPG13-shBCAT1-H3K9pr-r1_nodup_sorted_unique"
  "shBCAT1-H3K9pr-r2_60M|H3K9pr|shBCAT1|r2|ucsc|DIPG13-shBCAT1-H3K9pr-r2_60M_nodup_sorted_unique"
)
MARKS=(H3K27ac ENL H3K9pr)

GB2UCSC="CM000663.2:chr1 CM000664.2:chr2 CM000665.2:chr3 CM000666.2:chr4 \
CM000667.2:chr5 CM000668.2:chr6 CM000669.2:chr7 CM000670.2:chr8 CM000671.2:chr9 \
CM000672.2:chr10 CM000673.2:chr11 CM000674.2:chr12 CM000675.2:chr13 \
CM000676.2:chr14 CM000677.2:chr15 CM000678.2:chr16 CM000679.2:chr17 \
CM000680.2:chr18 CM000681.2:chr19 CM000682.2:chr20 CM000683.2:chr21 \
CM000684.2:chr22 CM000685.2:chrX"

# ============================================================ STAGE 0 counts ==
# Fragment totals from the 10 kb pass. Fail loudly rather than silently
# depth-matching to a missing (=0) target.
COUNTS="$OUT/qc/fragment_counts.tsv"
if [[ ! -s "$COUNTS" ]]; then
  tmp="$COUNTS.partial"; : > "$tmp"
  for s in "${SAMPLES[@]}"; do
    IFS='|' read -r label mark cond rep naming src <<< "$s"
    f="$QCDIR/counts/${label}.10kb.tsv.gz"
    [[ -s "$f" ]] || { echo "[fatal] no 10 kb counts for $label ($f) — run scripts/count_10kb_bins.sh first" >&2; exit 1; }
    n=$(gzip -dc "$f" | awk '{t+=$3} END{printf "%d", t}')
    (( n > 0 )) || { echo "[fatal] zero fragments for $label" >&2; exit 1; }
    printf '%s\t%s\t%s\t%s\t%d\n' "$label" "$mark" "$cond" "$rep" "$n" >> "$tmp"
    echo "[stage0] $label -> $n fragments"
  done
  mv "$tmp" "$COUNTS"
fi
column -t "$COUNTS"

get_n()   { awk -F'\t' -v l="$1" '$1==l{print $5}' "$COUNTS"; }
labels_of_mark() { awk -F'\t' -v m="$1" '$2==m{print $1}' "$COUNTS"; }
labels_of_pair() { awk -F'\t' -v m="$1" -v r="$2" '$2==m && $4==r{print $1}' "$COUNTS"; }

min_of() {
  local m="" l n
  for l in "$@"; do
    n=$(get_n "$l")
    [[ -n "$n" ]] || { echo "[fatal] no fragment count for '$l'" >&2; return 1; }
    if [[ -z "$m" ]] || (( n < m )); then m=$n; fi
  done
  echo "$m"
}

# ======================================================== STAGE 1 downsample ==
# One pass: region-restrict, rebuild header (drop non-primary @SQ), rename
# GenBank accessions, MAPQ-filter, subsample. Regions are emitted in header
# order so the result stays coordinate-sorted.
make_ds() {                        # $1 label  $2 naming  $3 src  $4 target  $5 tag
  local label=$1 naming=$2 src=$3 target=$4 tag=$5
  local bam="$BAMDIR/${src}.bam"
  local out="$OUT/02_ds/${label}.${tag}.bam"
  local bw="$OUT/03_bw/${label}.${tag}.cpm.bw"
  [[ -s "$bw" ]]      && { echo "[skip] ds $label.$tag (bigWig exists)"; return; }
  [[ -s "$out.bai" ]] && { echo "[skip] ds $label.$tag"; return; }

  local cur frac keep map
  cur=$(get_n "$label")
  frac=$(awk -v t="$target" -v c="$cur" 'BEGIN{printf "%.6f", (t<c? t/c : 1)}')

  # NB: awk -v cannot carry a literal newline, so these lists are space-separated.
  if [[ $naming == gb ]]; then
    keep=$(for p in $GB2UCSC; do printf '%s ' "${p%%:*}"; done)
    map="$GB2UCSC"
  else
    keep=$(printf 'chr%s ' {1..22} X)
    map=""
  fi
  # region list in header order -> output stays coordinate-sorted
  local regions
  regions=$(samtools view -H "$bam" \
    | awk -v K="$keep" 'BEGIN{n=split(K,a," "); for(i=1;i<=n;i++) if(a[i]!="") k[a[i]]=1}
          $1=="@SQ"{ for(i=2;i<=NF;i++) if($i ~ /^SN:/){ s=substr($i,4); if(s in k) print s } }' \
    | tr '\n' ' ')
  [[ -n "${regions// /}" ]] || { echo "[fatal] no primary contigs matched in $src header" >&2; return 1; }

  local sarg=""
  if awk -v f="$frac" 'BEGIN{exit !(f<0.999)}'; then sarg="-s ${SEED}.${frac#0.}"; fi
  echo "[stage1] $label.$tag  $cur -> $target  (frac=$frac)  naming=$naming"

  # shellcheck disable=SC2086
  samtools view -h -q 30 $sarg "$bam" $regions \
    | awk -v MAP="$map" -v REG="$regions" 'BEGIN{FS=OFS="\t"
             n=split(MAP,kv," "); for(i=1;i<=n;i++){split(kv[i],p,":"); m[p[1]]=p[2]}
             n=split(REG,r," ");  for(i=1;i<=n;i++) if(r[i]!="") keep[r[i]]=1 }
           /^@SQ/ { sn=""; for(i=2;i<=NF;i++) if($i ~ /^SN:/){ sn=substr($i,4);
                      if(sn in m) $i="SN:" m[sn] }
                    if(sn=="" || !(sn in keep)) next
                    print; next }
           /^@/   { print; next }
           { $3 = ($3 in m) ? m[$3] : $3
             # A mate on a contig we just dropped from the header would make
             # samtools warn on every such line ("unrecognized mate reference").
             # Those reads are never proper pairs, so they are excluded by
             # --samFlagInclude 66 downstream; blank the mate reference rather
             # than emit an unresolvable name.
             if($7 == "=")            { }
             else if($7 in m)         { $7 = m[$7] }
             else if(!($7 in keep))   { $7 = "*"; $8 = 0 }
             print }' \
    | samtools view -b -@ "$THREADS" -o "$out.partial" -
  mv "$out.partial" "$out"
  samtools index -@ "$THREADS" "$out"
}
# get_n too: make_ds calls it inside the xargs subshell, where an unexported
# function would be silently undefined and the fragment total would come back empty.
export -f make_ds get_n
export OUT BAMDIR THREADS SEED GB2UCSC COUNTS

# Decide, per mark, which tiers are needed.
#   tierA : all four libraries of the mark matched to the group minimum
#           (absolute / pooled comparisons, and the minus-IgG signal tracks)
#   tierB : matched within replicate pair — only when that buys real depth
#           back (pair min > TIER_TOL x group min), for replicate-blocked deltas
# NB: macOS ships bash 3.2 -- no associative arrays, no mapfile. tierB pairs are
# tracked in a plain space-delimited string.
JOBSPEC="$OUT/qc/ds_jobs.tsv"; : > "$JOBSPEC"
TIERB_PAIRS=" "
has_tierb() { case "$TIERB_PAIRS" in *" $1:$2 "*) return 0;; *) return 1;; esac; }

for mark in "${MARKS[@]}"; do
  labs=$(labels_of_mark "$mark")
  # shellcheck disable=SC2086
  gmin=$(min_of $labs)
  echo "[tier] $mark group min = $gmin"
  for l in $labs; do
    src=$(printf '%s\n' "${SAMPLES[@]}" | awk -F'|' -v l="$l" '$1==l{print $6}')
    nm=$(printf '%s\n'  "${SAMPLES[@]}" | awk -F'|' -v l="$l" '$1==l{print $5}')
    printf '%s\t%s\t%s\t%s\t%s\n' "$l" "$nm" "$src" "$gmin" tierA >> "$JOBSPEC"
  done
  for rep in r1 r2; do
    pl=$(labels_of_pair "$mark" "$rep")
    npl=$(printf '%s\n' $pl | grep -c . || true)
    [[ $npl -eq 2 ]] || { echo "[warn] $mark $rep: $npl libraries, no tierB"; continue; }
    # shellcheck disable=SC2086
    pmin=$(min_of $pl)
    if awk -v p="$pmin" -v g="$gmin" -v t="$TIER_TOL" 'BEGIN{exit !(p > t*g)}'; then
      echo "[tier] $mark $rep pair min = $pmin (> ${TIER_TOL}x group min) -> tierB"
      TIERB_PAIRS="$TIERB_PAIRS$mark:$rep "
      for l in $pl; do
        src=$(printf '%s\n' "${SAMPLES[@]}" | awk -F'|' -v l="$l" '$1==l{print $6}')
        nm=$(printf '%s\n'  "${SAMPLES[@]}" | awk -F'|' -v l="$l" '$1==l{print $5}')
        printf '%s\t%s\t%s\t%s\t%s\n' "$l" "$nm" "$src" "$pmin" tierB >> "$JOBSPEC"
      done
    else
      echo "[tier] $mark $rep pair min = $pmin (<= ${TIER_TOL}x group min) -> tierA serves as tierB"
    fi
  done
done

# Largest target first: the big downsampled BAMs are built while the disk is
# still empty, and the small ones land at the end once the bigWigs have taken
# their share. Sorting is not cosmetic here -- see the disk note below.
sort -t$'\t' -k4,4nr "$JOBSPEC" -o "$JOBSPEC"
echo "[plan] $(wc -l < "$JOBSPEC" | tr -d ' ') track jobs (largest first):"
column -t "$JOBSPEC"
[[ "${PLAN_ONLY:-0}" == 1 ]] && { echo "[plan-only] stopping before any work"; exit 0; }

# ========================================================== STAGE 2 coverage ==
cover() {                          # $1 stem (label.tag)
  local stem=$1
  local bam="$OUT/02_ds/${stem}.bam"
  local bw="$OUT/03_bw/${stem}.cpm.bw"
  [[ -s "$bw" ]] && { echo "[skip] cover $stem"; return; }
  [[ -s "$bam" ]] || { echo "[fatal] missing $bam" >&2; return 1; }
  echo "[stage2] $stem"
  bamCoverage -b "$bam" -o "$bw.partial" \
    --binSize "$BINSIZE" --normalizeUsing CPM \
    --extendReads --samFlagInclude 66 --samFlagExclude 1804 \
    --blackListFileName "$BLACKLIST" \
    -p "$THREADS"
  mv "$bw.partial" "$bw"
  if [[ "$KEEP_DS" != "1" ]]; then rm -f "$bam" "$bam.bai"; echo "[clean] removed $bam"; fi
}
export -f cover; export BINSIZE BLACKLIST KEEP_DS

# ---------------------------------------------------- fused downsample+cover --
# Downsampling every sample first would need ~23 GB of intermediate BAM at once;
# this volume has ~16 GiB free. Fusing the two stages means a downsampled BAM
# only ever exists between its own creation and its own bigWig, so peak usage is
# JOBS intermediates instead of 16.
build() {                          # $1 label $2 naming $3 src $4 target $5 tag
  local free_kb
  free_kb=$(df -k "$OUT" | awk 'NR==2{print $4}')
  if [[ $free_kb -lt 5242880 ]]; then      # 5 GiB
    echo "[fatal] only $((free_kb/1048576)) GiB free — refusing to start $1.$5" >&2
    return 1
  fi
  make_ds "$1" "$2" "$3" "$4" "$5" || return 1
  cover "$1.$5"
}
export -f build

awk -F'\t' '{print $1"\t"$2"\t"$3"\t"$4"\t"$5}' "$JOBSPEC" \
  | xargs -P "$JOBS" -n5 bash -c 'build "$0" "$1" "$2" "$3" "$4"'

# ========================================================== STAGE 3 compare ==
# IgG CPM tracks come from phase 1 (both IgG libraries matched to their own
# group minimum there). Symlinked in so 03_bw is a single self-contained set.
for c in shNT shBCAT1; do
  src="$PHASE1_BW/${c}-IgG-r1.ds.cpm.bw"
  if [[ -s "$src" ]]; then ln -sf "$src" "$OUT/03_bw/${c}-IgG-r1.ds.cpm.bw"
  else echo "[warn] missing $src — minus-IgG tracks will be skipped"; fi
done
for f in "$PHASE1_BW"/*H3K4me3*.cpm.bw; do
  [[ -e "$f" ]] && ln -sf "$f" "$OUT/03_bw/$(basename "$f")"
done

cmp2() {                           # $1 bw1  $2 bw2  $3 op  $4 out
  local out="$OUT/04_cmp/$4"
  [[ -s "$out" ]] && { echo "[skip] cmp $4"; return; }
  for b in "$1" "$2"; do
    [[ -e "$OUT/03_bw/$b" ]] || { echo "[warn] cmp $4: missing $b — skipped"; return 0; }
  done
  echo "[stage3] $4"
  bigwigCompare -b1 "$OUT/03_bw/$1" -b2 "$OUT/03_bw/$2" \
    --operation "$3" --binSize "$BINSIZE" --pseudocount 1 \
    --skipZeroOverZero --blackListFileName "$BLACKLIST" \
    -p "$THREADS" -o "$out.partial"
  mv "$out.partial" "$out"
}

# (i) IP minus matched IgG — tierA, i.e. all libraries of a mark at equal depth
while IFS=$'\t' read -r label naming src target tag; do
  [[ $tag == tierA ]] || continue
  c=shNT; [[ $label == shBCAT1* ]] && c=shBCAT1
  cmp2 "${label}.tierA.cpm.bw" "${c}-IgG-r1.ds.cpm.bw" subtract "${label}.tierA.minusIgG.bw"
done < <(sort -u "$JOBSPEC")

# (ii) shBCAT1 vs shNT within mark and replicate — the promoter-delta input.
# Copy number is shared between conditions and cancels in this ratio.
for mark in "${MARKS[@]}"; do
  for rep in r1 r2; do
    tag=tierA; has_tierb "$mark" "$rep" && tag=tierB
    nt=$(labels_of_pair "$mark" "$rep" | grep '^shNT'    || true)
    bc=$(labels_of_pair "$mark" "$rep" | grep '^shBCAT1' || true)
    [[ -n $nt && -n $bc ]] || { echo "[warn] $mark $rep incomplete pair — no delta"; continue; }
    cmp2 "${bc}.${tag}.cpm.bw" "${nt}.${tag}.cpm.bw" log2 "delta_${mark}_${rep}.log2.bw"
  done
done

echo
echo "[done] CPM tracks in $OUT/03_bw ; comparisons in $OUT/04_cmp"
echo "[note] ds BAMs $( [[ "$KEEP_DS" == 1 ]] && echo retained || echo 'deleted (KEEP_DS=1 to keep)' )"
df -h "$ROOT" | awk 'NR==2{print "[disk] free: "$4}'
