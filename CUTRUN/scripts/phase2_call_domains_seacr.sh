#!/usr/bin/env bash
# =============================================================================
# PHASE 2 (SEACR) — interpretation domains per mark. Replaces the MACS2 version.
#
# WHY SEACR, AND WHY DEPTH-CAPPED
# --------------------------------
# SEACR is the right caller for CUT&RUN sparse signal. It is NOT, however,
# cheaper at full pooled depth -- measured on chr20, H3K4me3 at full depth covers
# 39% of the chromosome in 120,405 contiguous islands (~5.6M genome-wide), and
# SEACR's R stage does read.table() on one row PER ISLAND. ENL at 52M fragments
# is worse still (~70% covered). Beyond the cost, that is outside SEACR's design
# regime: its block model assumes gapped, sparse coverage, and at 39-70%
# occupancy the islands merge into the background.
#
# So the peak-calling input is capped at PEAK_N fragments per pool (default 5M),
# where coverage falls to a few percent -- sparse, gapped, SEACR's native case,
# and a few minutes of R instead of tens.
#
#   Measured (chr20, H3K4me3):  full 39.0% covered / 120,405 islands
#                               25%  15.3% covered /  60,492 islands
#                                5%   4.0% covered /  19,785 islands
#
# This costs nothing that matters: domains are used ONLY to decide membership --
# which regions a delta is interpreted on. Every reported magnitude still comes
# from the full-depth counts (Phase 3). Peak fidelity is not the deliverable.
#
# Subsampling is applied per library at p = PEAK_N / pooled_total, so a pool
# keeps its replicate proportions, with a fixed seed for reproducibility.
#
# Control: matched IgG, also capped at PEAK_N, so target and control enter at the
# same depth and SEACR runs in "non" mode -- its own "norm" mode would rescale
# the control to the target and is explicitly not recommended when the data are
# already normalised to each other.
#
# Peaks are called PER CONDITION and UNIONED. Union, not intersection:
# intersecting keeps only regions surviving in both conditions, selecting against
# exactly the regions that change most and biasing every delta toward zero.
#
# RUN DETACHED:
#   nohup ./scripts/phase2_call_domains_seacr.sh > logs/phase2seacr_$(date +%F_%H%M).log 2>&1 &
# =============================================================================
set -euo pipefail

ROOT="${ROOT:-$HOME/Desktop/BCAT1/BCAT1/CUTRUN}"
BAMDIR="${BAMDIR:-$ROOT/hg38}"
QC="${QC:-$ROOT/analysis_QC_normalization/multimark_2026-07}"
OUT="${OUT:-$QC/domains}"
SEACR="${SEACR:-$HOME/tools/SEACR/SEACR_1.3.sh}"
GENOME="${GENOME:-$QC/annotation/genome.txt}"

PEAK_N="${PEAK_N:-5000000}"     # fragments per pool fed to the caller
MODE="${MODE:-relaxed}"         # relaxed: threshold between knee and peak
SEED="${SEED:-11}"
THREADS="${THREADS:-6}"

ENVBIN="${ENVBIN:-$HOME/miniconda3/envs/k4width/bin}"
export PATH="$ENVBIN:$PATH"
[[ -x "$SEACR" ]] || { echo "[fatal] SEACR not at $SEACR" >&2; exit 1; }
command -v bedtools >/dev/null || { echo "[fatal] bedtools missing" >&2; exit 1; }
# SEACR calls Rscript; the conda env's R has no edgeR/etc but SEACR needs only base R
command -v Rscript >/dev/null || { echo "[fatal] Rscript missing" >&2; exit 1; }

mkdir -p "$OUT"/{bg,peaks} "$ROOT/logs"
BLBED="$OUT/hg38-blacklist.v2.bed"
[[ -s "$BLBED" ]] || gunzip -c "$ROOT/hg38-blacklist.v2.bed.gz" > "$BLBED"
MANIFEST="$QC/sample_manifest.tsv"
PROMCAT="$QC/promcat"

GB2UCSC="CM000663.2:chr1 CM000664.2:chr2 CM000665.2:chr3 CM000666.2:chr4 \
CM000667.2:chr5 CM000668.2:chr6 CM000669.2:chr7 CM000670.2:chr8 CM000671.2:chr9 \
CM000672.2:chr10 CM000673.2:chr11 CM000674.2:chr12 CM000675.2:chr13 \
CM000676.2:chr14 CM000677.2:chr15 CM000678.2:chr16 CM000679.2:chr17 \
CM000680.2:chr18 CM000681.2:chr19 CM000682.2:chr20 CM000683.2:chr21 \
CM000684.2:chr22 CM000685.2:chrX"

# cond may be "ALL", meaning both conditions pooled into one call
pool_total () {                  # $1 mark  $2 cond -> summed fragment total
  local t=0 n
  while IFS=$'\t' read -r label mk cd rep naming src; do
    [[ "$mk" == "$1" ]] || continue
    [[ "$2" == "ALL" || "$cd" == "$2" ]] || continue
    n=$(awk -F'\t' '$1=="total"{print $2}' "$PROMCAT/${label}.category.tsv")
    t=$(( t + n ))
  done < <(tail -n +2 "$MANIFEST")
  echo "$t"
}

# ------------------------------------------------ capped, sorted bedgraph -----
make_bg () {                     # $1 mark  $2 cond
  local mark=$1 cond=$2
  local bg="$OUT/bg/${cond}-${mark}.bg"
  [[ -s "$bg" ]] && { echo "[skip] bg $cond-$mark"; return; }
  local tot p
  tot=$(pool_total "$mark" "$cond")
  (( tot > 0 )) || { echo "[fatal] no fragments for $cond-$mark" >&2; return 1; }
  p=$(awk -v n="$PEAK_N" -v t="$tot" 'BEGIN{printf "%.8f", (n<t? n/t : 1)}')
  echo "[bg] $cond-$mark  pool=$tot  keep p=$p  -> target ~$PEAK_N fragments"

  local tmp="$OUT/bg/${cond}-${mark}.frag.tmp"
  : > "$tmp"
  while IFS=$'\t' read -r label mk cd rep naming src; do
    [[ "$mk" == "$mark" ]] || continue
    [[ "$cond" == "ALL" || "$cd" == "$cond" ]] || continue
    local regions map
    if [[ $naming == gb ]]; then
      map="$GB2UCSC"; regions=$(for q in $GB2UCSC; do printf '%s ' "${q%%:*}"; done)
    else
      map=""; regions=$(printf 'chr%s ' {1..22} X)
    fi
    # shellcheck disable=SC2086
    samtools view -@ "$THREADS" -q 30 -f 66 -F 1804 "$BAMDIR/$src" $regions \
      | awk -v MAP="$map" -v P="$p" -v S="$SEED" 'BEGIN{OFS="\t"; srand(S)
             n=split(MAP,kv," "); for(i=1;i<=n;i++){split(kv[i],q,":"); m[q[1]]=q[2]}}
           { if(rand()>P) next
             t=$9; if(t<0) t=-t
             if(t<=0) next
             c=($3 in m)?m[$3]:$3
             s=($4<$8?$4:$8)-1
             print c, s, s+t }' >> "$tmp"
    echo "    + $label"
  done < <(tail -n +2 "$MANIFEST")

  echo "    $(wc -l < "$tmp" | tr -d ' ') fragments -> bedgraph"
  sort -k1,1 -k2,2n -S 2G "$tmp" \
    | bedtools genomecov -bg -i - -g "$GENOME" > "$bg.partial"
  mv "$bg.partial" "$bg"
  rm -f "$tmp"
  echo "    $(wc -l < "$bg" | tr -d ' ') bedgraph intervals"
}

# ONE control: both IgG libraries pooled. Calling per condition against its own
# IgG was structurally unfair -- shBCAT1-IgG fragments are 26% longer, so at equal
# fragment counts it covered 775 Mb against shNT-IgG's 608 Mb (27% more), which in
# "non" mode raises the shBCAT1 threshold and suppressed its peaks in every mark
# (H3K4me3 71,814 -> 19,556; H3K27ac 92,359 -> 331). Same fragment-length defect
# that removed IgG from the differential estimator in PLAN_direct_ratio_migration.
make_bg IgG ALL

summary="$OUT/domain_summary.tsv"
printf 'mark\tbasis\tdomains\ttotal_bp\tmedian_width\n' > "$summary"

# ONE domain set per mark, called on BOTH conditions combined.
#
# Not a union of per-condition calls: those thresholds proved unstable on this
# data (ENL inverted, 1,775 shNT vs 81,991 shBCAT1, on a mark whose enrichment is
# 1.64x against IgG's 1.18x). A single call on the combined signal is symmetric in
# the two conditions by construction, so it cannot bias the contrast, and it does
# not inherit either condition's threshold accident. Delta is still measured from
# full-depth per-library counts in Phase 3; these regions only decide membership.
for mark in H3K4me3 H3K27ac H3K9pr ENL; do
  echo "=== $mark ==="
  pk="$OUT/peaks/ALL-${mark}"
  if [[ ! -s "${pk}_final.bed" ]]; then
    make_bg "$mark" ALL
    echo "[seacr] $mark (combined conditions, $MODE)"
    bash "$SEACR" "$OUT/bg/ALL-${mark}.bg" "$OUT/bg/ALL-IgG.bg" \
         non "$MODE" "$pk" > "$pk.seacr.log" 2>&1 || {
           echo "[fatal] SEACR failed for $mark; see $pk.seacr.log" >&2; exit 1; }
    raw="${pk}.${MODE}.bed"
    [[ -s "$raw" ]] || { echo "[fatal] no SEACR output $raw" >&2; exit 1; }
    awk 'BEGIN{OFS="\t"} $1 ~ /^chr([0-9]+|X)$/ {print $1,$2,$3}' "$raw" \
      | sort -k1,1 -k2,2n | bedtools merge -i - \
      | bedtools subtract -a - -b "$BLBED" > "${pk}_final.bed"
    rm -f "$OUT/bg/ALL-${mark}.bg"
  fi
  n=$(wc -l < "${pk}_final.bed" | tr -d ' ')
  tb=$(awk '{t+=$3-$2} END{printf "%d", t+0}' "${pk}_final.bed")
  mw=$(awk '{print $3-$2}' "${pk}_final.bed" | sort -n \
       | awk -v n="$n" 'NR==int((n+1)/2){v=$1} END{print v+0}')
  printf '%s\tcombined\t%s\t%s\t%s\n' "$mark" "$n" "$tb" "${mw:-0}" >> "$summary"
  awk -v m="$mark" 'BEGIN{OFS="\t"}{print $1,$2,$3,m"_dom_"NR}' "${pk}_final.bed" \
    > "$OUT/${mark}.union.bed"
  echo "  $mark: $n domains, $(awk -v b="$tb" 'BEGIN{printf "%.1f", b/1e6}') Mb, median ${mw:-0} bp"
done

echo
column -t "$summary"
echo "[done] domains in $OUT"
df -h "$ROOT" | awk 'NR==2{print "[disk] free: "$4}'
