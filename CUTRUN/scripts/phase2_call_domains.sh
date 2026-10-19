#!/usr/bin/env bash
# =============================================================================
# PHASE 2 — define the interpretation domain for each mark.
#
# Why this exists: the direct shBCAT1/shNT ratio replaced IgG for measuring
# CHANGE, but a ratio cannot say whether a region has any signal at all --
# delta ~ 0 at a silent locus and delta ~ 0 at a strong shared peak are
# indistinguishable. So IgG keeps its OTHER job: a peak-calling control.
# No IgG value enters any reported magnitude; it only decides membership,
# where the 26% shBCAT1-vs-shNT IgG fragment-length skew shifts a threshold
# negligibly instead of biasing an effect size.
#
# Peaks are called PER CONDITION against the matched IgG, then UNIONED across
# conditions. Union, not intersection: intersecting would keep only regions that
# survived in both conditions, which selects against exactly the regions that
# change most and would bias every downstream delta toward zero.
#
# Input is pooled fragments (both replicates, full depth) as 3-column BED, which
# MACS2 reads with `-f BEDPE`. Fragment coordinates come straight from the pairs
# (min(POS,PNEXT)-1 .. +|TLEN|), the same definition used by every count in this
# project. Building fragment BEDs avoids regenerating the ~35 GB of filtered BAMs
# that were deleted after the track build.
#
# Peak style per mark, matched to its signal character (see QC report 6d):
#   H3K4me3, H3K27ac  -> narrow   (punctate; 6.3x and 3.5x promoter enrichment)
#   H3K9pr, ENL       -> broad    (dispersed; 4.1x and 1.6x)
#
# RUN DETACHED:
#   nohup ./scripts/phase2_call_domains.sh > logs/phase2dom_$(date +%F_%H%M).log 2>&1 &
# =============================================================================
set -euo pipefail

ROOT="${ROOT:-$HOME/Desktop/BCAT1/BCAT1/CUTRUN}"
BAMDIR="${BAMDIR:-$ROOT/hg38}"
QC="${QC:-$ROOT/analysis_QC_normalization/multimark_2026-07}"
OUT="${OUT:-$QC/domains}"

ENVBIN="${ENVBIN:-$HOME/miniconda3/envs/k4width/bin}"
export PATH="$ENVBIN:$PATH"
command -v macs2 >/dev/null || { echo "[fatal] macs2 not found" >&2; exit 1; }

THREADS="${THREADS:-6}"
QNARROW="${QNARROW:-0.01}"
QBROAD="${QBROAD:-0.05}"

mkdir -p "$OUT"/{frag,peaks} "$ROOT/logs"
BLBED="$OUT/hg38-blacklist.v2.bed"
[[ -s "$BLBED" ]] || gunzip -c "$ROOT/hg38-blacklist.v2.bed.gz" > "$BLBED"

MANIFEST="$QC/sample_manifest.tsv"
GB2UCSC="CM000663.2:chr1 CM000664.2:chr2 CM000665.2:chr3 CM000666.2:chr4 \
CM000667.2:chr5 CM000668.2:chr6 CM000669.2:chr7 CM000670.2:chr8 CM000671.2:chr9 \
CM000672.2:chr10 CM000673.2:chr11 CM000674.2:chr12 CM000675.2:chr13 \
CM000676.2:chr14 CM000677.2:chr15 CM000678.2:chr16 CM000679.2:chr17 \
CM000680.2:chr18 CM000681.2:chr19 CM000682.2:chr20 CM000683.2:chr21 \
CM000684.2:chr22 CM000685.2:chrX"

# --------------------------------------------------- pooled fragment BED ------
make_frag () {                    # $1 mark  $2 condition -> $OUT/frag/<cond>-<mark>.bed
  local mark=$1 cond=$2
  local out="$OUT/frag/${cond}-${mark}.bed"
  [[ -s "$out" ]] && { echo "[skip] frag $cond-$mark"; return; }
  echo "[frag] $cond-$mark"
  : > "$out.partial"
  while IFS=$'\t' read -r label mk cd rep naming src; do
    [[ "$mk" == "$mark" && "$cd" == "$cond" ]] || continue
    local regions map
    if [[ $naming == gb ]]; then
      map="$GB2UCSC"; regions=$(for p in $GB2UCSC; do printf '%s ' "${p%%:*}"; done)
    else
      map=""; regions=$(printf 'chr%s ' {1..22} X)
    fi
    # shellcheck disable=SC2086
    samtools view -@ "$THREADS" -q 30 -f 66 -F 1804 "$BAMDIR/$src" $regions \
      | awk -v MAP="$map" 'BEGIN{OFS="\t"; n=split(MAP,kv," ")
             for(i=1;i<=n;i++){split(kv[i],p,":"); m[p[1]]=p[2]}}
           { c=($3 in m)?m[$3]:$3
             t=$9; if(t<0) t=-t
             if(t<=0) next
             s=($4<$8?$4:$8)-1
             print c, s, s+t }' >> "$out.partial"
    echo "    + $label"
  done < <(tail -n +2 "$MANIFEST")
  mv "$out.partial" "$out"
  echo "    $(wc -l < "$out" | tr -d ' ') fragments"
}

# IgG controls are reused by every mark -> build once, keep
make_frag IgG shNT
make_frag IgG shBCAT1

summary="$OUT/domain_summary.tsv"
[[ -s "$summary" ]] || printf 'mark\tstyle\tcondition\tpeaks\ttotal_bp\tmedian_width\n' > "$summary"

for mark in H3K4me3 H3K27ac H3K9pr ENL; do
  case "$mark" in
    H3K4me3|H3K27ac) style=narrow; EXTRA="-q $QNARROW" ;;
    *)               style=broad;  EXTRA="--broad --broad-cutoff $QBROAD" ;;
  esac
  echo "=== $mark ($style) ==="
  for cond in shNT shBCAT1; do
    pk="$OUT/peaks/${cond}-${mark}"
    if [[ ! -s "${pk}_final.bed" ]]; then
      make_frag "$mark" "$cond"
      echo "[macs2] $cond-$mark ($style)"
      # shellcheck disable=SC2086
      macs2 callpeak -t "$OUT/frag/${cond}-${mark}.bed" \
                     -c "$OUT/frag/${cond}-IgG.bed" \
                     -f BEDPE -g hs -n "${cond}-${mark}" --outdir "$OUT/peaks" \
                     --keep-dup all $EXTRA 2> "$OUT/peaks/${cond}-${mark}.macs2.log"
      raw="$OUT/peaks/${cond}-${mark}_peaks.narrowPeak"
      [[ $style == broad ]] && raw="$OUT/peaks/${cond}-${mark}_peaks.broadPeak"
      # drop blacklist and non-primary contigs, merge touching calls
      awk 'BEGIN{OFS="\t"} $1 ~ /^chr([0-9]+|X)$/ {print $1,$2,$3}' "$raw" \
        | sort -k1,1 -k2,2n | bedtools merge -i - \
        | bedtools subtract -a - -b "$BLBED" > "${pk}_final.bed"
    fi
    n=$(wc -l < "${pk}_final.bed" | tr -d ' ')
    # portable width stats: no gawk asort, no in-memory sort in awk
    total_bp=$(awk '{t+=$3-$2} END{printf "%d", t+0}' "${pk}_final.bed")
    median_w=$(awk '{print $3-$2}' "${pk}_final.bed" | sort -n \
               | awk -v n="$n" 'NR==int((n+1)/2){print; exit}')
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$mark" "$style" "$cond" "$n" \
           "$total_bp" "${median_w:-0}" >> "$summary"
    echo "  $cond: $n peaks, $(awk -v b="$total_bp" 'BEGIN{printf "%.1f", b/1e6}') Mb, median ${median_w:-0} bp"
  done

  # union across conditions = the interpretation domain
  cat "$OUT/peaks/shNT-${mark}_final.bed" "$OUT/peaks/shBCAT1-${mark}_final.bed" \
    | sort -k1,1 -k2,2n | bedtools merge -i - \
    | awk -v m="$mark" 'BEGIN{OFS="\t"}{print $1,$2,$3,m"_dom_"NR}' > "$OUT/${mark}.union.bed"
  echo "  UNION: $(wc -l < "$OUT/${mark}.union.bed" | tr -d ' ') domains"

  # target fragment BEDs are large; drop once peaks exist (IgG kept)
  rm -f "$OUT/frag/shNT-${mark}.bed" "$OUT/frag/shBCAT1-${mark}.bed"
done

echo
column -t "$summary"
echo "[done] domains in $OUT"
df -h "$ROOT" | awk 'NR==2{print "[disk] free: "$4}'
