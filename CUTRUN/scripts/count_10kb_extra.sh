#!/usr/bin/env bash
# =============================================================================
# 10 kb counts for samples added after the main pass (CTCF, H3K4me1).
#
# Deliberately separate from count_10kb_bins.sh rather than an edit to it:
# that script REWRITES sample_manifest.tsv at the end, and the SEACR domain job
# reads that manifest while running. This one only writes into counts/ and its
# own manifest.
#
# The binning rule is copied EXACTLY from count_10kb_bins.sh -- leftmost position
# of read1, `int((POS-1)/BIN)*BIN` -- not the fragment midpoint used by
# count_regions.py. The existing 18-column matrix was built with the read1 rule,
# and a matrix whose columns disagree about fragment placement is worse than one
# built on a slightly cruder rule consistently.
#
# Both new BAMs are already UCSC-named and MAPQ30-filtered, so unlike the H3K27ac
# batch they need no renaming and -q 30 is a no-op re-application.
# =============================================================================
set -euo pipefail

ROOT="${ROOT:-$HOME/Desktop/BCAT1/BCAT1/CUTRUN}"
BAMDIR="${BAMDIR:-$ROOT/hg38}"
OUT="${OUT:-$ROOT/analysis_QC_normalization/multimark_2026-07}"
BIN=10000

ENVBIN="${ENVBIN:-$HOME/miniconda3/envs/k4width/bin}"
export PATH="$ENVBIN:$PATH"
THREADS="${THREADS:-6}"

mkdir -p "$OUT/counts"

# label|mark|condition|rep|naming|source basename
SAMPLES=(
  "CTCF|CTCF|shNT|r1|ucsc|CTCF_nodup_sorted_unique"
  "H3K4me1|H3K4me1|shNT|r1|ucsc|K4me1_nodup_sorted_unique"
)
CHRS=$(printf 'chr%s ' {1..22} X)

for s in "${SAMPLES[@]}"; do
  IFS='|' read -r label mark cond rep naming src <<< "$s"
  bam="$BAMDIR/${src}.bam"
  [[ -s "$bam" ]] || { echo "[fatal] missing $bam" >&2; exit 1; }
  [[ -s "$bam.bai" ]] || { echo "[index] $label"; samtools index -@ "$THREADS" "$bam"; }

  out="$OUT/counts/${label}.10kb.tsv.gz"
  if [[ -s "$out" ]]; then echo "[skip] $label"; continue; fi
  echo "[count] $label"
  # shellcheck disable=SC2086
  samtools view -@ "$THREADS" -q 30 -f 66 -F 1804 "$bam" $CHRS \
    | awk -v BIN="$BIN" '{ b = int(($4-1)/BIN)*BIN; k[$3"\t"b]++ }
          END{ for(x in k) print x"\t"k[x] }' \
    | sort -k1,1 -k2,2n | gzip -c > "$out.partial"
  mv "$out.partial" "$out"
  echo "[done] $label  $(gzip -dc "$out" | awk '{t+=$3} END{print t}') fragments"
done

MAN="$OUT/sample_manifest_extra.tsv"
{ echo -e "label\tmark\tcondition\trep\tnaming\tsource_bam\tassay\tnote"
  echo -e "CTCF\tCTCF\tshNT\tr1\tucsc\tCTCF_nodup_sorted_unique.bam\tCUT&RUN\tshNT DIPG-XIII; Diagenode C15410210; single library; no matched IgG"
  echo -e "H3K4me1\tH3K4me1\tshNT\tr1\tucsc\tK4me1_nodup_sorted_unique.bam\tCUT&RUN\tshNT DIPG-XIII; CST 5326; single library; no matched IgG"
} > "$MAN"
echo "[ok] counts in $OUT/counts ; manifest $MAN"
