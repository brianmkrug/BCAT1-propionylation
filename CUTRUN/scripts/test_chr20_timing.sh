#!/usr/bin/env bash
# =============================================================================
# chr20 timing gate — run BEFORE the genome-wide Phase 1.
#
# Reproduces the real Stage-1 + Stage-4 configuration (header surgery, corrected
# flags) on ONE sample, times bamCoverage on chr20, and extrapolates.
# chr20 ~= 2% of the genome, so multiply by ~50.
#
# Uses ABSOLUTE binary paths — no `conda activate` (per plan §Phase 0: activating
# inside a Claude Code tool call can hang the shell-snapshot step).
#
# Run in a real Terminal:
#   bash ${BCAT1_REPO:-.}/CUTRUN/scripts/test_chr20_timing.sh
# =============================================================================
set -euo pipefail

ENVBIN="${ENVBIN:-$HOME/miniconda3/envs/k4width/bin}"
ROOT="${ROOT:-$HOME/Desktop/BCAT1/BCAT1/CUTRUN}"
BAMDIR="$ROOT/hg38"
WORK="$ROOT/H3K4me3/timing_test"
SAMPLE="${SAMPLE:-DIPG13-shNT-H3K4me3-r1_nodup_sorted_unique}"
THREADS="${THREADS:-8}"
BINSIZE=50

SAMTOOLS="$ENVBIN/samtools"
BAMCOVERAGE="$ENVBIN/bamCoverage"

mkdir -p "$WORK"
BL="$ROOT/H3K4me3/reference/hg38-blacklist.v2.bed"
[[ -s "$BL" ]] || { mkdir -p "$(dirname "$BL")"; gunzip -c "$ROOT/hg38-blacklist.v2.bed.gz" > "$BL"; }

echo "[arch] host=$(uname -m)  translated(1=Rosetta)=$(sysctl -n sysctl.proc_translated 2>/dev/null || echo 0)"
echo "[arch] samtools: $(file -b "$SAMTOOLS" | sed 's/,.*//')"
echo "[cfg]  sample=$SAMPLE  threads=$THREADS  binSize=$BINSIZE"
echo

# --- how bad is the header? (the documented deepTools slowdown) --------------
NSQ=$("$SAMTOOLS" view -H "$BAMDIR/$SAMPLE.bam" | grep -c '^@SQ')
echo "[header] @SQ records in source: $NSQ  (expect ~711; 686 are decoy/alt)"

# --- source index is required for a region-restricted view -------------------
if [[ ! -s "$BAMDIR/$SAMPLE.bam.bai" ]]; then
  echo "[index] building missing source index ..."
  time "$SAMTOOLS" index -@ "$THREADS" "$BAMDIR/$SAMPLE.bam"
fi

# --- Stage-1 equivalent, chr20 only, with header surgery ---------------------
MAIN="$WORK/${SAMPLE}.chr20.bam"
if [[ ! -s "$MAIN.bai" ]]; then
  echo "[stage1] header surgery + chr20 extract ..."
  time { "$SAMTOOLS" view -h "$BAMDIR/$SAMPLE.bam" chr20 \
    | awk 'BEGIN{FS=OFS="\t"} /^@SQ/ { if ($2=="SN:chr20") print; next } { print }' \
    | "$SAMTOOLS" view -b -@ "$THREADS" -o "$MAIN" - ; }
  "$SAMTOOLS" index -@ "$THREADS" "$MAIN"
fi
echo "[header] @SQ after surgery: $("$SAMTOOLS" view -H "$MAIN" | grep -c '^@SQ')"
echo

# --- Stage-4 equivalent: the step whose runtime we care about ----------------
echo "[stage4] bamCoverage on chr20 (this is the timed step) ..."
START=$(date +%s)
"$BAMCOVERAGE" -b "$MAIN" -o "$WORK/chr20.cpm.bw" \
  --binSize "$BINSIZE" --normalizeUsing CPM \
  --extendReads --samFlagInclude 66 --samFlagExclude 1804 \
  --blackListFileName "$BL" -p "$THREADS"
END=$(date +%s)
ELAPSED=$((END-START))

echo
echo "==================== RESULT ===================="
printf "chr20 bamCoverage: %d s\n" "$ELAPSED"
printf "extrapolated genome-wide (x50): %d s  (~%d min)\n" "$((ELAPSED*50))" "$(((ELAPSED*50)/60))"
printf "12 samples at 3 concurrent (~x4 wall): ~%d min\n" "$((((ELAPSED*50)/60)*4))"
echo "==============================================="
echo "Gate: if chr20 alone exceeds a few minutes the bottleneck is I/O or"
echo "emulation, not CPU — confirm BAMs are on local storage before the full run."
