#!/usr/bin/env bash
# =============================================================================
# Phase 1 / 1b — chromosome-filtered, depth-matched CPM tracks
# DIPGXIII shNT vs shBCAT1 CUT&RUN : H3K4me3, RNAPII-Ser2P (Nat + Xlink), IgG
#
# Replaces the slow single-threaded `bamCompare` approach with:
#   filter -> count -> downsample -> bamCoverage (CPM) -> bigwigCompare
#
# Designed for Apple Silicon (M4 Pro, 24 GB). Restartable: every stage skips
# work that is already done, so you can Ctrl-C and re-run safely.
#
# RUN DETACHED, NOT INSIDE CLAUDE CODE:
#   cd ${BCAT1_REPO:-.}/CUTRUN
#   chmod +x scripts/run_phase1_tracks.sh
#   nohup ./scripts/run_phase1_tracks.sh > logs/phase1_$(date +%F_%H%M).log 2>&1 &
#   tail -f logs/phase1_*.log
# =============================================================================

set -euo pipefail

# ---------------------------------------------------------------- config ----
ROOT="${ROOT:-$HOME/Desktop/BCAT1/BCAT1/CUTRUN}"
BAMDIR="${BAMDIR:-$ROOT/hg38}"
OUT="${OUT:-$ROOT/H3K4me3/phase1_tracks}"   # analysis outputs live under CUTRUN/H3K4me3/
BLACKLIST_GZ="${BLACKLIST_GZ:-$ROOT/hg38-blacklist.v2.bed.gz}"

BINSIZE="${BINSIZE:-50}"
SEED="${SEED:-42}"                 # fixed downsampling seed — do not change mid-project

# Put the analysis env on PATH here, rather than relying on the caller's shell.
# The plan warns against `conda activate` inside a Claude Code tool call (it can
# hang the shell-snapshot step), and a bare `nohup ./run_phase1_tracks.sh` from a
# terminal without the env active would otherwise die at stage 1 with
# "samtools: command not found".
ENVBIN="${ENVBIN:-$HOME/miniconda3/envs/k4width/bin}"
[[ -x "$ENVBIN/samtools" ]] || { echo "[fatal] no samtools at $ENVBIN — set ENVBIN=" >&2; exit 1; }
export PATH="$ENVBIN:$PATH"

# Fragment definition, matched to the QC report (one fragment per read pair):
#   -f 66   = properly paired (2) AND first in pair (64)
#   -F 1804 = drop unmapped, mate-unmapped, secondary, QC-fail, duplicate
SAM_INCLUDE=66
SAM_EXCLUDE=1804

# ------------------------------------------------------------ core budget ----
# On Apple Silicon prefer performance cores; E-cores add little and cause
# scheduler thrash. M4 Pro = 8 or 10 P-cores + 4 E-cores.
if PCORES=$(sysctl -n hw.perflevel0.logicalcpu 2>/dev/null); then :; else
  PCORES=$(sysctl -n hw.logicalcpu 2>/dev/null || nproc)
fi
JOBS="${JOBS:-3}"                              # concurrent samples
THREADS="${THREADS:-$(( PCORES / JOBS ))}"
(( THREADS < 2 )) && THREADS=2
echo "[cfg] P-cores=$PCORES  concurrent jobs=$JOBS  threads/job=$THREADS  (total ~$((JOBS*THREADS)))"
echo "[cfg] deepTools parallelism scales sub-linearly; 3 jobs x N threads beats 1 job x 3N."

# ------------------------------------------------- provenance: architecture ----
# Recorded, not enforced. Rosetta costs roughly 15-30% on these tools — worth
# knowing when you interpret timings, not worth blocking a run over.
echo "[arch] host=$(uname -m)  python=$(python -c 'import platform;print(platform.machine())' 2>/dev/null || echo '?')"
echo "[arch] this shell translated by Rosetta: $(sysctl -n sysctl.proc_translated 2>/dev/null || echo 0) (1=yes)"
echo "[arch] samtools binary: $(file -b "$(command -v samtools)" 2>/dev/null | sed 's/,.*//')"
echo "[arch] conda subdir(s) in env: $(conda list --explicit 2>/dev/null | grep -oE 'osx-(64|arm64)' | sort -u | tr '\n' ' ')"

mkdir -p "$OUT"/{01_main,02_ds,03_bw,04_cmp,qc} "$ROOT/logs"

BLACKLIST="$OUT/hg38-blacklist.v2.bed"
[[ -s "$BLACKLIST" ]] || gunzip -c "$BLACKLIST_GZ" > "$BLACKLIST"

CHRS=$(printf 'chr%s ' {1..22} X)              # chrY/chrM excluded, per QC report

# ------------------------------------------------------------- manifest ----
# label|prep|condition|mark|source basename (without .bam)
SAMPLES=(
  "shNT-H3K4me3-r1|na|shNT|H3K4me3|DIPG13-shNT-H3K4me3-r1_nodup_sorted_unique"
  "shNT-H3K4me3-r2|na|shNT|H3K4me3|DIPG13-shNT-H3K4me3-r2_nodup_sorted_unique"
  "shBCAT1-H3K4me3-r1|na|shBCAT1|H3K4me3|DIPG13-shBCAT1-H3K4me3-r1_nodup_sorted_unique"
  "shBCAT1-H3K4me3-r2|na|shBCAT1|H3K4me3|DIPG13-shBCAT1-H3K4me3-r2_nodup_sorted_unique"
  "shNT-Pol2-Nat-r1|nat|shNT|Pol2|DIPG13-shNT-RNAp2s2-Nat-r1_nodup_sorted_unique"
  "shNT-Pol2-Nat-r2|nat|shNT|Pol2|DIPG13-shNT-RNAp2s2-Nat-r2_nodup_sorted_unique"
  "shBCAT1-Pol2-Nat-r1|nat|shBCAT1|Pol2|DIPG13-shBCAT1-RNAp2s2-Nat-r1_nodup_sorted_unique"
  "shBCAT1-Pol2-Nat-r2|nat|shBCAT1|Pol2|DIPG13-shBCAT1-RNAp2s2-Nat-r2_nodup_sorted_unique"
  "shNT-Pol2-Xlink-r1|xlink|shNT|Pol2|DIPG13-shNT-RNAp2s2-Xlink-r1_nodup_sorted_unique"
  "shBCAT1-Pol2-Xlink-r1|xlink|shBCAT1|Pol2|DIPG13-shBCAT1-RNAp2s2-Xlink-r1_nodup_sorted_unique"
  "shNT-IgG-r1|na|shNT|IgG|DIPG13-shNT-IgG-r1_nodup_sorted_unique"
  "shBCAT1-IgG-r1|na|shBCAT1|IgG|DIPG13-shBCAT1-IgG-r1_nodup_sorted_unique"
)

# =============================================================================
# STAGE 1 — restrict to chr1-22,X AND rebuild the header
# -----------------------------------------------------------------------------
# Your BAMs carry 711 @SQ lines (686 _random/alt/decoy contigs). deepTools
# chunks its work off the header, so region-filtering the READS alone does not
# help — the contigs must leave the header too.
#
# NOTE: `samtools reheader` does NOT remap reference IDs, so dropping @SQ lines
# with it would silently corrupt the BAM. The SAM round-trip below is safe:
# samtools rebuilds tid from the header on SAM->BAM conversion.
# =============================================================================
stage1() {
  local label=$1 src=$2
  local out="$OUT/01_main/${label}.main.bam"
  [[ -s "$out.bai" ]] && { echo "[skip] stage1 $label"; return; }
  echo "[stage1] $label"
  # A region-restricted `samtools view` REQUIRES an index on the source BAM.
  # Most of these BAMs shipped without one, which would otherwise abort the run
  # (set -e + pipefail) on the first unindexed sample.
  [[ -s "$BAMDIR/${src}.bam.bai" ]] || {
    echo "[stage1] indexing source $src (missing .bai)"
    samtools index -@ "$THREADS" "$BAMDIR/${src}.bam"
  }
  samtools view -h "$BAMDIR/${src}.bam" $CHRS \
    | awk 'BEGIN{FS=OFS="\t"}
           /^@SQ/ { if ($2 ~ /^SN:chr([0-9]+|X)$/) print; next }
           { print }' \
    | samtools view -b -@ "$THREADS" -o "$out" -
  samtools index -@ "$THREADS" "$out"
}
export -f stage1; export OUT BAMDIR CHRS THREADS

printf '%s\n' "${SAMPLES[@]}" \
  | awk -F'|' '{print $1"\t"$5}' \
  | xargs -P "$JOBS" -n2 bash -c 'stage1 "$0" "$1"'

# =============================================================================
# STAGE 2 — count fragments on the analysed chromosomes
# Depth matching is computed on chr1-22,X only, so it reflects the data that
# actually enters the analysis.
# =============================================================================
COUNTS="$OUT/qc/fragment_counts.tsv"
if [[ ! -s "$COUNTS" ]]; then
  # Write to a temp and promote only on success. An interrupted stage 2 would
  # otherwise leave a PARTIAL non-empty file that this guard then skips, so
  # get_n() would return "" for the missing labels and every downsample target
  # would be computed from garbage.
  tmp="$COUNTS.partial"
  echo -e "label\tprep\tcondition\tmark\tfragments" > "$tmp"
  for s in "${SAMPLES[@]}"; do
    IFS='|' read -r label prep cond mark src <<< "$s"
    n=$(samtools view -c -@ "$THREADS" -f $SAM_INCLUDE -F $SAM_EXCLUDE "$OUT/01_main/${label}.main.bam")
    echo -e "${label}\t${prep}\t${cond}\t${mark}\t${n}" >> "$tmp"
    echo "[stage2] $label -> $n fragments"
  done
  mv "$tmp" "$COUNTS"
fi
column -t "$COUNTS"

get_n() { awk -F'\t' -v l="$1" '$1==l{print $5}' "$COUNTS"; }

# =============================================================================
# STAGE 3 — downsample to matched depth
# -----------------------------------------------------------------------------
# Groups (see plan §Phase 1 / 1b):
#   H3K4me3        : all four to the group minimum (shBCAT1-r2 is ~1.9x deeper)
#   IgG            : both to the group minimum (shNT-IgG is ~2x deeper)
#   Pol2 Nat TierA : all four to the group minimum  -> pooled analyses
#   Pol2 Nat TierB : matched WITHIN replicate pair  -> differential, replicate-blocked
#   Pol2 Xlink     : both to the group minimum
# =============================================================================
min_of() {
  # Hardened: `n` was global, and an empty count (label absent from the counts
  # file) evaluates to 0 in [[ ]] arithmetic, silently making the minimum 0 and
  # the downsample fraction 0 — i.e. an empty BAM. Fail loudly instead.
  local m="" l n
  for l in "$@"; do
    n=$(get_n "$l")
    [[ -n "$n" ]] || { echo "[fatal] no fragment count for '$l' in $COUNTS" >&2; return 1; }
    if [[ -z "$m" ]] || (( n < m )); then m=$n; fi
  done
  echo "$m"
}

downsample() {           # $1 label  $2 target  $3 out-tag
  local label=$1 target=$2 tag=$3
  local in="$OUT/01_main/${label}.main.bam"
  local out="$OUT/02_ds/${label}.${tag}.bam"
  [[ -s "$out.bai" ]] && { echo "[skip] ds $label.$tag"; return; }
  local cur; cur=$(get_n "$label")
  local frac; frac=$(awk -v t="$target" -v c="$cur" 'BEGIN{printf "%.6f", (t<c? t/c : 1)}')
  echo "[stage3] $label.$tag  $cur -> $target  (frac=$frac)"
  if awk -v f="$frac" 'BEGIN{exit !(f>=0.999)}'; then
    ln -sf "$in" "$out"; ln -sf "$in.bai" "$out.bai"
  else
    samtools view -b -@ "$THREADS" -s "${SEED}.${frac#0.}" -o "$out" "$in"
    samtools index -@ "$THREADS" "$out"
  fi
}

K4=(shNT-H3K4me3-r1 shNT-H3K4me3-r2 shBCAT1-H3K4me3-r1 shBCAT1-H3K4me3-r2)
IGG=(shNT-IgG-r1 shBCAT1-IgG-r1)
NAT=(shNT-Pol2-Nat-r1 shNT-Pol2-Nat-r2 shBCAT1-Pol2-Nat-r1 shBCAT1-Pol2-Nat-r2)
XL=(shNT-Pol2-Xlink-r1 shBCAT1-Pol2-Xlink-r1)

K4_MIN=$(min_of "${K4[@]}");  IGG_MIN=$(min_of "${IGG[@]}")
NAT_MIN=$(min_of "${NAT[@]}"); XL_MIN=$(min_of "${XL[@]}")
R1_MIN=$(min_of shNT-Pol2-Nat-r1 shBCAT1-Pol2-Nat-r1)
R2_MIN=$(min_of shNT-Pol2-Nat-r2 shBCAT1-Pol2-Nat-r2)

echo "[targets] H3K4me3=$K4_MIN  IgG=$IGG_MIN  Pol2Nat(TierA)=$NAT_MIN  Pol2Nat r1=$R1_MIN r2=$R2_MIN  Xlink=$XL_MIN"

for l in "${K4[@]}";  do downsample "$l" "$K4_MIN"  ds; done
for l in "${IGG[@]}"; do downsample "$l" "$IGG_MIN" ds; done
for l in "${XL[@]}";  do downsample "$l" "$XL_MIN"  ds; done
for l in "${NAT[@]}"; do downsample "$l" "$NAT_MIN" tierA; done
downsample shNT-Pol2-Nat-r1    "$R1_MIN" tierB
downsample shBCAT1-Pol2-Nat-r1 "$R1_MIN" tierB
downsample shNT-Pol2-Nat-r2    "$R2_MIN" tierB
downsample shBCAT1-Pol2-Nat-r2 "$R2_MIN" tierB

# =============================================================================
# STAGE 4 — bamCoverage -> per-sample CPM bigWigs
# -----------------------------------------------------------------------------
# Corrected flags vs the earlier bamCompare attempt:
#   -p                      was defaulting to 1 (single-threaded)
#   --normalizeUsing CPM    valid here; bamCoverage has no scaleFactorsMethod
#                           conflict (in bamCompare it requires
#                           --scaleFactorsMethod None)
#   --extendReads           use whole fragments
#   --samFlagInclude 66     one fragment per PAIR — deepTools otherwise counts
#                           each mate independently, double-counting fragments
#   --minMappingQuality     omitted: `samtools view -q 30` was already applied
#                           at alignment (see @PG in the BAM header)
# =============================================================================
cover() {                # $1 bam-stem in 02_ds  -> 03_bw/<stem>.cpm.bw
  local stem=$1
  local bam="$OUT/02_ds/${stem}.bam"
  local bw="$OUT/03_bw/${stem}.cpm.bw"
  [[ -s "$bw" ]] && { echo "[skip] cover $stem"; return; }
  echo "[stage4] $stem"
  bamCoverage -b "$bam" -o "$bw" \
    --binSize "$BINSIZE" --normalizeUsing CPM \
    --extendReads --samFlagInclude 66 --samFlagExclude 1804 \
    --blackListFileName "$BLACKLIST" \
    -p "$THREADS"
}
export -f cover; export BINSIZE BLACKLIST

ls "$OUT"/02_ds/*.bam | xargs -n1 basename | sed 's/\.bam$//' \
  | xargs -P "$JOBS" -I{} bash -c 'cover "$@"' _ {}

# =============================================================================
# STAGE 5 — bigwigCompare
# -----------------------------------------------------------------------------
# (i)  IP - IgG            : CPM-subtracted signal, for defining where signal is
# (ii) shBCAT1 vs shNT     : direct condition ratio. Copy number is shared
#                            between conditions and cancels here, which also
#                            sidesteps the missing crosslinked IgG (plan §4.4).
# =============================================================================
cmp2() {                 # $1 bw1  $2 bw2  $3 op  $4 out
  local out="$OUT/04_cmp/$4"
  [[ -s "$out" ]] && { echo "[skip] cmp $4"; return; }
  echo "[stage5] $4"
  bigwigCompare -b1 "$OUT/03_bw/$1" -b2 "$OUT/03_bw/$2" \
    --operation "$3" --binSize "$BINSIZE" --pseudocount 1 \
    --skipZeroOverZero --blackListFileName "$BLACKLIST" \
    -p "$THREADS" -o "$out"
}

# (i) IP minus matched IgG
for l in "${K4[@]}"; do
  c=shNT; [[ $l == shBCAT1* ]] && c=shBCAT1
  cmp2 "${l}.ds.cpm.bw" "${c}-IgG-r1.ds.cpm.bw" subtract "${l}.minusIgG.bw"
done
for l in "${NAT[@]}"; do
  c=shNT; [[ $l == shBCAT1* ]] && c=shBCAT1
  cmp2 "${l}.tierA.cpm.bw" "${c}-IgG-r1.ds.cpm.bw" subtract "${l}.tierA.minusIgG.bw"
done
for l in "${XL[@]}"; do
  c=shNT; [[ $l == shBCAT1* ]] && c=shBCAT1
  cmp2 "${l}.ds.cpm.bw" "${c}-IgG-r1.ds.cpm.bw" subtract "${l}.minusIgG.bw"
done

# (ii) direct shBCAT1 vs shNT, within mark/prep/replicate
cmp2 shBCAT1-H3K4me3-r1.ds.cpm.bw shNT-H3K4me3-r1.ds.cpm.bw log2 delta_H3K4me3_r1.log2.bw
cmp2 shBCAT1-H3K4me3-r2.ds.cpm.bw shNT-H3K4me3-r2.ds.cpm.bw log2 delta_H3K4me3_r2.log2.bw
cmp2 shBCAT1-Pol2-Nat-r1.tierB.cpm.bw shNT-Pol2-Nat-r1.tierB.cpm.bw log2 delta_Pol2_nat_r1.log2.bw
cmp2 shBCAT1-Pol2-Nat-r2.tierB.cpm.bw shNT-Pol2-Nat-r2.tierB.cpm.bw log2 delta_Pol2_nat_r2.log2.bw
cmp2 shBCAT1-Pol2-Xlink-r1.ds.cpm.bw shNT-Pol2-Xlink-r1.ds.cpm.bw log2 delta_Pol2_xlink_r1.log2.bw

echo
echo "[done] tracks in $OUT/03_bw and $OUT/04_cmp"
echo "[next] QC: plotFingerprint / multiBigwigSummary on 03_bw, then Phase 2."
