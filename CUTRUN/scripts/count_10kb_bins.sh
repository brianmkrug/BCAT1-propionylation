#!/usr/bin/env bash
# =============================================================================
# Harmonized 10 kb fragment counts for the multi-mark QC pass
#   marks retained : H3K4me3, H3K9pr, ENL, H3K27ac, IgG
#   marks dropped  : RNAPII-Ser2P (Nat + Xlink) -- failed IP enrichment
#                    (1.0-2.0x IgG, see H3K4me3/reports/phase6_findings_report.md)
#
# Why recount instead of reusing analysis_QC_normalization/data/counts_10kb_matrix.npz:
#   the four new H3K27ac BAMs differ from the 2025 set in two ways that must be
#   harmonized BEFORE any cross-sample comparison, and the fix has to be applied
#   identically to every column of the matrix:
#     1. reference naming -- H3K27ac BAMs were aligned to a GenBank-named GRCh38
#        (CM000663.2 = chr1 ... CM000685.2 = chrX, J01415.2 = chrM). Mapping below
#        was derived by contig-length matching against the chr-named BAMs; it is
#        1:1 and unambiguous for all 24 primaries.
#     2. MAPQ -- the two H3K27ac r1 libraries were filtered at `samtools view -bq 1`
#        whereas every other library in the project (including H3K27ac r2) was
#        filtered at `-q 30`. Counting at MAPQ>=30 is a no-op for the pre-filtered
#        libraries and brings the r1 K27ac libraries onto the same footing.
#
# Fragment definition (identical for every sample):
#   -f 66   properly paired AND first-in-pair  -> one fragment per read pair
#   -F 1804 drop unmapped / mate-unmapped / secondary / QC-fail / duplicate
#   -q 30   MAPQ floor (see above)
#   bin assigned from the leftmost position of read1; at 10 kb the ~130 bp
#   read1-vs-fragment-midpoint offset is immaterial.
#
# Grid: the exact bin grid of the 2026-07-08 matrix -- chr1-22 + chrX, 10 kb,
# origin 0, 303,114 bins. chrY/chrM/scaffolds excluded.
#
# Output: one gzipped 3-column TSV per sample (chrom, start, fragments) in
# $OUT/counts/, plus a fragment-total TSV. Restartable: finished samples skip.
#
# RUN DETACHED:
#   cd ${BCAT1_REPO:-.}/CUTRUN
#   nohup ./scripts/count_10kb_bins.sh > logs/count10kb_$(date +%F_%H%M).log 2>&1 &
# =============================================================================
set -euo pipefail

ROOT="${ROOT:-$HOME/Desktop/BCAT1/BCAT1/CUTRUN}"
BAMDIR="${BAMDIR:-$ROOT/hg38}"
OUT="${OUT:-$ROOT/analysis_QC_normalization/multimark_2026-07}"
BIN=10000

ENVBIN="${ENVBIN:-$HOME/miniconda3/envs/k4width/bin}"
[[ -x "$ENVBIN/samtools" ]] || { echo "[fatal] no samtools at $ENVBIN" >&2; exit 1; }
export PATH="$ENVBIN:$PATH"

if PCORES=$(sysctl -n hw.perflevel0.logicalcpu 2>/dev/null); then :; else
  PCORES=$(sysctl -n hw.logicalcpu 2>/dev/null || nproc)
fi
JOBS="${JOBS:-3}"
THREADS="${THREADS:-$(( PCORES / JOBS ))}"; (( THREADS < 2 )) && THREADS=2
echo "[cfg] P-cores=$PCORES jobs=$JOBS threads/job=$THREADS bin=$BIN"

mkdir -p "$OUT/counts" "$ROOT/logs"

# ----------------------------------------------------------------- manifest --
# label|mark|condition|rep|naming|source basename (no .bam)
#   naming = ucsc  -> BAM already uses chr1..chr22,chrX
#            gb    -> BAM uses CM0006xx.2 GenBank accessions (needs translation)
SAMPLES=(
  "shNT-H3K4me3-r1|H3K4me3|shNT|r1|ucsc|DIPG13-shNT-H3K4me3-r1_nodup_sorted_unique"
  "shNT-H3K4me3-r2|H3K4me3|shNT|r2|ucsc|DIPG13-shNT-H3K4me3-r2_nodup_sorted_unique"
  "shBCAT1-H3K4me3-r1|H3K4me3|shBCAT1|r1|ucsc|DIPG13-shBCAT1-H3K4me3-r1_nodup_sorted_unique"
  "shBCAT1-H3K4me3-r2|H3K4me3|shBCAT1|r2|ucsc|DIPG13-shBCAT1-H3K4me3-r2_nodup_sorted_unique"
  "shNT-H3K9pr-r1|H3K9pr|shNT|r1|ucsc|DIPG13-shNT-H3K9pr-r1_nodup_sorted_unique"
  "shNT-H3K9pr-r2|H3K9pr|shNT|r2|ucsc|DIPG13-shNT-H3K9pr-r2_nodup_sorted_unique"
  "shBCAT1-H3K9pr-r1|H3K9pr|shBCAT1|r1|ucsc|DIPG13-shBCAT1-H3K9pr-r1_nodup_sorted_unique"
  "shBCAT1-H3K9pr-r2_60M|H3K9pr|shBCAT1|r2|ucsc|DIPG13-shBCAT1-H3K9pr-r2_60M_nodup_sorted_unique"
  "shNT-ENL-r1|ENL|shNT|r1|ucsc|DIPG13-shNT-ENL-r1_nodup_sorted_unique"
  "shNT-ENL-r2|ENL|shNT|r2|ucsc|DIPG13-shNT-ENL-r2_nodup_sorted_unique"
  "shBCAT1-ENL-r1|ENL|shBCAT1|r1|ucsc|DIPG13-shBCAT1-ENL-r1_nodup_sorted_unique"
  "shBCAT1-ENL-r2|ENL|shBCAT1|r2|ucsc|DIPG13-shBCAT1-ENL-r2_nodup_sorted_unique"
  "shNT-H3K27ac-r1|H3K27ac|shNT|r1|gb|shNT-H3K27ac-r1_nodup_sorted_unique"
  "shNT-H3K27ac-r2|H3K27ac|shNT|r2|gb|shNT-H3K27ac-r2_nodup_sorted_unique"
  "shBCAT1-H3K27ac-r1|H3K27ac|shBCAT1|r1|gb|shBCAT1-H3K27ac-r1_nodup_sorted_unique"
  "shBCAT1-H3K27ac-r2|H3K27ac|shBCAT1|r2|gb|shBCAT1-H3K27ac-r2_nodup_sorted_unique"
  "shNT-IgG-r1|IgG|shNT|r1|ucsc|DIPG13-shNT-IgG-r1_nodup_sorted_unique"
  "shBCAT1-IgG-r1|IgG|shBCAT1|r1|ucsc|DIPG13-shBCAT1-IgG-r1_nodup_sorted_unique"
)

# --------------------------------------------------- chromosome name mapping --
# GenBank accession -> UCSC, verified by contig length against the chr-named BAMs.
GB2UCSC="CM000663.2:chr1 CM000664.2:chr2 CM000665.2:chr3 CM000666.2:chr4 \
CM000667.2:chr5 CM000668.2:chr6 CM000669.2:chr7 CM000670.2:chr8 CM000671.2:chr9 \
CM000672.2:chr10 CM000673.2:chr11 CM000674.2:chr12 CM000675.2:chr13 \
CM000676.2:chr14 CM000677.2:chr15 CM000678.2:chr16 CM000679.2:chr17 \
CM000680.2:chr18 CM000681.2:chr19 CM000682.2:chr20 CM000683.2:chr21 \
CM000684.2:chr22 CM000685.2:chrX"

CHRS_UCSC=$(printf 'chr%s ' {1..22} X)
CHRS_GB=$(for p in $GB2UCSC; do printf '%s ' "${p%%:*}"; done)

# A region-restricted `samtools view` requires an index; several of the 2025 BAMs
# (ENL, H3K9pr) shipped without one. Build them up front, serially, so the
# parallel counting stage below never races two indexers onto the same file.
for s in "${SAMPLES[@]}"; do
  IFS='|' read -r _ _ _ _ _ src <<< "$s"
  bam="$BAMDIR/${src}.bam"
  [[ -s "$bam" ]] || { echo "[fatal] missing $bam" >&2; exit 1; }
  [[ -s "$bam.bai" ]] || { echo "[index] $src"; samtools index -@ "$PCORES" "$bam"; }
done

count_one() {                     # $1 label  $2 naming  $3 source basename
  local label=$1 naming=$2 src=$3
  local out="$OUT/counts/${label}.10kb.tsv.gz"
  [[ -s "$out" ]] && { echo "[skip] $label"; return; }
  local bam="$BAMDIR/${src}.bam"
  [[ -s "$bam" ]]      || { echo "[fatal] missing $bam" >&2; return 1; }
  [[ -s "$bam.bai" ]]  || { echo "[fatal] missing index $bam.bai" >&2; return 1; }

  local regions map
  if [[ $naming == gb ]]; then regions="$CHRS_GB"; map="$GB2UCSC"; else regions="$CHRS_UCSC"; map=""; fi

  echo "[count] $label ($naming)"
  # shellcheck disable=SC2086
  samtools view -@ "$THREADS" -q 30 -f 66 -F 1804 "$bam" $regions \
    | awk -v BIN="$BIN" -v MAP="$map" 'BEGIN{
          n=split(MAP,kv," ");
          for(i=1;i<=n;i++){ split(kv[i],p,":"); m[p[1]]=p[2] }
        }
        { c = (($3 in m) ? m[$3] : $3);
          b = int(($4-1)/BIN)*BIN;          # SAM POS is 1-based
          k[c"\t"b]++ }
        END{ for(x in k) print x"\t"k[x] }' \
    | sort -k1,1 -k2,2n \
    | gzip -c > "$out.partial"
  mv "$out.partial" "$out"
  echo "[done] $label"
}
export -f count_one
export OUT BAMDIR BIN THREADS CHRS_UCSC CHRS_GB GB2UCSC

printf '%s\n' "${SAMPLES[@]}" \
  | awk -F'|' '{print $1"\t"$5"\t"$6}' \
  | xargs -P "$JOBS" -n3 bash -c 'count_one "$0" "$1" "$2"'

# ------------------------------------------------------------ manifest dump --
MAN="$OUT/sample_manifest.tsv"
{ echo -e "label\tmark\tcondition\trep\tnaming\tsource_bam"
  for s in "${SAMPLES[@]}"; do
    IFS='|' read -r l m c r n src <<< "$s"; echo -e "$l\t$m\t$c\t$r\t$n\t${src}.bam"
  done
} > "$MAN"

echo "[ok] per-sample counts in $OUT/counts ; manifest $MAN"
