#!/usr/bin/env bash
# =============================================================================
# Disjoint genomic partition for the coverage-category breakdown.
#
# Four categories, assigned by priority so every base belongs to exactly one:
#   1 promoter    TSS +/- 1 kb   (H3K4me3/reference/core.bed, the same windows the
#                                 rest of the project uses -- 36,689 promoters,
#                                 merged here because many overlap)
#   2 exonic      GENCODE v38 exon union, minus promoter
#   3 intronic    gene-body union, minus exons, minus promoter
#   4 intergenic  everything else on chr1-22 + chrX
#
# Priority matters: a promoter base that is also exonic counts as promoter, so
# "is H3K27ac promoter-focused" is not answered by double-counting.
#
# Output: partition.bed  (chrom, start, end, category_code) sorted, disjoint,
# covering chr1-22,X exactly once.
# =============================================================================
set -euo pipefail

ROOT="${ROOT:-$HOME/Desktop/BCAT1/BCAT1/CUTRUN}"
GTF="${GTF:-$ROOT/gencode.v38.annotation.gtf.gz}"
CORE="${CORE:-$ROOT/H3K4me3/reference/core.bed}"
OUT="${OUT:-$ROOT/analysis_QC_normalization/multimark_2026-07/annotation}"

ENVBIN="${ENVBIN:-$HOME/miniconda3/envs/k4width/bin}"
export PATH="$ENVBIN:$PATH"
command -v bedtools >/dev/null || { echo "[fatal] bedtools not on PATH" >&2; exit 1; }

mkdir -p "$OUT"
cd "$OUT"

# chrom sizes for chr1-22,X, taken from a project BAM header so they match the
# alignments exactly rather than a downloaded file that might differ.
if [[ ! -s genome.txt ]]; then
  samtools view -H "$ROOT/hg38/DIPG13-shNT-IgG-r1_nodup_sorted_unique.bam" \
    | awk '$1=="@SQ"{ sn=""; ln="";
             for(i=2;i<=NF;i++){ if($i ~ /^SN:/) sn=substr($i,4); if($i ~ /^LN:/) ln=substr($i,4) }
             if(sn ~ /^chr([0-9]+|X)$/) print sn"\t"ln }' \
    | sort -k1,1 > genome.txt
fi
echo "[genome] $(wc -l < genome.txt | tr -d ' ') chromosomes"

keep_main () { awk 'BEGIN{OFS="\t"} $1 ~ /^chr([0-9]+|X)$/ && $2 < $3'; }

# ---- source features ---------------------------------------------------------
[[ -s promoter.bed ]] || \
  sort -k1,1 -k2,2n "$CORE" | keep_main | cut -f1-3 | bedtools merge -i - > promoter.bed

[[ -s exon_union.bed ]] || \
  gzip -dc "$GTF" | awk 'BEGIN{OFS="\t"} $3=="exon"{print $1,$4-1,$5}' \
    | keep_main | sort -k1,1 -k2,2n | bedtools merge -i - > exon_union.bed

[[ -s gene_union.bed ]] || \
  gzip -dc "$GTF" | awk 'BEGIN{OFS="\t"} $3=="gene"{print $1,$4-1,$5}' \
    | keep_main | sort -k1,1 -k2,2n | bedtools merge -i - > gene_union.bed

# ---- disjoint categories -----------------------------------------------------
[[ -s exonic.bed ]] || bedtools subtract -a exon_union.bed -b promoter.bed > exonic.bed
[[ -s intronic.bed ]] || bedtools subtract -a gene_union.bed -b exon_union.bed \
  | bedtools subtract -a - -b promoter.bed > intronic.bed
[[ -s intergenic.bed ]] || cat promoter.bed exonic.bed intronic.bed \
  | sort -k1,1 -k2,2n | bedtools merge -i - \
  | bedtools complement -i - -g genome.txt > intergenic.bed

# ---- single partition file ---------------------------------------------------
{ awk 'BEGIN{OFS="\t"}{print $1,$2,$3,1}' promoter.bed
  awk 'BEGIN{OFS="\t"}{print $1,$2,$3,2}' exonic.bed
  awk 'BEGIN{OFS="\t"}{print $1,$2,$3,3}' intronic.bed
  awk 'BEGIN{OFS="\t"}{print $1,$2,$3,4}' intergenic.bed
} | sort -k1,1 -k2,2n > partition.bed

# ---- verify it really is a partition -----------------------------------------
echo "[check] disjoint and complete?"
awk 'BEGIN{OFS="\t"}
  { if($1==pc){ if($2 < pe) ov++; if($2 > pe) gap++ } ; pc=$1; pe=$3; n++ }
  END{ printf "  intervals=%d  overlaps=%d  gaps=%d\n", n, ov+0, gap+0 }' partition.bed

awk '{cov[$4]+=$3-$2; tot+=$3-$2} END{
  printf "  promoter   %12d bp  %5.2f%%\n", cov[1], 100*cov[1]/tot
  printf "  exonic     %12d bp  %5.2f%%\n", cov[2], 100*cov[2]/tot
  printf "  intronic   %12d bp  %5.2f%%\n", cov[3], 100*cov[3]/tot
  printf "  intergenic %12d bp  %5.2f%%\n", cov[4], 100*cov[4]/tot
  printf "  total      %12d bp\n", tot }' partition.bed
awk '{s+=$2} END{printf "  genome.txt total %d bp\n", s}' genome.txt
echo "[ok] $OUT/partition.bed"
