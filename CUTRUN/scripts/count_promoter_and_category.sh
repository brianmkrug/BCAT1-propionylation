#!/usr/bin/env bash
# =============================================================================
# One pass per library, producing BOTH:
#   (a) per-promoter fragment counts over the 36,689 core.bed windows (TSS +/- 1 kb)
#   (b) fragment counts in the 4 disjoint genomic categories
#       (promoter / exonic / intronic / intergenic)
#
# Same fragment definition as every other count in this project:
#   -f 66 -F 1804, MAPQ >= 30, chr1-22 + chrX, GenBank->UCSC rename for H3K27ac.
#
# Fragments are assigned by MIDPOINT, not by read1 5' end. At 10 kb that choice was
# immaterial; over a 2 kb promoter window a ~130 bp offset is not, so the midpoint
# is computed from TLEN:
#     TLEN > 0 : mid = POS-1 + TLEN/2
#     TLEN < 0 : mid = POS-1 + reflen + TLEN/2        (read1 is the rightmost mate)
#
# Categories are disjoint, so each fragment contributes to exactly one -- no
# double counting. Promoter windows DO overlap (up to 7 deep), so a fragment can
# count toward several promoters; that is correct for per-promoter quantification
# and is why the two structures are searched separately.
#
# RUN DETACHED:
#   cd ${BCAT1_REPO:-.}/CUTRUN
#   nohup ./scripts/count_promoter_and_category.sh > logs/promcat_$(date +%F_%H%M).log 2>&1 &
# =============================================================================
set -euo pipefail

ROOT="${ROOT:-$HOME/Desktop/BCAT1/BCAT1/CUTRUN}"
BAMDIR="${BAMDIR:-$ROOT/hg38}"
QCDIR="${QCDIR:-$ROOT/analysis_QC_normalization/multimark_2026-07}"
ANN="${ANN:-$QCDIR/annotation}"
CORE="${CORE:-$ROOT/H3K4me3/reference/core.bed}"
OUT="${OUT:-$QCDIR/promcat}"

ENVBIN="${ENVBIN:-$HOME/miniconda3/envs/k4width/bin}"
[[ -x "$ENVBIN/samtools" ]] || { echo "[fatal] no samtools at $ENVBIN" >&2; exit 1; }
export PATH="$ENVBIN:$PATH"

if PCORES=$(sysctl -n hw.perflevel0.logicalcpu 2>/dev/null); then :; else
  PCORES=$(sysctl -n hw.logicalcpu 2>/dev/null || nproc)
fi
JOBS="${JOBS:-3}"
THREADS="${THREADS:-$(( PCORES / JOBS ))}"; (( THREADS < 2 )) && THREADS=2
echo "[cfg] jobs=$JOBS threads/job=$THREADS"

mkdir -p "$OUT" "$ROOT/logs"
[[ -s "$ANN/partition.bed" ]] || { echo "[fatal] run scripts/build_genome_partition.sh first" >&2; exit 1; }

# promoter windows, sorted, with their ORIGINAL core.bed line number as id so the
# output joins back to the reference without ambiguity
PROMS="$OUT/promoters_sorted.bed"
[[ -s "$PROMS" ]] || awk 'BEGIN{OFS="\t"}{print $1,$2,$3,NR}' "$CORE" | sort -k1,1 -k2,2n > "$PROMS"

GB2UCSC="CM000663.2:chr1 CM000664.2:chr2 CM000665.2:chr3 CM000666.2:chr4 \
CM000667.2:chr5 CM000668.2:chr6 CM000669.2:chr7 CM000670.2:chr8 CM000671.2:chr9 \
CM000672.2:chr10 CM000673.2:chr11 CM000674.2:chr12 CM000675.2:chr13 \
CM000676.2:chr14 CM000677.2:chr15 CM000678.2:chr16 CM000679.2:chr17 \
CM000680.2:chr18 CM000681.2:chr19 CM000682.2:chr20 CM000683.2:chr21 \
CM000684.2:chr22 CM000685.2:chrX"

MANIFEST="$QCDIR/sample_manifest.tsv"
[[ -s "$MANIFEST" ]] || { echo "[fatal] missing $MANIFEST" >&2; exit 1; }

count_one () {
  local label=$1 naming=$2 src=$3
  local outp="$OUT/${label}.promoter.tsv.gz"
  local outc="$OUT/${label}.category.tsv"
  [[ -s "$outp" && -s "$outc" ]] && { echo "[skip] $label"; return; }
  local bam="$BAMDIR/${src}"
  [[ -s "$bam" ]] || { echo "[fatal] missing $bam" >&2; return 1; }

  local regions map
  if [[ $naming == gb ]]; then
    map="$GB2UCSC"
    regions=$(for p in $GB2UCSC; do printf '%s ' "${p%%:*}"; done)
  else
    map=""
    regions=$(printf 'chr%s ' {1..22} X)
  fi

  echo "[count] $label ($naming)"
  # shellcheck disable=SC2086
  samtools view -@ "$THREADS" -q 30 -f 66 -F 1804 "$bam" $regions \
    | awk -v MAP="$map" -v PART="$ANN/partition.bed" -v PROM="$PROMS" \
          -v OUTC="$outc.partial" '
      function bsearch(lo, hi, m,   a, b, mid) {
        # last interval whose start <= m, within [lo,hi)
        a = lo; b = hi - 1
        while (a <= b) { mid = int((a+b)/2)
          if (ps[mid] <= m) a = mid + 1; else b = mid - 1 }
        return b            # -1 if none
      }
      BEGIN{
        FS=OFS="\t"
        n = split(MAP, kv, " ")
        for (i=1; i<=n; i++) { split(kv[i], p, ":"); ren[p[1]] = p[2] }

        # ---- partition: flat arrays, per-chrom [lo,hi) ----
        np = 0
        while ((getline line < PART) > 0) {
          split(line, f, "\t")
          ps[np] = f[2]+0; pe[np] = f[3]+0; pc[np] = f[4]+0
          if (!(f[1] in plo)) plo[f[1]] = np
          phi[f[1]] = np + 1
          np++
        }
        close(PART)

        # ---- promoters: separate flat arrays, per-chrom [lo,hi) ----
        nq = 0
        while ((getline line < PROM) > 0) {
          split(line, f, "\t")
          qs[nq] = f[2]+0; qe[nq] = f[3]+0; qid[nq] = f[4]+0
          if (!(f[1] in qlo)) qlo[f[1]] = nq
          qhi[f[1]] = nq + 1
          nq++
        }
        close(PROM)
        MAXW = 2000     # every core.bed window is exactly 2 kb
      }
      {
        c = ($3 in ren) ? ren[$3] : $3
        t = $9 + 0
        if (t > 0) mid = $4 - 1 + int(t/2)
        else       mid = $4 - 1 + length($10) + int(t/2)

        # ---- category (disjoint: exactly one hit) ----
        if (c in plo) {
          # swap the partition arrays into the search by copying bounds
          lo = plo[c]; hi = phi[c]
          a = lo; b = hi - 1; got = -1
          while (a <= b) { m2 = int((a+b)/2)
            if (ps[m2] <= mid) { got = m2; a = m2 + 1 } else b = m2 - 1 }
          if (got >= 0 && mid < pe[got]) cat[pc[got]]++
        }

        # ---- promoters (overlapping: walk back while still in range) ----
        if (c in qlo) {
          lo = qlo[c]; hi = qhi[c]
          a = lo; b = hi - 1; got = -1
          while (a <= b) { m2 = int((a+b)/2)
            if (qs[m2] <= mid) { got = m2; a = m2 + 1 } else b = m2 - 1 }
          for (i = got; i >= lo && qs[i] > mid - MAXW; i--)
            if (mid < qe[i]) hits[qid[i]]++
        }
        tot++
      }
      END{
        printf "category\tfragments\n"                        > OUTC
        printf "promoter\t%d\n",   cat[1]+0                   > OUTC
        printf "exonic\t%d\n",     cat[2]+0                   > OUTC
        printf "intronic\t%d\n",   cat[3]+0                   > OUTC
        printf "intergenic\t%d\n", cat[4]+0                   > OUTC
        printf "total\t%d\n",      tot+0                      > OUTC
        close(OUTC)
        for (k in hits) print k, hits[k]
      }' \
    | sort -k1,1n | gzip -c > "$outp.partial"

  mv "$outc.partial" "$outc"
  mv "$outp.partial" "$outp"
  echo "[done] $label  $(awk -F'\t' '$1=="total"{print $2}' "$outc") fragments"
}
export -f count_one
export OUT BAMDIR ANN PROMS THREADS GB2UCSC

awk -F'\t' 'NR>1{print $1"\t"$5"\t"$6}' "$MANIFEST" \
  | xargs -P "$JOBS" -n3 bash -c 'count_one "$0" "$1" "$2"'

echo "[ok] per-promoter counts and category counts in $OUT"
