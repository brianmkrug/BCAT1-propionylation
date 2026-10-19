# xHSJ019 xenograft — H3K27M mutant vs BCAT1/H3-KO: DE + GSEA

**Comparison:** K27M (n=3: 167, 172, 616) vs KO (n=3: 613, 614, 618). Positive log2FC / NES = higher in **K27M**.
**Pipeline:** Salmon transcript quant (hg38, Ensembl/GENCODE ENST) → gene-level counts → PyDESeq2 → gseapy prerank (ranked by DESeq2 Wald *stat*).

## Methods
- **Tx→gene:** built from Ensembl GRCh38 release-110 GTF, matched on version-stripped ENST IDs (93.3% of quantified transcripts mapped). Gene-level counts = sum of Salmon estimated reads (`NumReads`) per gene, rounded → 35,604 genes.
- **DESeq2:** design `~condition`, contrast K27M vs KO. Prefilter total count ≥10 → 18,918 genes tested. Median-of-ratios normalization (size factors 0.41–2.47, reflecting the wide library-size range).
- **GSEA:** prerank, 1000 permutations, gene-set size 5–500, seed 42, ranked by Wald statistic. Gene sets: **MSigDB Hallmark 2020** (50) and a **curated KEGG "Metabolism" set** (78 pathways spanning carbohydrate, energy, lipid, nucleotide, amino-acid, glycan, cofactor/vitamin, terpenoid and xenobiotic metabolism; `KEGG_metabolism_expanded.gmt`).

## Data quality
PCA cleanly separates K27M from KO along PC1 (K27M negative, KO positive). BCAT1 is strongly up in K27M (log2FC +1.62, padj ≈ 0) — consistent with H3K27M biology.

## Differential expression
3,178 genes at padj<0.05; with |log2FC|>1: **971 up in K27M, 1,062 up in KO**. Top hits include ANXA1, IRX4, PKIB, GNG12, LPL, CD109 (up) and PEG3, STOM, SIX1, FOXD1 (up in KO).

## GSEA — Hallmark (42 sets at FDR<0.25)
K27M is dominated by **proliferation** (E2F Targets NES 3.13, G2-M Checkpoint 2.80, Myc Targets V1/V2, Mitotic Spindle, DNA Repair), **inflammatory/immune signaling** (TNFα-NF-κB 2.57, IL-6/JAK/STAT3, IFN-γ/α, IL-2/STAT5), **EMT** (2.37), and a coherent **anabolic-metabolic block**: **mTORC1 Signaling (2.17), Cholesterol Homeostasis (2.18), Hypoxia, Glycolysis (1.64), Fatty Acid Metabolism (1.53)**. Notably BCAT1 sits in the mTORC1 leading edge. Only Pancreas Beta Cells and KRAS Signaling Dn are enriched in KO.

## GSEA — KEGG metabolism
The metabolic signal is concentrated in the **mevalonate → sterol/cholesterol biosynthesis axis**:

| Pathway | NES | FDR |
|---|---|---|
| Steroid biosynthesis | 2.03 | 0.005 |
| Terpenoid backbone biosynthesis (mevalonate) | 1.98 | 0.010 |
| One carbon pool by folate | 1.73 | 0.098 |
| Selenocompound metabolism | 1.63 | 0.161 |
| Glutathione metabolism | 1.45 | 0.41 (nom p 0.04) |

Leading edges: **Steroid biosynthesis** — SQLE, FDFT1, CYP51A1, MSMO1, DHCR24, NSDHL, EBP, LSS, SC5D; **Terpenoid backbone** — HMGCS1, MVD, FDPS, IDI1, PMVK, ACAT2. BCAA-relevant **Valine/leucine/isoleucine degradation** trends up (NES 1.24, nom p 0.14) but is not FDR-significant; fatty-acid degradation, cholesterol metabolism and one-carbon/folate all trend up. OXPHOS and TCA are flat/slightly negative.

## Interpretation (for the BCAT1 program)
K27M cells run a highly proliferative, anabolic program: mTORC1 activation coupled to de-novo **cholesterol/mevalonate synthesis**, one-carbon/folate and glutathione (redox) metabolism, glycolysis and lipogenic fatty-acid metabolism — with BCAT1 itself elevated. This is consistent with a BCAA→acyl-CoA/lipogenic axis feeding sterol and fatty-acid biosynthesis. The mevalonate/SREBP output is the clearest, most robust metabolic difference.

## Files
- `gene_counts_raw.csv` — gene × sample raw counts
- `deseq2_K27M_vs_KO.csv` — full DE table (baseMean, log2FoldChange, stat, pvalue, padj)
- `ranked_genes.rnk` — GSEA input ranking
- `gsea_hallmark_results.csv`, `gsea_kegg_metab_results.csv` — full GSEA tables
- `KEGG_metabolism_expanded.gmt` — the curated metabolism gene sets
- `figures/` — PCA, volcano, Hallmark & KEGG NES barplots

*Caveats: n=3/group; gene-level counts use simple `NumReads` summation (no tximport length offset); GSEA FDR for the metabolism set reflects the smaller 78-pathway universe.*
