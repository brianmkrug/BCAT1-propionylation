# Methods

## Single-cell datasets and identification of malignant cells
Processed single-cell/single-nucleus RNA-seq count matrices were obtained for two pediatric
brain-tumor cohorts, both generated on the 10x Genomics Chromium 3′ platform and processed with
Cell Ranger (Zheng et al., 2017). Gliomas were taken from Jessa et al. (2022) (Gene Expression
Omnibus GSE210567): H3 K27M-mutant diffuse midline gliomas (DMG), H3.3-G34R/V and H3-wildtype
hemispheric high-grade gliomas, and posterior-fossa group-A ependymomas. Medulloblastomas (MB)
were taken from Riemondy et al. (2022) (GSE155446), spanning the WNT, SHH, Group 3 and Group 4
subgroups. For gliomas, malignant cells were defined by the authors' `Malignant_normal_consensus`
annotation; for MB, by the `coarse_cell_type = "malignant"` annotation. Only malignant cells were
used in all downstream analyses.

## Tumor entity/subgroup definitions
Glioma samples were assigned to entities using the molecular subgroup annotation in Supplementary
Table 1 of Jessa et al. (2022), yielding six groups: H3.1/2K27M-pons (n = 7 samples),
H3.3K27M-pons (n = 6), H3.3K27M-thalamus (n = 4), H3.3G34R/V-cortex (n = 10), H3-wildtype-cortex
(n = 4) and PFA ependymoma (n = 5). A single H3.3K27M-cortex sample was excluded as an
under-powered outlier (n = 1). MB samples were assigned to WNT (n = 1), SHH (n = 9), Group 3 and
Group 4 at the patient level; one sample of ambiguous Group 3/4 identity was excluded. The MB-WNT
group derives from a single patient and is shown for completeness only.

## Pseudobulk differential expression (basis for both panels)
Per sample, raw UMI counts were summed across malignant cells to form gene-by-sample pseudobulk
profiles (gene identifiers collapsed to symbols by summation). Differential expression between
pooled H3K27M-DMG (H3.1/2K27M-pons + H3.3K27M-pons + H3.3K27M-thalamus; n = 17 samples) and
non-DMG gliomas (H3.3G34R/V-cortex + H3-wildtype-cortex; n = 14) was computed with DESeq2 (Love,
Huber & Anders, 2014), as implemented in the Python `inmoose` package, using the negative-binomial
Wald test with default dispersion and log-fold-change shrinkage. Ependymoma and MB were excluded
from this contrast, so the test reflects within-study, glioma-internal differences. Genes detected
with ≥10 counts in ≥3 samples were retained (23,866 genes). P-values were adjusted by the
Benjamini–Hochberg procedure; genes with adjusted p < 0.05 and |log2 fold change| > 1 were
considered significant.

## Volcano plot of metabolic genes (Panel: volcano_metabolic_K27M_vs_nonDMG)
Each gene was plotted as log2 fold change (H3K27M-DMG vs non-DMG glioma) versus –log10 adjusted
p-value. Metabolic genes were annotated using the 85 curated KEGG metabolic pathways distributed
with scMetabolism (Wu et al., 2021; derived from KEGG, Kanehisa & Goto, 2000) and collapsed into
eight super-categories (Lipid, Amino acid, Carbohydrate, Nucleotide, Energy, Glycan,
Cofactor/Vitamin, Xenobiotic); genes belonging to multiple pathways were assigned to a single
super-category by a fixed priority order. Metabolic genes were colored by super-category and
non-metabolic genes shown in grey. Threshold lines mark adjusted p = 0.05 and |log2FC| = 1. For
legibility, only genes with –log10 adjusted p > 4 were labeled, with non-overlapping placement by
`adjustText`; the y-axis was clipped at 16 (a small number of highly significant non-metabolic
genes extend beyond and are drawn at the axis top). Of 1,667 KEGG-metabolic genes tested, 148 were
significant, with Lipid the largest and most up-regulated category in DMG.

## Lipid driver genes and expression dot plot (Panel: dotplot_panel_top20lipid_grouped_v2)
Lipid-metabolism genes were defined as the union of the KEGG lipid-related pathways (fatty-acid
biosynthesis/elongation/degradation, unsaturated fatty acids, glycerolipid, glycerophospholipid,
ether lipid, sphingolipid, glycosphingolipid [lacto/neolacto, globo, ganglio], sterol/steroid,
bile acid, arachidonic/linoleic/α-linolenic acid, ketone bodies and GPI-anchor biosynthesis). From
the differential-expression results above, we selected lipid genes that were (i) significantly
up-regulated in H3K27M-DMG (adjusted p < 0.05, log2FC > 0) and (ii) more highly expressed in DMG
than in both non-DMG gliomas and MB (mean pseudobulk log2 CPM), then ranked them by DESeq2
baseMean (mean normalized expression) and retained the top 20. Genes were grouped for display by
lipid pathway.

Expression was visualized across all ten tumor subgroups (six glioma + four MB) at single-cell
resolution. For each gene and subgroup, dot color encodes the mean library-size–normalized log
expression across malignant cells (log1p of counts per 10,000; per-cell total counts computed from
the full count matrix), and dot size encodes the percentage of malignant cells with non-zero
counts. A shared color scale was used across cohorts. Two genes absent from the MB count matrix
(B3GALT5, DGKE) are marked as not-detected in the MB columns.

## Statistical considerations
Differential expression was performed on pseudobulk profiles so that the biological sample (not the
individual cell) is the unit of replication, avoiding pseudoreplication inflation from correlated
cells within a tumor. DESeq2's empirical-Bayes dispersion shrinkage stabilizes variance estimates
at the modest per-group sample sizes, and multiple testing was controlled genome-wide by
Benjamini–Hochberg FDR. The differential-expression contrast (and therefore the volcano and the
gene ranking) was restricted to gliomas from a single study (Jessa et al., 2022), so it is not
confounded by cross-study batch. In contrast, the dot plot juxtaposes two independently generated
datasets (Jessa 2022 gliomas; Riemondy 2022 MB); although both used the same platform (10x
Chromium 3′, Cell Ranger), which limits technical divergence, cross-study batch effects cannot be
excluded, and glioma-versus-MB differences in the dot plot should be interpreted as descriptive /
hypothesis-generating rather than as a formal statistical comparison. The MB-WNT column (single
patient) is not statistically interpretable.

## Software and data availability
Analyses used Python (numpy, pandas, scipy); DESeq2 via `inmoose`; figures in matplotlib with
label placement by `adjustText`. Metabolic gene sets: scMetabolism / KEGG. Vector figures were
exported as SVG (text retained as editable text) and PDF (embedded TrueType). Source data:
GEO GSE210567 (gliomas) and GSE155446 (medulloblastoma).

## References
1. Jessa S, et al. K27M in canonical and noncanonical H3 variants occurs in distinct
   oligodendroglial cell lineages in brain midline gliomas. *Nat Genet.* 2022;54(12):1865–1880.
   doi:10.1038/s41588-022-01205-w.
2. Riemondy KA, et al. Neoplastic and immune single-cell transcriptomics define subgroup-specific
   intra-tumoral heterogeneity of childhood medulloblastoma. *Neuro Oncol.* 2022;24(2):273–286.
   doi:10.1093/neuonc/noab135.
3. Zheng GXY, et al. Massively parallel digital transcriptional profiling of single cells.
   *Nat Commun.* 2017;8:14049. doi:10.1038/ncomms14049.
4. Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data
   with DESeq2. *Genome Biol.* 2014;15(12):550. doi:10.1186/s13059-014-0550-8.
5. Wu Y, et al. Spatiotemporal immune landscape of colorectal cancer liver metastasis at
   single-cell level (scMetabolism). *Cancer Discov.* 2022;12(1):134–153.
   doi:10.1158/2159-8290.CD-21-0316.
6. Kanehisa M, Goto S. KEGG: Kyoto Encyclopedia of Genes and Genomes. *Nucleic Acids Res.*
   2000;28(1):27–30. doi:10.1093/nar/28.1.27.
