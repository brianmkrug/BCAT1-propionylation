# Metabolic pathway ssGSEA across pediatric brain tumor entities (Jessa et al. 2022)

First-pass analysis: curated human metabolic pathway activity (KEGG-85) across tumor-cell
pooled pseudobulk, comparing H3K27M-DMG with other pediatric brain tumor subgroups.

## Data source
- Processed scRNA/snRNA-seq from Jessa et al., *Nature Genetics* 2022 (s41588-022-01205-w).
- GEO **GSE210567** (scRNA-seq subseries): 37 samples, per-sample Cell Ranger counts matrices
  + `metadata.NG2022.tsv` (cell annotations). Raw human tumor data is EGA-controlled; these
  open processed matrices are sufficient for pseudobulk.
- Sample molecular subgroup from Supplementary Table 1 (`1-Sample_summary`).
- Gene sets: **85 KEGG metabolic pathways** (Xiao et al. 2019 / scMetabolism `KEGG_metabolism_nc.gmt`).

## Pipeline
1. Downloaded & extracted GSE210567 (37 samples x {matrix, features, barcodes, metadata}).
2. Tumor vs non-tumor split using `Malignant_normal_consensus_Jessa2022` (per-cell call).
3. Tumor-cell pseudobulk = **sum of raw counts across malignant cells**, per sample and pooled
   per entity; CPM-normalized, log2(CPM+1). Genes collapsed Ensembl->symbol.
4. **ssGSEA** (gseapy, rank normalization) on log2CPM pseudobulk -> pathway NES.
5. Row z-scored NES heatmap; K27M-DMG-vs-other contrast table.

## Entities (malignant cells pooled)
| entity | samples | malignant cells |
|---|---|---|
| HGG-H3.3G34-cortex | 10 | 38,212 |
| DMG-H3.3K27M-thalamus | 4 | 34,759 |
| DMG-H3.3K27M-pons | 6 | 27,287 |
| Ependymoma-PFA | 5 | 23,519 |
| DMG-H3.1/2K27M-pons | 7 | 9,853 |
| HGG-H3WT-cortex | 4 | 8,784 |
| HGG-H3.3K27M-cortex | 1 | 1,197 |

Total: 143,611 malignant, 33,301 normal cells. 80/85 pathways scored (5 near-absent human
pathways dropped: BCAA biosynthesis, D-Gln/D-Glu, D-Arg/D-Orn, biotin, lipoic acid).

## Validation
Marker log2CPM confirms correct entity wiring: EZHIP/CXorf67 specific to Ependymoma-PFA;
FOXG1 high in hemispheric G34/H3WT; OLIG1/2 + PDGFRA high in oligodendroglial-lineage K27M-DMG.

## Caveats
- **H3.3K27M-cortex is n=1** (1,197 cells) — treat as anecdotal, noisy on the heatmap.
- ssGSEA absolute NES is pathway-size dependent; compare *relative* (z-scored) values across entities.
- Pooled pseudobulk masks intratumoral heterogeneity; per-sample scores (`ssgsea_persample_NES.tsv`)
  are provided for variance/statistical follow-up.

## Files
- `heatmap_metabolic_ssgsea_entities.png` — main heatmap (row z-scored NES).
- `ssgsea_scores_table.tsv` — NES per entity + K27M-DMG-vs-other contrast (sorted).
- `ssgsea_perentity_NES.tsv`, `ssgsea_perentity_rowZ.tsv`, `ssgsea_persample_NES.tsv`.
- `pseudobulk_malignant_perentity_{counts,log2CPM}.tsv`, `pseudobulk_malignant_persample_{counts,log2CPM}.tsv`.
- `qc_tumor_normal_counts.csv`, `manifest.csv`.
- `scripts/` — build_manifest.py, build_pseudobulk.py, aggregate.py, run_ssgsea.py, plot_and_table.py.
- `kegg85.gmt` — gene sets used.

---

# Revision v2 — per-sample statistics, H3.3K27M-cortex excluded

Dropped the single **H3.3K27M-cortex** sample (P-1569_S-1569; n=1) to avoid sparse-data bias.
Analysis now uses **36 samples across 6 entities** (K27M-DMG: 17 samples = H3.1/2-pons 7 +
H3.3-pons 6 + H3.3-thalamus 4; other: 19 = G34 10, H3WT 4, PFA 5).

## Statistics (entity as group factor, on per-sample ssGSEA NES)
Run on the **per-sample** scores (`ssgsea_persample_NES_no_cortex.tsv`), not pooled pseudobulk,
so each entity has replicates. All p-values BH-FDR corrected across the 80 pathways.
- **Kruskal-Wallis** omnibus across the 6 entities — 10 pathways FDR<0.05.
- **limma** moderated omnibus group-effect F (inmoose port, eBayes variance shrinkage) — 44 FDR<0.05
  (more powerful/parametric than KW).
- **K27M-DMG vs rest**: Wilcoxon rank-sum (17 vs 19 samples) — 9 FDR<0.05;
  limma moderated t (equal-weighted group-mean contrast) — 4 FDR<0.05.

## K27M-DMG-vs-rest metabolic signature (concordant across tests)
Higher in K27M-DMG (positive logFC, FDR<0.05 in Wilcoxon and/or limma):
- **Sphingolipid metabolism** (limma FDR 0.004, Wilcoxon FDR 0.001) — sig in **both**
- **Glycosaminoglycan biosynthesis – chondroitin/dermatan sulfate** — sig in **both**
- **Glycerolipid metabolism** — sig in **both**
- **Ether lipid metabolism** — sig in **both**
- Glycosaminoglycan biosynthesis – heparan sulfate/heparin; Nicotinate/nicotinamide (Wilcoxon)

Lower in K27M-DMG: GAG degradation, drug metabolism–CYP450, ubiquinone/terpenoid-quinone biosynthesis.

Interpretation: a coherent **membrane-lipid + glycosaminoglycan biosynthetic program**
(sphingolipid / glycerolipid / ether-lipid + chondroitin-dermatan/heparan sulfate) is elevated in
K27M-DMG relative to hemispheric G34/H3WT gliomas and PFA ependymoma. Note the PFA-ependymoma
column is globally high across many oxidative/amino-acid pathways (distinct entity, expected).
KW is more conservative than limma; the DMG-vs-rest hits above are the method-robust core.

## New / updated files
- `heatmap_metabolic_ssgsea_entities_no_cortex.png` — revised heatmap (6 entities, sig-annotated).
- `dmg_vs_rest_contrast.png` — top pathways distinguishing K27M-DMG (limma logFC, sig marked).
- `stats_entity_pathway.tsv` — full stats table (KW, Wilcoxon, limma t + omnibus F, FDRs, effect sizes).
- `ssgsea_persample_NES_no_cortex.tsv`, `ssgsea_entitymean_NES_no_cortex.tsv`.
- `scripts/stats_analysis.py`, `scripts/limma_patch.py`, `scripts/final_plots.py`.

---

# Revision v3 — medulloblastoma added as an entity group

## Dataset selection (platform-matched to Jessa 2022)
Goal: add MB across all 4 subgroups using data generated like Jessa (10X Chromium 3' + Cell Ranger)
to avoid technology-driven artifacts.
- **Chosen: Riemondy et al., Neuro-Oncology 2022 — GSE155446** (PMID 34077540). 10X Chromium 3',
  Cell Ranger, NovaSeq 6000; 28+ childhood MB (WNT/SHH/GP3/GP4); open processed counts +
  `coarse_cell_type` (malignant call) + subgroup. 95 MB counts matrix — compact and directly usable.
- **Rejected: Hovestadt et al., Nature 2019 — GSE119926.** Smart-seq2 (plate-based, full-length) —
  a platform MISMATCH that would confound the comparison with technical artifacts.

## MB processing (same pipeline)
Genes×cells counts streamed in row-chunks (1.5 GB CSV) and matrix-multiplied by a malignant-cell→
sample indicator to get per-sample malignant pseudobulk (18,868 symbols × 30 samples; 220M counts).
Per-sample subgroup assigned at patient level (GP3/4 cells folded in); 1 fully-intermediate sample
(GP3/4) excluded; 1 WNT sample retained but flagged (n=1). CPM→log2 → ssGSEA (KEGG-85).
Marker validation confirmed labels: SHH=GLI1/ATOH1/HHIP; Group3=MYC/NRL/IMPG2 (photoreceptor);
Group4=KCNA1/EOMES/UNC5D; WNT=WIF1/DKK2/LEF1.

## Combined analysis
Per-sample ssGSEA NES from both studies concatenated (rank-based within-sample, so cross-study-safe
at the scoring step): **77 shared pathways × 65 samples across 10 entities** (6 glioma + MB-WNT/SHH/
Group3/Group4). Entity heatmap = row z-scored mean NES.

Signal (limma + Kruskal-Wallis, FDR<0.05): omnibus group effect in ~75/77 pathways.
- **Glioma > MB**: chondroitin/dermatan & heparan GAG biosynthesis, ganglio-glycosphingolipid,
  fatty-acid biosynthesis, glycerolipid/ether-lipid, O-/N-glycan biosynthesis (membrane-lipid/glycan
  biosynthetic program).
- **MB > glioma**: oxidative phosphorylation, TCA cycle, pyruvate metabolism, glutathione,
  cysteine/methionine (oxidative/mitochondrial program).

## IMPORTANT caveat — cross-study batch
Glioma (Jessa) and MB (Riemondy) come from different studies/labs. Same platform (10X/Cell Ranger)
minimizes but does NOT remove batch; the glioma-vs-MB axis is much larger in magnitude (logFC ~0.1–0.2)
than within-study contrasts (DMG-vs-rest ~0.03–0.06), so **glioma-vs-MB differences are partly
technical and should be treated as hypothesis-generating**. The trustworthy comparisons remain the
within-study ones: DMG-vs-other-glioma (Jessa) and MB subgroup-vs-subgroup (Riemondy). WNT is n=1.

## New files (v3)
- `heatmap_all_entities_glioma_MB.png` — 10-entity metabolic heatmap.
- `glioma_vs_mb_contrast.png` — top glioma-vs-MB pathways (batch-caveated).
- `combined_persample_NES.tsv`, `combined_entitymean_NES.tsv`, `combined_sample_entity.csv`.
- `stats_combined_entity_pathway.tsv` — KW, limma omnibus-F, glioma-vs-MB (Wilcoxon + limma).
- `mb_pseudobulk_malignant_persample_{counts,log2CPM}.tsv`, `mb_ssgsea_persample_NES.tsv`,
  `mb_sample_subgroup.csv`.
- `scripts/mb_part_pb.py`, `scripts/mb_finish.py`, `scripts/combine_and_stats.py`, `scripts/plot_combined.py`.

---

# Revision v4 — leading-edge driver genes (lipid pathways) + violin by subgroup

Identified genes driving the K27M-DMG enrichment of the three DMG-up lipid pathways via **GSEA
prerank** on the within-glioma DMG-vs-other-glioma per-gene limma-t ranking (batch-free), then
extracted the leading edge for each set.

Top leading-edge drivers (gene: DMG-vs-other-glioma limma t):
- **Sphingolipid metabolism**: UGT8 (t=9.1), GAL3ST1 (t=6.2), NEU4, SMPD3, UGCG — UGT8/GAL3ST1 are
  galactosylceramide-synthase / cerebroside-sulfotransferase (myelin sphingolipids), fitting the
  oligodendroglial lineage of K27M-DMG.
- **Glycerolipid metabolism**: DGAT2 (t=6.3), DGKG (t=4.1), MBOAT2, LPIN1, DGKE.
- **Glycosphingolipid biosynthesis (lacto/neolacto)**: ST8SIA1 (GD3 synthase, t=3.2), ST3GAL6, FUT7.

`violin_leadingedge_genes_by_subgroup.png`: per-**malignant-cell** log-normalized expression
(log1p CP10k) of 6 chosen drivers (2 per pathway), split by the 6 glioma subgroups; DMG subgroups
in red, mean dot + %-expressing annotated. UGT8/GAL3ST1 are near-exclusively expressed in DMG
(e.g. UGT8 15–28% of DMG cells vs 0–3% in other gliomas). Per-cell values extracted directly from
the Jessa Cell Ranger matrices for malignant cells.

New files: `violin_leadingedge_genes_by_subgroup.png`, `leadingedge_expression_by_subgroup.csv`
(mean lognorm + %expressing per gene per subgroup), `leading_edge_genes.json`,
`scripts/leading_edge.py`, `scripts/extract_cells.py`, `scripts/plot_violin.py`.

---

# Revision v5 — dot-plot of driver genes across ALL subgroups (glioma + MB)

Added the 4 MB subgroups to the leading-edge driver-gene view and switched to a **dot-plot**
(dot size = % malignant cells expressing; color = mean log-norm expression), which handles the
zero-inflation of these enzymes better than violins and shows the scope of differences compactly.

MB per-cell values computed from GSE155446 raw counts: per-cell library sizes (Σ counts) obtained
by streaming the counts matrix; the 5 available driver genes extracted and normalized as log1p(CP10k)
over malignant cells, grouped by subgroup. **GAL3ST1 is not present in the MB counts matrix** (marked
with × / NA).

Readout: the sphingolipid/glycosphingolipid drivers (UGT8, ST8SIA1, ST3GAL6) are largest and warmest
in K27M-DMG (pons/thalamus); MB subgroups are mostly low, except **MB-SHH** which shows moderate UGT8
(14% cells) and ST8SIA1 (18%). DGKG (glycerolipid) is broadly expressed across all entities.
Glioma-vs-MB magnitudes remain cross-study (batch caution).

New files: `dotplot_leadingedge_all_subgroups.png`, `dotplot_expression_all_subgroups.csv`
(mean lognorm + %expressing per gene per subgroup, glioma + MB), `scripts/dotplot.py`,
`scripts/mb_colsums.py`.

---

# Revision v6 — BCAT1/BCAT2 added to driver dot-plot

Added BCAT1 and BCAT2 (BCAA transaminases) at the top of the driver-gene dot-plot, above the lipid-
pathway leading-edge genes. Per-cell expression extracted the same way (glioma from Jessa Cell Ranger
matrices; MB from GSE155446 raw counts, log1p CP10k over malignant cells). **BCAT2 is absent from the
MB counts matrix** (marked × / NA, like GAL3ST1).

Readout: BCAT1 is broadly expressed in gliomas (26–31% of malignant cells across DMG + hemispheric
HGG) but near-absent in ependymoma-PFA (2%); among MB it is restricted to MB-SHH (39%) and low in
WNT/Group3/Group4 (~4%). BCAT2 is lower/more uniform across gliomas (8–17%).

New files: `dotplot_BCAT_leadingedge_all_subgroups.png`,
`dotplot_expression_all_subgroups_withBCAT.csv`, `scripts/dotplot_bcat.py`, `scripts/extract_bcat.py`.

---

# Revision v7 — pseudobulk DE (H3K27M vs non-K27M glioma) + robust lipid gene swap

**DE test:** pooled H3K27M gliomas (n=17: H3.1/2K27M-pons, H3.3K27M-pons, H3.3K27M-thalamus) vs
non-K27M gliomas (n=14: H3.3G34 + H3WT cortex; EPN excluded as non-glioma, MB excluded). Malignant
per-sample pseudobulk counts, **DESeq2** (inmoose port), Wald test, BH-FDR. 23,866 genes tested.
**2,632 genes padj<0.05 (1,309 UP in H3K27M, 1,323 down).** Full ranked list:
`DE_H3K27M_vs_nonK27M_glioma_ranked.csv`; lipid-only up-DEGs: `DE_lipid_up_in_K27M.tsv`.

Top up-in-H3K27M DEGs (padj≈0): SFRP2, FAM162B, SEMA3D, ACADL, KCNB1, DNER, ACAN, LY6H, C1QL1, TNS3
(OPC/oligodendroglial-lineage program).

**DGAT2 / GAL3ST1 were lowly expressed** (baseMean 44 / 120). Replaced in the dot-plot with two
robustly-expressed, K27M-enriched (vs other gliomas AND MB) lipid genes:
- **SCD5** (stearoyl-CoA desaturase 5; unsat-FA): baseMean 3090, expressed in 53–71% of DMG cells.
- **PTGDS** (prostaglandin-D2 synthase; oligodendrocyte-associated): baseMean 1521, LFC 2.6,
  28–41% of DMG cells vs 4–5% in hemispheric HGG (sharpest glioma-internal gap).

Other strong robust candidates (see `DE_lipid_up_in_K27M.tsv`): DGKI, AGPAT4, ST3GAL5 (GM3 synthase),
ST6GALNAC3, ACSL4, PLA2G4A, and **NEU4** (most K27M-specific: near-zero in other gliomas + MB).

New files: `dotplot_BCAT_robustlipid_all_subgroups.png`, `dotplot_expression_all_subgroups_v2.csv`,
`DE_H3K27M_vs_nonK27M_glioma_ranked.csv`, `DE_lipid_up_in_K27M.tsv`, `scripts/de_k27m.py`,
`scripts/dotplot_v2.py`, `scripts/extract_new.py`.

---

# Revision v8 — expanded panel: top-10 robust lipid DEGs

Larger dot-plot panel: BCAT1/BCAT2 + the **top-10 robustly-expressed lipid DEGs** up in H3K27M
(significant padj<0.05, expressed above BOTH non-K27M gliomas and MB, ranked by baseMean):
SCD5, DGKI, DGKG, ST8SIA1, PTGDS, AGPAT4, FDFT1, ST3GAL5, ST6GALNAC3, ACSL4 (all present in MB).
`top10_robust_lipid_DEGs.tsv` lists baseMean/LFC/padj/per-group log2CPM; per-cell dot data in
`dotplot_panel_top10lipid.csv`.

Pathways represented: unsaturated-FA (SCD5), glycerolipid/phospholipid (DGKI, DGKG, AGPAT4),
glycosphingolipid (ST8SIA1, ST3GAL5, ST6GALNAC3), arachidonic/prostaglandin (PTGDS), sterol (FDFT1),
fatty-acid activation/degradation (ACSL4). All show large, warm dots across the K27M-DMG columns.

New files: `dotplot_panel_top10lipid_all_subgroups.png`, `dotplot_panel_top10lipid.csv`,
`top10_robust_lipid_DEGs.tsv`, `scripts/dotplot_panel.py`, `scripts/extract_panel.py`.

---

# Revision v9 — top-20 lipid DEG panel, grouped by pathway

Extended the dot-plot to BCAT1/2 + the **top-20 robustly-expressed lipid DEGs** up in H3K27M
(padj<0.05, > non-K27M glioma & > MB, ranked by baseMean): SCD5, DGKI, DGKG, ST8SIA1, PTGDS, AGPAT4,
FDFT1, ST3GAL5, ST6GALNAC3, ACSL4, PIGM, UGT8, PLA2G4A, NEU4, B3GALT5, DGKE, CYP2U1, CYP39A1, ENPP2,
DHCR24. Genes are **grouped by lipid pathway** with labeled color boxes on the left: BCAA transaminase,
Fatty acid, Glycero-/phospholipid, Sphingolipid, Glycosphingolipid, Eicosanoid/arachidonic,
Sterol/bile acid, GPI-anchor. Dot spacing increased to avoid overlap. B3GALT5 and DGKE are absent
from the MB matrix (× / NA).

Files: `dotplot_panel_top20lipid_grouped.png`, `dotplot_panel_top20lipid.csv`,
`top20_robust_lipid_DEGs.tsv`, `scripts/dotplot_panel20.py`, `scripts/extract_panel2.py`.

---

# Revision v10 — volcano plot: metabolic genes among DMG-vs-non-DMG DEGs

Volcano of the pseudobulk DE (H3K27M-DMG n=17 vs non-DMG glioma [G34+H3WT] n=14; EPN/MB excluded),
with KEGG metabolic genes colored by super-category (Lipid, Amino acid, Carbohydrate, Nucleotide,
Energy, Glycan, Cofactor/Vitamin, Xenobiotic) over a grey background of all genes. Thresholds
padj<0.05 & |log2FC|>1; labels via adjustText; y-axis zoomed to 16 (few non-metabolic genes extend
higher). 1,667 KEGG-metabolic genes tested; 148 significant.

Readout — metabolic genes ARE well represented among DEGs, and **Lipid is the largest and most
up-skewed metabolic category** (50 sig, 35 up in DMG: ACADL, DGAT2, GAL3ST1, UGT8, AKR1C1/2, ACSL4,
SCD5, ST8SIA1, PTGDS, FDFT1, DGKI; down: ELOVL2, GPX8). Glycan also up-skewed (16/22 up: GALNT/CHST/
C1GALT1/ST-transferases). Amino-acid genes skew down (GAD2, SHMT1, ANPEP). Full list:
`significant_metabolic_DEGs.csv`.

New files: `volcano_metabolic_K27M_vs_nonDMG.png` + `.svg`, `significant_metabolic_DEGs.csv`,
`gene_metabolic_category.json`, `scripts/volcano.py`.

---

# Revision v11 — lipid-only volcano

Second volcano variant: only **Lipid**-category genes colored (red); all other genes medium grey
(more visible than v10's light grey). All 50 lipid genes crossing padj<0.05 & |log2FC|>1 are labeled
with adjustText repulsion + leader lines (35↑ in DMG, 15↓). Files: `volcano_lipid_K27M_vs_nonDMG.png`
+ `.svg`, `scripts/volcano_v2.py`.

---

# Revision v12 — Illustrator-editable vector outputs

Patched all 16 figure-generating scripts to keep text editable in Adobe Illustrator: inserted
`plt.rcParams['svg.fonttype']='none'` and `plt.rcParams['pdf.fonttype']=42` immediately after the
pyplot import (plus a savefig shim so each figure is written as PNG + SVG + PDF). Re-ran every script
and regenerated all outputs. SVGs now store text as live `<text>` with font-family (not outlined
paths); PDFs embed TrueType (`/FontFile2`, no Type3), so labels/genes/axes are selectable and editable.
All 14 figures now have matching .png/.svg/.pdf. (`final_plots.py` also made idempotent.)

---

# References (data sources)

1. **Gliomas (H3K27M-DMG, H3.3G34, H3WT, PFA ependymoma):** Jessa S, Mohammadnia A, Harutyunyan AS,
   et al. K27M in canonical and noncanonical H3 variants occurs in distinct oligodendroglial cell
   lineages in brain midline gliomas. *Nature Genetics.* 2022;54(12):1865–1880.
   doi:10.1038/s41588-022-01205-w. Data: GEO GSE210567 (scRNA-seq).

2. **Medulloblastoma (WNT, SHH, Group 3, Group 4):** Riemondy KA, Venkataraman S, Willard N, et al.
   Neoplastic and immune single-cell transcriptomics define subgroup-specific intra-tumoral
   heterogeneity of childhood medulloblastoma. *Neuro-Oncology.* 2022;24(2):273–286.
   doi:10.1093/neuonc/noab135. PMID: 34077540. Data: GEO GSE155446.
