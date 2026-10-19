# BCAT1 knockout brain RNA-seq — differential expression & cell-type GSEA

Analysis of bulk RNA-seq from *Bcat1* knockout (KO) vs wild-type (WT) mouse brain
(hindbrain) at two ages, **P6** and **P90**, run separately per age. Goal: identify
gene-expression changes caused by BCAT1 loss and interpret them as shifts in cell-type
identity/composition using published cell-type signatures.

## Data

Salmon quant (`salmon_quant/BCAT1-KO/`), mm10 / GENCODE vM25 index, transcript-level.

| Age | KO replicates | WT replicates |
|-----|---------------|---------------|
| P6  | KOr2, KOr3, KOr4 (3) | WTr1–WTr4 (4) |
| P90 | KOr1–KOr4 (4) | WTr1–WTr4 (4) |

(P6 KOr1 was not present in the folder.)

## Methods

1. **Transcript→gene**: GENCODE vM25 GTF built ENSMUST→gene-symbol map (99.4% of quantified
   transcripts matched). Estimated counts (`NumReads`) summed to gene level.
2. **Differential expression**: DESeq2 (PyDESeq2 0.5.2), design `~condition`, contrast
   **KO vs WT**, run independently for each age. Genes with ≥10 total reads kept.
   Ranking metric for GSEA = Wald statistic (`stat`).
3. **GSEA**: preranked (gseapy 1.3.0, 1000 permutations, set size 10–500) against each
   signature collection. Positive NES = enriched in KO; negative NES = depleted in KO.
4. **Cross-species mapping**: human signatures translated to mouse symbols via the MGI
   mouse–human ortholog table (HOM_MouseHumanSequence), so all sets were scored against
   the mouse ranked list.

## Knockout QC — passes

*Bcat1* is the single strongest down-regulated gene at both ages, confirming sample labels:

| Age | Bcat1 log2FC (KO/WT) | padj |
|-----|----------------------|------|
| P6  | −3.05 | 5e-175 |
| P90 | −3.79 | 8e-95 |

Residual ~10% signal in KO is typical (partial transcript / retained exons). Global PCA
(`figures/QC_PCA.png`) shows KO/WT partially separating — the whole-tissue transcriptional
effect is modest, which is why preranked GSEA (sensitive to coordinated small shifts) is the
appropriate lens.

Numbers of strongly DE genes (padj<0.05, |log2FC|>1): **P6 = 32** (9 up / 23 down);
**P90 = 15** (10 up / 5 down).

## Signature collections (one heatmap each)

- **Kleinman Braindex — mouse** (`mm_sym`, 290 developmental brain-atlas cell-type clusters)
- **Kleinman Braindex — human** (`hg_sym`, 290 clusters → mouse orthologs)
- **MSigDB m8** — mouse cell-type signatures (233 sets)
- **MSigDB c8** — human cell-type signatures (866 sets → mouse)
- **MSigDB MH** — mouse Hallmark (50 sets)
- **MSigDB H** — human Hallmark (50 sets → mouse)

Braindex cluster labels: prefix `F-`/`P-` = forebrain/pons region; `e10..p6` = embryonic/
postnatal age of the reference cell; suffix `-P` = proliferating/cycling version of the
cell state. See Jessa et al. *Nat. Genet.* 2022 (Braindex) and Jessa et al. *Cancer Res.*
2025 (FOXR2 supplement) for full cluster definitions.

## Headline findings (consistent across independent collections)

- **Ependymal identity is gained in KO at both ages.** Braindex EPEN clusters (F-p6_EPEN,
  P-p6_EPEN2, P-p3_EPEN1) and `DESCARTES_ORGANOGENESIS_EPENDYMAL_CELL` (m8) are positively
  enriched at P6 and P90.
- **Oligodendrocyte-lineage maturation is reduced in KO at P90.** Braindex NFOL (newly-formed
  oligodendrocyte), MOL (mature oligodendrocyte) and `TABULA_MURIS_..._OLIGODENDROCYTE`
  signatures are depleted at P90.
- **Age-dependent reversal of proliferating progenitor programs.** Many `-P` (cycling
  progenitor/RGC/OPC-P/NEURP) Braindex signatures are *depleted* at P6 but *enriched* at P90
  in KO — an age-flipped shift in the balance of cycling neural progenitor states.

Interpretation is left to you; NES/FDR for every set is in `out/GSEA_all_results_long.csv`.

## Metabolism-focused GSEA (P6)

To go beyond broad Hallmark classes, preranked GSEA was run on the **P6** ranked list against
two metabolism-specific collections:

- **KEGG mouse metabolism** — 95 pathways downloaded fresh from the KEGG REST API
  (`rest.kegg.jp`, mmu metabolism category, IDs mmu00xxx/01xxx), native mouse symbols, no
  ortholog step. GMT: `gmt/KEGG_metabolism_mouse.gmt`.
- **Curated 47-pathway metabolic panel** — the same panel used in the Velmeshev snRNA-seq
  metabolic-landscape analysis (`human_data/metabolic_panel.json`, KEGG + Hallmark + custom
  modules), translated to mouse via the MGI ortholog table. GMT:
  `gmt/Curated_metabolic_panel_mouse.gmt`.

Both collections give the same picture (28/95 KEGG and 29/47 curated sets at FDR<0.25):

- **Up in KO** (energy / central-carbon metabolism): oxidative phosphorylation, TCA cycle,
  glycolysis/gluconeogenesis, carbon metabolism, pyruvate, pentose phosphate, butanoate, and
  fatty-acid synthesis — a broad up-shift in mitochondrial/energy metabolism.
- **Down in KO**: **Valine/leucine/isoleucine (BCAA) degradation** and 2-oxocarboxylic-acid
  metabolism (the transamination substrates/products of BCAT1) are down — the direct
  metabolic signature of losing the BCAA transaminase — together with pantothenate/CoA
  biosynthesis, lysine degradation, cysteine/methionine metabolism, mitochondrial biogenesis,
  one-carbon/folate and mTORC1.

The coordinated down-regulation of BCAA catabolism is a mechanistic validation independent of
the *Bcat1* transcript knockdown itself. Full ranked table: `out/Metabolism_GSEA_P6_ranked.csv`.

Figures: `figures/metab_KEGG_p6_ranked.png`, `figures/metab_curated_p6_ranked.png`
(NES ranked; red = up in KO, blue = down; `*` = FDR q<0.25; ◀ marks BCAA-related pathways).

### P6 vs P90 side-by-side, and driver genes

The same GSEA was run on **P90** (`gsea/*_p90.csv`) and placed side-by-side
(`figures/metab_KEGG_P6vsP90_heatmap.png`). Two distinct temporal patterns emerge:

- **Transient P6 energy up-shift.** OXPHOS, TCA, glycolysis, carbon metabolism and fatty-acid
  biosynthesis are strongly up in KO at P6 but *reverse to neutral/down* by P90 — a
  compensatory energy-metabolism response specific to the early postnatal brain.
- **Durable, strengthening amino-acid catabolism deficit.** BCAA (Val/Leu/Ile) degradation,
  2-oxocarboxylic-acid metabolism, Cys/Met metabolism, pantothenate/CoA biosynthesis and
  amino-acid/cofactor biosynthesis are down at both ages and **more strongly down at P90**
  (BCAA degradation NES −1.4 → −2.0, FDR<0.05 at P90). This is the persistent metabolic
  signature of losing the BCAA transaminase.

**Leading-edge driver genes** for the top 10 P6 pathways are in
`out/P6_leading_edge_top_pathways.json` (per pathway) and `out/P6_leading_edge_genes_long.csv`
(long form, 146 unique genes). The 46 most-changed drivers (ranked by |P6 Wald stat|) are
shown as a within-age z-scored expression heatmap across all P6 and P90 samples
(`figures/metab_driver_genes_zscore_heatmap.png`, matrix in `out/metab_driver_genes_zscore.csv`).
Up-in-KO drivers are OXPHOS/TCA/glycolysis genes (Ndufs7, Sdhb/Sdhd, Uqcrc1, Cyc1, Aldoc,
Eno1/2, Gapdh, Idh3g…); down-in-KO drivers include *Bcat1* itself plus the KEGG "lysine
degradation" set — which is dominated by histone lysine methyltransferases (Kmt2a/c/d, Kmt5a/b,
Ehmt1, Setd7, Setdb1, Prdm2, Nsd2), of potential relevance to the chromatin/H3K9pr angle —
and BCAA/CoA genes (Pank3, Mmut, Acadsb, Aasdhppt). The KO/WT contrast is visually strong at
P6 and largely dissipates at P90 for these P6-defined drivers, consistent with the GSEA.

## Outlier correction — P6 WTr3 removed (`figures/outlier_corrected/`)

P6 sample **WTr3** was a clear outlier by PCA (far from the other WT samples) and in the
driver-gene z-scores. It was inflating the within-WT dispersion, which suppressed DESeq2 power.
Removing it and re-running the entire P6 side:

| P6 analysis | Samples | DEGs (padj<0.05, no LFC filter) | Up / Down |
|-------------|---------|--------------------------------|-----------|
| Original    | 4 WT / 3 KO | 647 | 247 / 400 |
| **WTr3 removed** | 3 WT / 3 KO | **1,887** | 905 / 982 |
| (P90, unchanged) | 4 WT / 4 KO | 25 | 12 / 13 |

Removing the outlier ~tripled the DEG count (647 → 1,887); 606/647 original DEGs are retained
(94% concordant) plus 1,281 newly significant. *Bcat1* unchanged (log2FC −3.05). PCA after
removal: PC1 rises to 40% and cleanly separates KO from WT (`QC_PCA_p6_noWTr3.png`).

**All downstream P6 analyses were re-run on the corrected data** (files/figures suffixed
`_noWTr3`): cell-type GSEA heatmaps (all 6 collections), metabolism ranked plots + KEGG
P6-vs-P90 heatmap, and the driver-gene z-score heatmap (P6 columns now WTr1/2/4 + KO). The
biological conclusions are unchanged but sharper — more sets reach FDR<0.25 and NES magnitudes
increase. The two temporal metabolic patterns (transient P6 energy up-shift; durable/strengthening
BCAA-catabolism deficit) and cell-type findings (ependymal up, oligodendrocyte down, progenitor
age-flip) all hold. Corrected GSEA tables: `gsea/*_p6_noWTr3.csv`; corrected DE table:
`out/deseq2_p6_noWTr3_KOvsWT.csv`.

The original 7-sample P6 outputs are retained alongside for comparison; use the `_noWTr3`
versions as the primary P6 result.

### Braindex cell-type leading-edge z-score heatmap

`figures/outlier_corrected/braindex_celltype_LE_zscore_heatmap.png` — driver genes behind the
Braindex mouse cell-type signals. Starting from Braindex mouse signatures that reach FDR<0.25 in
P6 (corrected) or P90 *and* whose names contain RGC / OPC / GLIP / NFOL / MOL (36 signatures),
their leading-edge genes were pooled into the 5 cell-type groups, deduped within group, and the
top 25 genes/group by magnitude of change are shown as within-age z-scores (P6 columns WTr1/2/4 +
KO; P90 all 8). Gene→signature membership and DE stats: `out/braindex_celltype_LE_genes.csv`.

Read: RGC/OPC/GLIP groups are dominated by the shared **cell-cycle/mitotic program** (Mki67,
Top2a, Ccnb1/Ccnd2, Cenpf, Kif11/20b/23, Nusap1, Hells, Smc4) — down in KO at P6, up at P90 (the
proliferating "-P" age-flip). NFOL and especially **MOL** (myelin genes Fa2h, Ugt8a, Elovl7,
Ermn, Gjc2, Erbb3) are down in KO, more strongly at P90 — the oligodendrocyte-maturation deficit.
Note the RGC/OPC/GLIP gene lists overlap heavily because proliferating progenitors of all three
lineages share the same mitotic leading edge.

**Cell-cycle vs lineage split** (`figures/outlier_corrected/braindex_celltype_LE_zscore_CC_vs_lineage.png`,
table `out/braindex_celltype_LE_CC_vs_lineage.csv`): each cell-type group is split into a
*cell-cycle* sub-block (italic labels) and a *lineage* sub-block. Cell-cycle membership = mouse
Hallmark E2F_TARGETS ∪ G2M_CHECKPOINT ∪ the canonical Tirosh/Seurat S+G2M gene lists. This
separates the two signals: (1) the **shared proliferation program** (Mki67, Top2a, Ccnb1, Cenpf,
Nusap1, Smc4, Kif20b…) is near-identical across RGC/OPC/GLIP "cyc" blocks — down in KO at P6, up at
P90 (age-flip); (2) the **lineage-specific** markers differ by type — OPC-lin (Bcan, Sall3, Podxl2,
Ly6h, Scg5) trends up in KO at P6, whereas NFOL-lin (Gpr17, Enpp6) and MOL-lin (Ugt8a, Ermn, Elovl7,
Erbb3, Tppp) are down in KO, strongest at P90 — the maturation deficit independent of proliferation.
NFOL/MOL have essentially no cell-cycle genes (0 and 1), confirming they are post-mitotic.

### Pathway-grouped metabolic gene heatmap (outlier-corrected)

`figures/outlier_corrected/metab_pathway_grouped_zscore_heatmap.png` (genes/stats in
`out/metab_pathway_grouped_genes.csv`) — corrected-P6/P90 z-score heatmap of DE-relevant
metabolic-pathway members (padj<0.25 at P6), **47 genes across five pathway groups** (colored
sidebar at left). Layout: a genotype annotation bar on top (blue = WT, yellow = KO, matching the
shBCAT1/NT scheme) with legend; a white gap separates the P6 and P90 sample blocks (two distinct
experiments); within-age z-score with the z-scale as a horizontal colorbar below; replicate
labels omitted. A right-hand panel shows the DE **log2FC (KO/WT)** at P6 and P90 per gene on a
yellow(up)/purple(down) scale (clip ±1; Bcat1 saturates). Bold + `*` = padj<0.05 at P6. OXPHOS is
trimmed to its 20 most significant members; the Pantothenate/CoA group was dropped; Bcat1 is
placed in BCAA degradation. First two groups use KEGG sets, the amino-acid groups use curated-panel
sets; genes are deduped across groups.

- **OXPHOS (37 genes)** — coordinately **up in KO at P6** (complex I `Nduf*`, complex III
  `Uqcr*/Cyc1`, complex IV `Cox*`, complex II `Sdhb/Sdhd`, V-ATPase `Atp6v*`); 16 at padj<0.05.
  The up-shift largely dissipates at P90. `Atp6v0a2` is the lone member down in KO.
- **2-Oxocarboxylic acid (9)** — *mixed, not uniformly down*: the TCA/transamination members
  (Got2, Ogdh, Aco2, Cs, Got1, Idh3g/3b, Pdha1) are **up** in KO at P6, while **Bcat1 itself is
  strongly down**. The pathway's negative GSEA NES is therefore driven mainly by the Bcat1
  knockdown (and Gpt) sitting at the bottom of the ranked list — not a coordinated pathway-wide
  decrease. Biologically this reads as loss of the BCAA transaminase with compensatory up-shift
  of the connected TCA/glutamate carbon-handling enzymes.
- **Pantothenate / CoA (7)** — the rate-limiting pantothenate kinases **Pank3 (padj<0.05) and
  Pank2 are down** in KO at P6, along with Aasdhppt, while Coasy/Aldh2/Aldh1b1 are up — so the
  "CoA biosynthesis down" signal is specifically at the pantothenate-kinase entry step.

The three added amino-acid groups (bottom of the heatmap): **BCAA degradation** (Mmut\* down;
Mccc1/Acadsb trend down; Acsf3/Acaa2 up), **Cys/Met metabolism** (Mat2a\* and Amd1\* — the
SAM-synthesis/polyamine branch — down; Mtr/Mat2b down; Mdh2\*/Ldhb/Mpst up), and **One-carbon /
folate** (Tyms\* down, with Shmt1/Mthfd2/Mthfr all trending down — the cleanest coordinately-down
block). The Mat2a/Amd1 + folate down-shift points to reduced SAM / methyl-donor capacity in KO at
P6 (relevant to the H3K9pr/methylation theme). All three effects are P6-specific (P90 columns wash
out). See `out/curated_pathways_leadingedge_p6.csv` and the forest plot
`figures/outlier_corrected/curated_pathways_leadingedge_forest_p6.png`.

Companion schematic: `figures/outlier_corrected/BCAA_2oxo_reaction_schematic_P6.png` draws the
BCAA→2-oxo-acid / TCA reaction network with each enzyme box colored by its corrected-P6 log2FC
(red up / blue down, ±0.6 clip; * = padj<0.05; Bcat1 value printed since it saturates the scale).
It makes the rerouting explicit: **Bcat1 −3.05** at the BCAA transamination entry, the BCKDH
oxidation step and Bcat2 essentially unchanged, and the connected TCA / 2-oxoglutarate carbon-
handling enzymes (Cs, Aco2\*, Idh3, Ogdh\*, Got1/2\*, Pdha1) coordinately **up** — with Gpt
(alanine transaminase) the one down. Values from `out/metab_pathway_grouped_genes.csv`.

## ssGSEA — per-sample cell-type scores, WT vs KO across ages

`figures/ssgsea_braindex_WTvsKO_lineplots.png` — single-sample ssGSEA (gseapy, rank
normalization) of the 290 Braindex mouse signatures across all P6+P90 samples (WTr3 excluded).
Per-sample score matrix: `out/ssgsea_braindex_scores.csv`; group means and change metrics:
`out/ssgsea_braindex_group_stats.csv`.

Line-plot format: y = ssGSEA score; x stratified by age (P6, P90); within each age the WT and
KO group means (±SE, with replicate dots) are connected by a solid line; faint dotted lines
trace each genotype across age. Nine example signatures are shown, chosen for the largest
changes between P6 and P90.

Each panel is annotated with a per-age **KO-vs-WT Welch t-test**, BH-FDR-corrected across all
290 signatures within each age (`out/ssgsea_braindex_KOvsWT_ttests.csv`; stars: *** q<0.001,
** q<0.01, * q<0.05, † q<0.25, ns otherwise). The genotype effect is weak by this test: at P6,
136/290 signatures reach q<0.25 but only 1 reaches q<0.05 (the example progenitor signatures
sit at q≈0.25, †); at P90 nothing survives FDR (n=4v4, effect ~0.05). This is the honest,
underpowered reality of a subtle whole-tissue genotype effect — the preranked GSEA is more
sensitive because it pools information across the ranked gene list rather than testing one
aggregate score with 3–4 replicates.

Key point: the **developmental (age) change dwarfs the genotype effect** — proliferating
progenitor signatures (RGC-P, GLIP-P, ASTP-P) score high at P6 and collapse by P90 (Δ≈−0.4),
while mature types rise (MOL, microglia, excitatory neuron). The KO-vs-WT differences are
smaller (~0.05) but consistent with the GSEA: KO sits *below* WT at P6 for proliferating
progenitors and for MOL/NFOL, *above* WT for OPC (P6) and ependymal (both ages). So BCAT1 loss
nudges cell-type composition on top of a much larger normal maturation trajectory.

## Files

- `figures/heatmap_<collection>.png` — top signatures by NES, P6 vs P90 (`*` = FDR q<0.25)
- `figures/QC_PCA.png`
- `out/deseq2_<age>_KOvsWT.csv` — full DE tables
- `out/rank_<age>.rnk` — ranked lists used for GSEA
- `gsea/<collection>_<age>.csv` — full GSEA output (NES, p, FDR, leading edge)
- `out/GSEA_all_results_long.csv` — all collections/ages combined
- `gmt/*.gmt` — signature sets as used (incl. ortholog-translated human sets)

*Note on environment*: R/DESeq2 were unavailable, so PyDESeq2 (a faithful DESeq2
re-implementation) and gseapy (fgsea-equivalent preranked GSEA) were used. Count
aggregation used summed estimated counts rather than a tximport length-offset; for a
publication-grade run consider re-checking top hits with tximport + DESeq2 in R.
