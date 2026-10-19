# post_gastrulation — BCAT1 in post-gastrulation human embryos

## Downloads (2026-10-09)
- xu2023_GSE157329/ — Xu et al. 2023 Nat Cell Biol, CS12–16 whole embryo, 10x. 185,140 cells, 32,351 genes.
  raw_counts.mtx.gz (genes x cells, integer UMIs) + cell_annotate (developmental system / final_annotation / embryo / stage / dissection_part).
- zeng2023_GSE155121/ — Zeng et al. 2023 Cell Stem Cell, PCW3–12, 10x. GSE155121_human_data_raw.h5ad.gz = 463,304 cells x 32,738 genes (hg19-style names), raw integer UMIs.
  obs carries ONLY `week_stage` (sample: W3-1, W4-1..3, W5-1..3, W6-1, W7-1, W8-1, W9-1..2, W12-1). NO cell-type labels.
  Sample tissue from GEO titles: W3 & W4-1 & W5-1 whole body; W4-2/3, W5-2/3, W6–W8 whole head; W9, W12 whole brain.
  GSM4695826 PCW3 matrix is the same W3 sample (redundant). Decompressed h5ad is 10 GB; the .h5ad (0 bytes) is a truncated placeholder from an interrupted decompress — ignore/delete.

## First pass (Xu, computed by streaming the mtx)
pbCPM = pseudobulk counts per million (sum of gene UMIs / sum of all UMIs in group). See xu2023_BCAT1_by_*.csv.

## Zeng 2023 annotation (2026-10-09; computed in cloud workspace, scripts in scripts_post_gastrulation/)
QC: n_genes>=500, pct_mt<=20, total UMI <= 99.5th pct per sample -> 446,689 of 463,304 cells.
Gene set: 2,768 genes (2,500 Zeng HVGs + Xu cell-type markers + manual lineage markers), matched by Ensembl ID.
1) Marker-based: PCA (50, fit on 120k subsample) -> Harmony by sample -> kNN(15) -> Leiden res 1.0 (40 clusters),
   hand-named from top cluster genes + lineage marker detection -> marker_fine / marker_coarse.
   OPC_like = PDGFRA>0 & OLIG2>0 within neural-progenitor clusters (brain OPC_like: OLIG1 96%, NKX2-2 80%, SOX10 58%, AQP4 1%).
   W4-W5 OPC_like cells are few and probably not true OPCs - interpret brain (W9-12) only.
2) Label transfer: PCA basis fit on Xu (each dataset z-scored with its own stats), Harmony by dataset, 30-NN vote from Xu
   -> xu_developmental_system / xu_final_annotation with conf_* = vote fraction. Xu leave-one-out system accuracy 0.85.
   Agreement with marker labels: neural progenitor 90%, neuron 87%, neural crest->schwann 95%, endothelium 98%, erythroid 96%.
   Fine-type (333 Xu types) confidence is low (median 0.23-0.53), and Xu has no glia/cortical neurons, so W9-12 fine
   transfers are not trustworthy; use marker_fine there.
Files: zeng2023_cell_annotations.csv.gz (per cell, incl. raw BCAT1/BCAT2/MYC UMIs and total UMI), zeng2023_BCAT1_*.csv, zeng2023_BCAT1_overview.png/pdf.
Caveats: one embryo per sample, sample = week = tissue confounded; W9-12 are brain-only and lower depth (median ~5k UMI); HBA/HBG ambient contamination is widespread.

## BCAT2 comparison (2026-10-09)
Same groupings, BCAT1 and BCAT2 side by side (+ BCAT1/BCAT2 pbCPM ratio): *_BCAT1_BCAT2_*.csv, figure zeng2023_BCAT1_vs_BCAT2.png/pdf (scripts_post_gastrulation/fig_bcat12.py).
