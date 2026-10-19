# Provenance of the figure code

Each plotted panel falls into one of four categories. The category for every panel is in
[`FIGURE_INDEX.md`](FIGURE_INDEX.md).

| Category | Meaning | Panels |
|---|---|---|
| **original** | The script that produced the file placed in the published figure. | 1A, 1B, 1D, 1E, 1F, 1G, 2B, 2E, 3B, 3C, 3H, 6A, 7D, 7E, 7I, 7J, S1A, S1B, S1C, S1D, S1E, S1F, S1G, S1H, S2A, S2B, S2D, S2E, S2F, S2G, S2K, S3A, S3B, S3D, S3E, S3F, S3G, S3H, S3I, S4G, S5D, S5K, S5L, S5M, S6B, S6C, S6D, S6E, S7B, S7D, S7E, S7F |
| **rebuilt ✓** | The original script was not preserved. The script here was rewritten from the source data and its output was checked value-by-value against the statistics printed on the published panel (p/q values, n, bar heights, correlation coefficients). Styling may differ slightly from the published panel; the numbers do not. | 3F, 3G, 4C, 4D, 4E, 5B, 6B, 6G, 7A, 7B, 7G, S4A, S4C, S4D, S4E, S5J, S6A, S6F, S7A |
| **rebuilt ✗** | A rewritten script whose output does **not** reproduce the published panel with the data currently available. Kept for transparency. | 6H, 6I, 7F, S6G |
| **no script** | The panel was made in an interactive session whose code was not kept, or the output file could not be traced. | 2F, 2G, 4B, 5A, 5E, 6J, S2H, S2I, S2J, S2L, S3C, S5A, S5G, S5H, S6H, S7C |

## Rebuilt ✗ — details

- **6H, 6I, S6G** (`CUTRUN/scripts/plot_delta_promoter_panels.py`). With the promoter delta table in
  hand, 6H gives 40 genes in the both-down quadrant (panel: 130) and 6I gives p = 1.8 × 10⁻⁴ (panel:
  2 × 10⁻⁷). The published panels were evidently made from a different promoter universe.
- **7F** (`CUTRUN/scripts/plot_mm10_FAmetab_H3K9pr.py`). Not yet compared against the published panel.

## No script — details

- **2F** — BCAT1/BCAT2 expression by glioma subgroup (boxplots).
- **2G** — Projected cell-type composition by glioma subgroup.
- **4B** — Acyl-CoA log2 fold changes (log2 ratio of mean total pools, Student's t-test); computed in a spreadsheet, no script.
- **5A** — BT245 histone acyl-PTM stoichiometry.
- **5E** — BT245 acyl-PTM changes on BCAT1 knockdown.
- **6J** — CUT&RUN genome tracks at HSD17B11.
- **S2H** — BCAT1 expression and detection across fine progenitor populations (Mannens et al.).
- **S2I** — BCAA-enzyme expression across three human brain references.
- **S2J** — BCAT1/BCAT2 gene detection by glioma subgroup.
- **S2L** — CAGE signal and motif positions at the BCAT1 promoter.
- **S3C** — P6 *Bcat1*-KO metabolic-pathway gene heatmap (gene list provided in `BCAT1_KO_analysis`).
- **S5A** — DIPGXIII histone acyl-PTM stoichiometry.
- **S5G** — DIPGXIII acyl-PTM changes on BCAT1 knockdown.
- **S5H** — H3K9 methylation states on BCAT1 knockdown.
- **S6H** — Promoter ΔENL vs ΔH3K9pr with SREBP-activation genes.
- **S7C** — H3K9pr change by genomic feature in BCAAr tumours.

All other panels without code are BioRender schematics or wet-lab images.

## Ambiguous attributions

Two or more scripts write the same output file with identical timestamps for **1E** and **S1G**; all
candidates are included and the index names the most likely one.

## Path changes made for publication

Scripts originally written in cloud analysis sessions hard-coded absolute paths
(`/sessions/<id>/mnt/BCAT1`, `/Users/.../BCAT1`). In this repository those prefixes were replaced by
`REPO_ROOT`, defined in one line at the top of each affected script. No other line was changed.
