# BCAT1 / histone propionylation in H3K27M glioma — analysis and figure code

Code accompanying:

> Krug B, Gatesman TA, Lopez Gutierrez R, *et al.*, Agnihotri S, Jabado N, Lu C.
> **Co-option of developmental branched-chain amino acid catabolism maintains glioma progenitor cell states through histone propionylation.**
> *Manuscript submitted, 2026.* (Citation and DOI will be added on publication.)

This repository contains the scripts that generate the computational panels of Figures 1–7 and
Supplementary Figures S1–S7, plus the upstream analysis steps they depend on (single-cell
integration and pseudotime, CUT&RUN normalization, embryonic atlas annotation, protein
structure/conservation). It does **not** contain data; see [`DATA_AVAILABILITY.md`](DATA_AVAILABILITY.md).

## Where to start

| File | What it is |
|---|---|
| [`FIGURE_INDEX.md`](FIGURE_INDEX.md) | Every figure panel → the script that draws it (also as [`FIGURE_INDEX.tsv`](FIGURE_INDEX.tsv)) |
| [`PROVENANCE.md`](PROVENANCE.md) | Which scripts are the originals, which were rebuilt and checked against the published panel, and which panels have no code |
| [`DATA_AVAILABILITY.md`](DATA_AVAILABILITY.md) | Every input dataset: our GEO deposits, public datasets reused, and data available on request |
| [`requirements.txt`](requirements.txt) / [`environment.yml`](environment.yml) | Python and R dependencies |
| `BCAT1_figure_palette.json` | Shared colour palette; scripts locate the repository root by this file |

## Repository layout

Directories mirror the analysis folders used during the project, so the relative paths inside
each script are unchanged.

| Directory | Analysis | Figures |
|---|---|---|
| `scXenograft/scripts/` | HSJ019 K27M vs K27M-KO xenograft scRNA-seq: QC, Harmony integration, diffusion pseudotime, cell-type labelling, gene heatmaps along pseudotime | 1A–B, S1A, S2E, S4G |
| `metabolon data/analysis/` | Metabolon profiling of isogenic cell cultures and xenografts | 1D–G, S1B–F, S3D, S3F |
| `jessa2022_metabolic_ssgsea/scripts/` | Lipid-gene expression across pediatric brain tumour scRNA-seq (Jessa 2022, Riemondy 2022) | S1G–H |
| `CRISPR screen/` | Metabolism-focused CRISPR screen ranking, KEGG enrichment, DepMap classification | 2B, S2A |
| `human_data/` | Human OPC→oligodendrocyte pseudotime (Velmeshev 2023) | 2E |
| `post_gastrulation/scripts_post_gastrulation/` | BCAT1/BCAT2 in post-gastrulation human embryos (Zeng 2023, Xu 2023) | S2F–G |
| `xHSJ019_K27M_vs_KO_analysis/` | Xenograft bulk RNA-seq metabolic gene heatmap | S2D |
| `CAGE/` | CAGE-seq promoter motif analysis | S2K |
| `plasmidsaurusRNAseq/` | shBCAT1 and isoleucine-restriction RNA-seq: GSEA and cross-contrast concordance | 3B–C, S3A, S3I, 7I–J, S7F |
| `BCAT1_KO_analysis/` | *Bcat1*-KO mouse brain RNA-seq (P6, P90) | 3F–G |
| `LFB/` | Luxol fast blue optical-density quantification | 3H |
| `manuscript bar plots/` | Quantification of immunoblots, clonogenic assays, α-KG assay, ENL immunofluorescence | S2B, S3B, S3E, S3G–H, S5D, S6B–E, S7E |
| `isotope/` | Acyl-CoA abundance and ¹³C-BCAA isotope tracing | 4C–E, S4A |
| `MDT patient data/` | Patient tumour / normal brain tissue metabolomics | S4C–E |
| `histoneMS/` | Histone PTM mass spectrometry and proteomics | 5B, S5J |
| `BCAT_protein_structure_conservation/scripts/` | BCAT1/BCAT2 alignment, AlphaFold superposition, phastCons | S5K–M |
| `CUTRUN/scripts/` | CUT&RUN counting, spike-in/IgG normalization, promoter and genome-wide analyses | 6A–B, 6G–I, S6A, S6F–G, 7F |
| `BCAAr_analysis/` | BCAA-restricted vs control diet allograft tumour RNA-seq | 7D–E, S7B |
| `diet data/`, `media/` | Diet composition, survival, medium vs plasma BCAA | 7A–B, 7G, S7A |
| `scRNA_mouse/gsea_analysis/scripts/` | Mouse H3K27M DMG allograft scRNA-seq (McNicholas 2023) | S7D |

Schematics (BioRender) and wet-lab images (immunoblots, IHC/IF, plates, histology) have no code.

## Running a script

1. Create the environment: `conda env create -f environment.yml && conda activate bcat1-figures`
   (or `pip install -r requirements.txt`; R packages are listed in `environment.yml`).
2. Place the input data at the paths the script reads (each script's docstring or first lines name
   its inputs; [`DATA_AVAILABILITY.md`](DATA_AVAILABILITY.md) says where each dataset comes from).
3. Run from the repository root, e.g.

   ```bash
   python "manuscript bar plots/plot_Clonogenic.py"
   python isotope/plot_propionylCoA_stacked.py
   ```

Scripts resolve paths relative to the repository root. Scripts that previously hard-coded an
absolute path now begin with a one-line `REPO_ROOT` definition; set `BCAT1_REPO=/path/to/repo` to
override it. R scripts use `Sys.getenv("BCAT1_REPO", ".")`, so run them from the repository root.
Outputs are written into the `figures/`, `plots/` or `out/` folder of each analysis directory.

The `post_gastrulation` scripts were run from inside `post_gastrulation/` with intermediate files in
`post_gastrulation/work/`; run them from that directory.

## Statistics

The statistical test for each panel is stated in the figure legends and repeated in each script's
docstring. In brief: Welch's two-sided t-test on log-transformed abundances (metabolomics, acyl-CoA,
immunoblot quantification), paired or ratio-paired t-tests where replicates are matched, DESeq2 /
PyDESeq2 with Benjamini–Hochberg adjustment (RNA-seq), pre-ranked GSEA (gseapy) with BH FDR,
Holm-corrected Welch tests on per-animal means (LFB, 3H), log-rank tests (survival), a one-sided
Mann–Whitney U test (6I) and two-sided Fisher's exact tests (super-ENL promoter sets, S6F).

## License

Code is released under the MIT License ([`LICENSE`](LICENSE)). Third-party datasets are subject to
their own terms.

## Contact

Brian Krug — corresponding author. Please open a GitHub issue for questions about the code.
