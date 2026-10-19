# BCAT1 vs BCAT2 — sequence alignment, predicted structure/disorder, vertebrate conservation

Generated 2026-09-29. All source data pulled live from UniProt, AlphaFold DB and the
UCSC Genome Browser API; nothing here is hand-entered.

## Figures

| file | contents |
|---|---|
| `Fig1_BCAT1_BCAT2_alignment_overview.pdf/.png` | (a) domain schematic of both proteins on alignment coordinates, (b) per-column identity/similarity, (c) 15-column sliding identity |
| `Fig1b_BCAT1_BCAT2_sequence_alignment.pdf/.png` | full residue-level alignment, 60 columns/block, identity + similarity shading, PLP lysine and CXXC marked |
| `Fig2_BCAT_structure_disorder.pdf/.png` | (a–c) AlphaFold cartoons coloured by pLDDT + superposition, (d–e) 2D per-residue tracks: DSSP secondary structure, pLDDT, metapredict V3 disorder, FoldIndex |
| `Fig3_BCAT_conservation_UCSC100way.pdf/.png` | (a–d) phyloP and phastCons over the complete CDS of each gene with the CDS exon map, (e) N-term vs core, (f) conservation by codon position, (g) both genes on shared alignment coordinates |

PDFs are vector with editable text (`pdf.fonttype 42`); PNGs are 450 dpi.
The 3-D cartoons in Fig 2a–c are raster (3Dmol.js, 4200 px) — everything else is vector.

## Key numbers

**Alignment** (Needleman–Wunsch, BLOSUM62, gap open −11 / extend −1, free end gaps)

- 377 aligned positions; **217 identical (57.6%)**, +61 similar → **73.7% similarity**
- BCAT1 residues **1–25** have no BCAT2 counterpart: the first 15 residues align only
  weakly, then BCAT1 carries a 14-column gap opposite the BCAT2 mitochondrial transit
  peptide. Windowed identity does not exceed 40% until alignment column 57.
- Catalytic PLP lysine aligns: BCAT1 K222 ↔ BCAT2 K229 (alignment column 236).
- CXXC redox motif aligns: BCAT1 C335-V-V-C338 ↔ BCAT2 C342-Q-V-C345 (column 349).

**Structure / disorder**

| | BCAT1 (P54687) | BCAT2 (O15382) |
|---|---|---|
| AlphaFold model | AF-P54687-F1 v6 | AF-O15382-F1 v6 |
| length | 386 aa | 392 aa |
| mean pLDDT | 94.4 | 93.1 |
| residues pLDDT < 70 | 23 (N-terminal 1–19) | 31 (N-terminal 1–29) |
| metapredict V3 IDR | **1–25** | **1–36** |
| DSSP H/E/coil | 0.24 / 0.28 / 0.47 | 0.22 / 0.28 / 0.49 |

Superposition of the two AlphaFold models (biotite `superimpose_homologs`, iterative
outlier removal): **Cα RMSD 0.35 Å over 235 equivalent residues** — the folded cores are
essentially identical. The only structural difference is the extended, low-confidence
N-terminal tail, which is present in both but sequence-divergent.

**Conservation** (UCSC hg38, 100-way vertebrate phyloP and phastCons, base resolution)

| | BCAT1 | BCAT2 |
|---|---|---|
| transcript | NM_005504.7 / ENST00000261192.12 (MANE Select) | NM_001190.4 / ENST00000316273.11 (MANE Select) |
| locus | chr12, − strand, 11 coding exons | chr19, − strand, 11 coding exons |
| coding bases | 1,161 (incl. stop) | 1,179 (incl. stop) |
| mean phyloP, whole CDS | 3.21 | 2.55 |
| mean phyloP, N-term vs core | **0.51 vs 3.39** (Mann–Whitney p = 1.3 × 10⁻¹¹) | **1.20 vs 2.69** (p = 3.6 × 10⁻⁶) |
| mean phastCons, N-term vs core | 0.31 vs 0.80 | 0.64 vs 0.69 |
| phyloP by codon position (1st/2nd/3rd) | 3.89 / 5.26 / 0.47 | 3.21 / 4.39 / 0.05 |

The BCAT1 N-terminal extension is therefore not only absent from BCAT2 — it is the least
constrained part of the BCAT1 coding sequence across vertebrates, while the catalytic core
is under strong purifying selection.

## Verification performed

- UCSC CDS, spliced in transcript orientation and reverse-complemented (both genes are on
  the minus strand), translates **exactly** to the UniProt canonical protein for both genes.
- CDS length = 3 × (protein length + 1) for both (stop codon included).
- 100% base coverage of phyloP and phastCons over every coding base (1,161 / 1,179).
- AlphaFold model sequence identical to the UniProt sequence for both entries.
- Codon-position signature (1st/2nd ≫ 3rd) confirms the reading frame and strand mapping.
- PLP lysine and CXXC positions re-derived from sequence, not assumed.

## Data tables

`BCAT1_per_residue_table.tsv`, `BCAT2_per_residue_table.tsv` — one row per residue:
DSSP state, pLDDT, relative solvent accessibility, metapredict V3 disorder, FoldIndex,
per-codon phyloP (mean and each codon position), per-codon phastCons, and the aligned
partner residue in the other paralogue with an identity flag.

`conservation_summary.json` — the summary statistics plotted in Fig 3e–f.

## Scripts (`scripts/`)

Run in this order; each writes into `data/`:

1. `fetch_cons.py` — MANE CDS structure + phyloP/phastCons from the UCSC REST API, with the translation check.
2. `struct_features.py` — pLDDT from the AlphaFold PDB, DSSP secondary structure and RSA, FoldIndex.
3. `superpose.py` — superimposes BCAT2 onto BCAT1 (biotite).
4. `render3d.py` — 3Dmol.js cartoons rendered headless (playwright + chromium).
5. `fig01_alignment.py`, `fig01b_alignment_sequence.py`, `fig02_structure.py`, `fig03_conservation.py`.

Colours follow `BCAT1_figure_palette.json` (`genes`: BCAT1 `#1b9e77`, BCAT2 `#7570b3`).

## Sources

- UniProt P54687 (BCAT1, cytosolic) and O15382 (BCAT2, mitochondrial).
- AlphaFold Protein Structure Database, model v6 (2025-08-01).
- UCSC Genome Browser REST API, hg38: `mane`, `phyloP100way`, `phastCons100way`
  (100-way vertebrate multiz alignment; track data dated 2015-05-08).
- metapredict V3 (Emenecker/Holehouse) for disorder; FoldIndex (Prilusky et al. 2005)
  as an orthogonal sequence-based check; DSSP for secondary structure.

## Caveats

- phyloP/phastCons are **nucleotide-level** scores from the genomic alignment. Third-codon
  positions are largely unconstrained, so per-codon means are dominated by the 1st/2nd
  positions — this is why Fig 3f is included. They are not a substitute for a protein-level
  multi-species alignment, which would be the next step if you want to show which individual
  residues are invariant across vertebrates.
- The AlphaFold N-terminal tails are low-confidence predictions (pLDDT < 70). Their extended
  conformation in the model should be read as "no confident structure", not as evidence of a
  genuinely extended rod.
- BCAT2's N-terminal 1–27 is a cleaved mitochondrial transit peptide, so part of the
  N-terminal divergence between the paralogues is targeting sequence rather than a retained
  domain. The BCAT1 1–25 extension has no such annotation.

---

## Addendum (2026-09-29): NLS prediction

`Fig4_BCAT_NLS_prediction.pdf/.png` and `NLS_prediction_notes.md` test whether the
BCAT1 N-terminal extension could act as a nuclear localization signal.

Short answer: no. cNLS Mapper finds no monopartite NLS anywhere in BCAT1 and only
score-2.0–2.6 bipartite hits (Kosugi's "cytoplasmic" band), none in residues 1–25;
PSORT II reports pat4/pat7/bipartite all absent with NLS score −0.47; NucPred scores
BCAT1 0.02 against a 0.8 threshold, below inert cytosolic negative controls. BCAT1
1–25 is net −2 and 8% K+R — acidic and cysteine-rich, not basic.

BCAT2 differs: cNLS Mapper scores a bipartite NLS at residue 25 at 5.8, the Human
Protein Atlas annotates BCAT2 as nucleoplasm + mitochondria, and a 2001 paper
(PMID 11574535) reported a BCATm splice variant acting as a nuclear co-repressor of
thyroid hormone receptor. See the notes file for the caveats.

`data/predictor_results.json` holds the verbatim server outputs;
`scripts/nls_scan.py` is the local pattern scan with its positive/negative controls.

---

## Composite publication figure (2026-09-29)

`Fig_composite_BCAT1_BCAT2.pdf` / `.svg` / `.png` — four panels in one figure,
built by `scripts/fig_composite.py` from the same underlying data:

- **A** residue-level BCAT1/BCAT2 alignment (was Fig 1b)
- **B** AlphaFold superposition (was Fig 2c)
- **C** BCAT1 phastCons across the complete CDS (was Fig 3b)
- **D** BCAT2 phastCons across the complete CDS (was Fig 3d)

7.2 × 10.3 in. PDF and SVG are vector with live text (`pdf.fonttype 42`,
`svg.fonttype none`) — every label is editable in Illustrator or Inkscape. The only
raster element is the 3-D cartoon in panel B (4200 px, embedded). PNG is a 450 dpi proof.

Captions are set inside the figure. Copy-paste version for a manuscript:

> **Figure X. BCAT1 and BCAT2 share a conserved catalytic core but differ entirely at the N-terminus.**
> **(A)** Global pairwise alignment of human BCAT1 (UniProt P54687, 386 aa) and BCAT2 (O15382, 392 aa),
> BCAT1 as reference. Needleman–Wunsch with BLOSUM62, gap open −11 / extend −1 and free end gaps:
> 217 of 377 aligned positions are identical (57.6%), 73.7% similar. Dark boxes, identical residues;
> light boxes, similar (BLOSUM62 > 0). Tinted residues mark the BCAT1 N-terminal extension (1–25) and
> the BCAT2 transit-peptide region (1–36); neither has a counterpart in the other paralogue, and
> windowed identity does not exceed 40% until alignment column 57. The catalytic PLP-binding lysine
> (BCAT1 K222 / BCAT2 K229) and the CXXC redox motif (BCAT1 C335–C338 / BCAT2 C342–C345) are aligned.
> **(B)** AlphaFold DB models AF-P54687-F1 and AF-O15382-F1 (v6) superimposed on structurally
> equivalent Cα atoms (iterative outlier removal): RMSD 0.35 Å over 235 residues. Folded cores, grey;
> BCAT1 1–25, green; BCAT2 1–30, purple. Both N-terminal tails are low-confidence predictions
> (pLDDT < 70) and are called intrinsically disordered (metapredict V3: BCAT1 1–25, BCAT2 1–36).
> **(C, D)** phastCons conservation (UCSC hg38, 100-way vertebrate alignment) across the complete
> coding sequence of each gene, shown as the per-codon mean of the three base scores against residue
> number. Dashed lines and the block track give CDS exon boundaries; shading marks the divergent
> N-terminal region. Mean phastCons in the N-terminal region versus the folded core is 0.31 vs 0.80
> for BCAT1 and 0.64 vs 0.69 for BCAT2; the corresponding phyloP difference is significant for both
> (Mann–Whitney p = 1.3 × 10⁻¹¹ and 3.6 × 10⁻⁶). The spliced CDS translates exactly to the UniProt
> canonical sequence for both transcripts, and phastCons covers 100% of coding bases.
