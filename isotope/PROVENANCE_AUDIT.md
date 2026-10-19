# `isotope/` — provenance audit (2026-08-28)

Single input file for all five panels:
`Copy of 240911_LuChao_Brian_BCAA_acylCoA_Trace_Peak Area Summary copy.xlsx`
No scripts, no intermediate CSVs, no notes exist in this folder.

## Workbook structure (14 sheets)

| sheet | role |
|---|---|
| `Sheet1_1` | **raw instrument export.** 225 samples x 70 isotopologue AUC columns (M0, 13C1…13C8 per acyl-CoA). Conditions: K27M-shNT, K27M-shBCAT1, K27M/KO in DM and HPLM, K27M/KO. Tracers: Leucine-C13, Isoleucine-C13, Valine-C13, Unlabeled. |
| `shBCAT1`, `HPLM`, `SCM` | per-experiment isotopologue percentages, natural-abundance-corrected (negative values present) |
| `isotopologue sums`, `Sheet2`, `Sheet1` | intermediate summations / paired peak-area blocks |
| `K27M-KO molar enrichments` | per-replicate molar enrichment %, 8 acyl-CoAs x tracer x genotype |
| `K27M-KO abundance labeling` | per-replicate propionyl-CoA abundance + label |
| `Plotting shBCAT1`, `Plot HPLM`, `Plot SCM`, `Plot log2FC`, `Plot K27m-KO abundance label` | tidy, plot-ready tables |

**Every derived sheet is hard-typed. There is not a single formula anywhere in the
workbook.** The derivation from raw AUC to molar enrichment was done by hand or pasted
as values, so the computation is not recorded in machine-readable form.

## What was verified

Reconstructed the plotting + statistics from the tidy sheets and compared with the
values printed on the figures.

| Panel | Source sheet | Result |
|---|---|---|
| **Fig 4C** molar enrichment | `K27M-KO molar enrichments` | **exact.** Welch t per metabolite, Benjamini–Hochberg across the 7 metabolites within each cell line x tracer. BT245/Ile: propionyl q=3.75e-08 (fig 4e-8), succinyl q=3.97e-04 (fig 4e-4). DIPGXIII/Ile: 3HMG 6.41e-03 (fig 0.0064), acetyl 6.41e-03 (fig 0.0064), propionyl 1.90e-03 (fig 0.0019), succinyl 7.29e-03 (fig 0.0073). |
| **Fig 4D** sh-NT vs sh-BCAT1 | `Plotting shBCAT1` | **exact.** Bars 1.00 / 0.15, label 0.24 / 0.02. Welch p = 0.0022 (sum) and 0.0013 (label) — both as printed. |
| **Fig 4E** SCM +/- | `Plot SCM` | **exact.** Bars 1.00 / 0.36, label 0.23 / 0.07. Welch p = 0.0023 and 0.0015 — both as printed. |
| **Fig S4B** HPLM K27M vs KO | `Plot HPLM`, **HPLM rows only** | **exact.** Bars 1.00 / 0.50, label 0.09 / 0.02. Welch p = 0.0040 and 0.0078 — both as printed. (The DMEM-F12 rows in the same sheet give 1.00/0.22 and p=0.0007 — not this panel.) |
| **Fig 4B** log2FC CoA species | `Plot log2FC` | plot reproducible from the tidy sheet, but the log2FC and p values **cannot be recomputed**: `Sheet2` has paired unlabelled numeric columns with no key and contains no KO rows. |

## Recovered recipe (for STAR Methods)

Panels 4D, 4E, S4B: each replicate's labelled and unlabelled propionyl-CoA is expressed
relative to the **mean total (labelled + unlabelled) of the control group** in that
comparison. Bars are group means of that normalized total; the black sub-bar is the mean
normalized labelled fraction; error bars are SEM; points are individual replicates
(n = 5 per group throughout). Two **Welch unequal-variance two-tailed t-tests** are
reported per comparison — one on the summed pool, one on the labelled fraction — and are
**not** corrected for multiple comparisons.

Panel 4C: Welch unequal-variance two-tailed t-test per acyl-CoA, Benjamini–Hochberg FDR
across the seven quantified metabolites within each cell line x tracer; q printed where
< 0.05. Note DIPGXIII/isoleucine butyryl-CoA reaches q = 0.033 but is not labelled on
the figure — confirm whether that is intentional.

## Gaps that remain

1. **Natural-abundance correction.** The corrected percentages appear already in
   `shBCAT1`/`HPLM`/`SCM`; the algorithm and software used are not recorded. Needed for
   Methods.
2. **Molar enrichment formula.** Not recorded. The values in
   `K27M-KO molar enrichments` cannot be re-derived from `Sheet1_1` without knowing it.
3. **Fig 4B derivation.** log2FC(K27M/KO) and its p value are not reconstructable from
   any sheet present. The underlying K27M-vs-KO abundance comparison may live in a file
   not in this folder.
4. **`SCM_propionylCoA_stacked` vs the "SCF(+)/SCF(-)" axis label** in Figure 4 — the
   workbook sheet is named SCM and reads "Stem cell factors (+)/(-)". The figure label
   should probably be SCM.

## Cell Press position

Cell Press requires that **all original code** be deposited with a DOI, and offers
"This paper does not report original code" as a valid section-2 statement. An analysis
done in Excel generates no original code, so its absence is not a policy breach. The
binding obligations for this folder are (a) the data — this workbook, or a tidy
derivative — deposited or available from the lead contact, and (b) a Quantification and
Statistical Analysis description of how enrichment was computed, which currently cannot
be written from the file alone.
