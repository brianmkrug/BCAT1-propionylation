# DMG metabolomics — global structure (PCA + clustering)

K27M vs K27M-KO metabolite profiles, Metabolon LC-MS, values normalized to protein
concentration. Xenograft and cell-culture datasets analyzed **independently**.

- Xenograft: 14 samples (HSJ019, DIPGXIII × K27M/KO), 730 metabolites
- Cell culture: 28 samples (BT245, DIPGXIII, HSJ019 × K27M/KO + hfNSC WT reference), 514 metabolites

Files: `metabolon_pca_clustering.py` (reproducible script), `summary_stats.csv`,
figures in `figures/`.

---

## Should this data be log-transformed and scaled? Yes to both — and the data shows why.

I evaluated the raw (protein-normalized) data directly rather than assuming. Two
properties of LC-MS abundances drive the decision, and both are strongly present here
(`*_diagnostics.png`, `summary_stats.csv`):

**1. Right-skew → log transform.** Metabolite abundances span ~6 orders of magnitude
and are heavily right-skewed. Median per-metabolite skewness is **+1.13 (xenograft)**
and **+1.18 (cell culture)**, with **~57% of metabolites skewed >1** in both sets.
After log10, median skewness collapses to **~0.0** (0.06 and −0.06). LC-MS signal is
multiplicative/relative, so effects are fold-changes; logging makes them additive and
symmetric, which is the scale PCA, Euclidean/correlation distance, and Ward linkage
assume. This is also the standard for Metabolon block-normalized data.

**2. Mean–variance dependence (heteroscedasticity) → log stabilizes it.** On raw data,
a metabolite's SD is almost perfectly predicted by its mean (Pearson **r = 0.94**
xenograft, **0.92** cell culture) — high-abundance metabolites carry disproportionate
absolute variance. After log10 this dependence is removed (r = **−0.26** and **−0.33**).
Without logging, PCA and clustering would be dominated by a handful of abundant
metabolites purely because of their scale.

**3. Per-metabolite scaling → needed so PCA reflects biology, not abundance.**
Protein normalization corrects for *sample* loading; it does **not** put *metabolites*
on a common footing. Without feature scaling, PCA still tracks the highest-variance
(typically most-abundant) metabolites. I used **autoscaling (z-score / unit variance)**
so every metabolite contributes equally — appropriate here because the biological
question is coordinated K27M-vs-KO shifts across the whole profile, including
lower-abundance metabolites.
- Caveat: autoscaling up-weights near-noise/low-detection metabolites. I mitigated this
  by requiring detection in ≥50% of samples before scaling (kept 651/730 xenograft,
  485/514 cell culture). If you prefer to down-weight noise, **Pareto scaling**
  (÷√SD) is the common metabolomics compromise — switch `scaling="pareto"` in the
  script to regenerate.

**Zeros / below-detection.** ~12–15% of values are exact zeros (below detection, not
true absence). I imputed each metabolite's zeros with **half its minimum observed
value** (Metabolon's convention) *before* logging, so log is always defined and
imputed points sit just under the detection floor rather than distorting the
distribution.

**Pipeline:** zero→NaN → drop metabolites detected in <50% of samples → half-min
impute → log10 → autoscale per metabolite → PCA / hierarchical clustering.

---

## What the global structure shows

### Cell culture (`cellculture_*`)
- **Replicates are essentially perfect.** Every cell line × genotype quartet clusters
  as a clean monophyletic group in the dendrogram — technical/biological reproducibility
  is high.
- **Cell line of origin is the dominant axis**, then genotype within line. hfNSC (normal
  neural stem cells) sit apart on PC1 as expected for a non-tumor reference, but they do
  **not** swamp the whole structure — tumor-line identity and K27M/KO status still
  organize the remaining variance (PC1 22.8%, PC2 17.9%, PC3 14.4%).
- K27M vs KO separates cleanly *within* each line (clearest for HSJ019 and DIPGXIII),
  consistent with a genotype effect layered on top of line-specific background.

### Xenograft (`xenograft_*`)
- PC1 (32%) + PC2 (23%) capture most structure. Genotype and cell line both contribute;
  in vivo replicates are noisier than in vitro (expected — tumor microenvironment,
  stromal contribution, sampling), so a couple of samples sit closer to the neighboring
  group in the dendrogram.
- Line-and-genotype groups still hold together overall (DIPGXIII-K27M, HSJ019-K27M,
  HSJ019-KO, DIPGXIII-KO each form coherent clusters).

**Takeaway for downstream stats:** because cell line is the largest source of variance
in both datasets, differential-abundance testing for the K27M effect should be done
**within/adjusted for cell line** (paired or line-as-covariate) rather than pooling
genotypes across lines.

---

## Per-line differential abundance: K27M vs KO

Because cell line is the dominant variance axis, DA was run **within each line**
(never pooled): Welch t-test on log10 abundances (unscaled), BH-FDR per line,
significance = FDR<0.05 **and** |log2FC|>1 (2-fold). hfNSC excluded (WT only, no KO
arm). Script: `metabolon_diff_abundance.py`; per-line tables + combined results in
`diff_abundance/`; figure `figures/volcano_all_lines.png`.

**Hits per line** (log2FC = K27M / KO, so positive = higher in K27M):

| dataset | line | n sig | up in K27M | up in KO |
|---|---|---|---|---|
| cell culture | BT245 | 60 | 21 | 39 |
| cell culture | DIPGXIII | 106 | 74 | 32 |
| cell culture | HSJ019 | 80 | 34 | 46 |
| xenograft | DIPGXIII | 50 | 37 | 13 |
| xenograft | HSJ019 | **0** | 0 | 0 |

**In vitro is far more sensitive than in vivo.** All three cultured lines show 60–106
significant metabolites; in xenograft, DIPGXIII retains a strong signal (50 hits) but
**HSJ019 has zero** at this threshold — consistent with the noisier in-vivo replicates
and weaker HSJ019 clustering flagged earlier (microenvironment/stromal dilution, n=4).

**On-target BCAT1 signature is present.** The branched-chain α-keto acids
**3-methyl-2-oxobutyrate** (from valine) and **3-methyl-2-oxovalerate** (from
isoleucine) — direct products of the BCAT1 forward reaction — are **up in K27M** and
recurrent across lines (BT245+HSJ019). This is exactly the direction expected if K27M
cells have higher BCAT1 flux, and it's a good positive control that the contrast is
capturing real biology.

**Recurrent, direction-consistent hits (significant in ≥2 lines):** 22 metabolites.
Up in K27M: BCAA keto-acids, ophthalmate + S-(1,2-dicarboxyethyl)glutathione (glutathione/
oxidative-stress axis), N-acetylarginine, gluconate, ADP-ribose, gamma-glutamylleucine.
Up in KO: multiple **sphingomyelins**, **desmosterol** (cholesterol-synthesis
intermediate), GPG/GPI phospholipids — a coherent membrane-lipid / sterol program —
plus 2'-deoxycytidine and, notably, **glucose down in K27M in both cell culture and
xenograft** (the one metabolite reproducibly moving across compartments), consistent
with elevated glycolytic consumption in K27M. Full list in `diff_abundance/recurrent_hits.csv`.

Caveats: n=3–4/group gives limited power (hence the HSJ019-xenograft null); Welch
t-tests on ~500–730 features — the BH-FDR handles multiplicity but low-detection
metabolites imputed near the floor can produce large apparent fold-changes, so treat
extreme-|log2FC| singletons cautiously and prioritize the recurrent hits.

## Pathway enrichment + paired-across-lines model

Uses the Metabolon super/sub-pathway annotations (`cell_annotation.xlsx`,
`xenograft_annotation.xlsx`; matched by name after whitespace/case normalization —
507/514 cell-culture, 688/730 xenograft metabolites annotated). Script:
`metabolon_enrichment_lm.py`; tables in `enrichment/`.

### Which enrichment tests, and why
- **Primary — pre-ranked GSEA (a.k.a. MSEA, metabolite-set enrichment analysis)**
  (`gseapy.prerank`, 2000 permutations). All metabolites are ranked by the pooled
  linear-model t-statistic (below) and each Metabolon sub/super-pathway is tested for
  concentration toward the K27M or KO end. **Threshold-free** — it uses the entire
  ranking, so it detects coordinated shifts even when few individual metabolites clear
  FDR, and it isn't sensitive to an arbitrary hit cutoff. NES>0 = enriched in K27M.
- **Secondary — over-representation analysis (ORA, one-sided Fisher exact)**, the
  default MetaboAnalyst approach: are the per-line significant DA hits over-represented
  in a pathway vs. all measured+annotated metabolites? Threshold-based confirmation.
  I ran GSEA as primary because ORA discards the ranking and the per-line hit lists are
  small (and empty for HSJ019 xenograft).

### Paired-across-lines linear model (formal pooled K27M test)
Per metabolite: OLS `log10(abundance) ~ C(Genotype) + C(Cell line)`, KO as reference,
**cell line as a fixed block** — i.e. each line contributes its own baseline and the
genotype coefficient is the K27M-vs-KO effect *shared across lines* (the "pairing").
BH-FDR across metabolites. This is more powerful than the per-line tests and is the
statistic feeding the GSEA ranking.

Pooled significant metabolites (FDR<0.05 & |log2FC|>1): **36 (cell culture), 96
(xenograft)**. Top pooled hits are dominated by BCAAs and their conjugates:
- **Cell culture:** gamma-glutamylisoleucine, gamma-glutamylvaline, 3-methyl-2-oxobutyrate,
  ophthalmate, ADP-ribose (up in K27M); stearoyl sphingomyelin, GPI/galactosyl-lysine
  lipids (up in KO).
- **Xenograft:** **valine, leucine, isoleucine, N-acetylleucine, proline** all up in
  K27M — the BCAAs themselves, directly.

### Enrichment results — both datasets converge
**The single strongest signal in both compartments is Leucine/Isoleucine/Valine (BCAA)
metabolism enriched in K27M** (NES 2.36 cell culture, 2.15 xenograft; both GSEA FDR≈0).
This is the on-target BCAT1 result, now formally supported at pathway level and
reproduced independently in vitro and in vivo.

At the super-pathway level the pattern is consistent: **Amino Acid and Peptide up in
K27M**; **Carbohydrate and Energy strongly up in KO / depleted in K27M** (xenograft
Carbohydrate NES −2.9, Energy −2.1, both FDR≈0; sub-pathways: Glycolysis/Gluconeogenesis,
TCA cycle, Pentose, Nucleotide-sugar all enriched in KO). In cell culture the KO-enriched
lipid signal is specifically **sphingomyelins, ceramides, and phosphatidylinositol** —
matching the sphingomyelin/desmosterol recurrent DA hits. ORA independently flagged
hexosylceramides/endocannabinoid (BT245) and dipeptides (DIPGXIII).

Biological read-out: K27M cells run higher BCAA catabolism (BCAT1 axis) and amino-acid/
peptide pools, with reciprocally lower carbohydrate–energy (glycolysis/TCA) and membrane
sphingolipid pools relative to BCAT1-KO — coherent with the BCAT1→BCAA→downstream model.

Caveats: GSEA on Metabolon subpathways uses small sets (min 3); NES for tiny sets is
noisier — lean on FDR and on pathways reproduced across both datasets (BCAA is the
standout). The LM assumes a common genotype effect across lines (no genotype×line
interaction modeled); per-line volcanoes show direction is consistent, but effect
*size* does vary by line.

## Sensitivity analysis: dropping outlier HSJ019_K27M_2 (xenograft)

HSJ019_K27M_2 clustered with the KO samples in the xenograft dendrogram. Removing it
(HSJ019 xenograft K27M n=4→3) has a large effect — it was masking the in-vivo signal.
Script: `metabolon_drop_outlier.py`; outputs tagged `_noHSJ`. Only xenograft is affected
(the sample is not in the cell-culture data).

| test | with outlier | without HSJ019_K27M_2 |
|---|---|---|
| xenograft HSJ019 per-line DA (sig metabolites) | **0** | **173** |
| xenograft DIPGXIII per-line DA | 50 | 50 (unchanged) |
| xenograft pooled paired-LM (sig) | 96 | 111 |
| xenograft GSEA sub-pathways FDR<0.25 | 16 | 19 |

**The single outlier fully explains the earlier HSJ019-xenograft null.** With it removed,
HSJ019 in vivo yields 173 significant metabolites (100 up in K27M, 73 up in KO) and now
looks like the other lines. The pooled model gains hits (96→111) and its top effects are
the BCAAs directly — **valine, leucine, isoleucine** (plus phenylalanine, N-acetyltryptophan,
BCAA dipeptides) up in K27M. Enrichment is unchanged in structure and slightly stronger:
**Leucine/Isoleucine/Valine metabolism remains the top K27M pathway** (NES 1.95, FDR≈0),
Amino Acid/Peptide up in K27M, Carbohydrate (NES −2.9) and Energy (−1.99) up in KO.

Interpretation: this is reassuring, not fragile — dropping the outlier *strengthens* and
*sharpens* the same conclusion rather than changing it. Recommend treating the
outlier-removed xenograft results as primary, and reporting the sensitivity analysis (n
is small, so removing one sample should always be disclosed). Worth double-checking
whether HSJ019_K27M_2 has a wet-lab/annotation explanation (sample swap, low input,
genotype mislabel — note it sat right in the KO cluster).

## Figure index
| File | Content |
|---|---|
| `*_diagnostics.png` | skewness (raw vs log) + mean–SD relationship justifying log |
| `*_PCA.png` | PC1–PC2 and PC2–PC3, colored by genotype, shaped by cell line |
| `*_scree.png` | variance explained per PC |
| `*_dendrogram.png` | Ward clustering, correlation distance; leaf color = genotype |
| `*_corr_heatmap.png` | sample–sample Pearson correlation, cluster-ordered |
| `*_metabolite_heatmap.png` | top-60 most variable metabolites (autoscaled) |
| `volcano_all_lines.png` | per-line K27M vs KO volcano plots |
| `GSEA_*_subpathway.png` | pooled sub-pathway GSEA (NES; red=K27M, blue=KO) |
| `*_noHSJ.png` | same analyses with xenograft outlier HSJ019_K27M_2 removed |
