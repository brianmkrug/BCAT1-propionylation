# Quantifying Luxol Fast Blue myelin signal in QuPath
### BCAT1 WT vs KO mouse brain — Klüver-Barrera (LFB + cresyl violet), whole-slide scans, manual regional ROIs

---

## 0. Short answer

Yes. QuPath (free, open-source — there is no paid version) has a genuinely robust brightfield
quantification path built on **Ruifrok–Johnston color deconvolution**, which is the correct
math for absorbing dyes like LFB. Klüver-Barrera is close to an ideal case: two chromogens
(LFB = blue, cresyl violet = purple), which is exactly what 2-stain deconvolution is designed
for. The workflow below separates LFB into its own optical-density channel, thresholds it with
a **pixel classifier** that is applied identically to every slide, and reports per-region
metrics.

**The single most important point:** LFB is *semi-quantitative*. Signal depends on section
thickness, staining batch, and especially the lithium-carbonate differentiation step. QuPath
will happily give you 4 decimal places of a number that is mostly a staining artifact if the
upstream histology isn't controlled. Section 1 matters more than sections 2–7.

---

## 1. Pre-analysis requirements (do not skip)

| Requirement | Why |
|---|---|
| Stain **all** WT and KO sections in a **single batch**, same jar, same day | LFB differentiation time is the dominant source of variance; it swamps genotype effects |
| Put WT and KO sections **on the same slide** where possible | Removes slide-level staining as a confounder entirely |
| Match rostrocaudal level between genotypes (record bregma coordinates) | Corpus callosum thickness and cortical WM density change steeply along the A–P axis |
| Uniform section thickness from the same microtome/cryostat session | OD scales with path length |
| Scan **all slides in one session**, fixed exposure/white balance, **auto-exposure OFF** | Deconvolution background and OD are meaningless across white-balance changes |
| Rename files to blinded codes before analysis; unblind only at the stats step | ROI drawing and threshold setting are subjective |
| ≥3 sections per region per animal, ≥4–5 animals per genotype | n = animal, not section (see §7) |

Also include an empty-glass area in the scan region — you need it to set the background RGB.

---

## 2. Project setup

1. `File → Project → Create project` → point at an empty folder (e.g. `BCAT1/LFB_QuPath/project/`).
2. `File → Project → Add images` → drag in all `.svs` / `.ndpi` / `.czi` / `.mrxs` files.
   - **Set image type: `Brightfield (other)`** — *not* `Brightfield (H&E)`. H&E presets will
     silently apply hematoxylin/eosin vectors that are wrong for LFB.
   - Leave "Auto-generate pyramids" as default; check **Import objects** off.
3. Confirm pixel size is read correctly (`Image` tab → Pixel width/height in µm). If blank, set
   it manually — every downstream µm-based setting depends on it.

---

## 3. Define the LFB and cresyl violet stain vectors (do this **once**)

You want *one* set of stain vectors applied to *every* image in the study. Per-image
auto-estimation is the classic mistake: it makes each slide's LFB channel mean something
slightly different, which biases group comparisons.

1. Open one representative, well-stained image (ideally a WT, mid-range staining intensity).
2. Draw a **small rectangle** over a region containing: strongly LFB+ white matter, clearly
   CV+ / LFB− grey matter neurons, and a patch of blank glass. Keep it small — large regions get
   downsampled and dilute the color information.
3. `Analyze → Preprocessing → Estimate stain vectors`.
4. In the dialog:
   - Press **Auto** to let QuPath fit the vectors from the pixel cloud.
   - **Leave "Exclude unrecognised colors (HE only)" UNCHECKED** — that option throws away
     anything not pink/blue-ish and will mangle LFB.
   - Check the two scatterplots: the two stain vectors should sit on the two visible arms of the
     pixel cloud. If one vector is chasing noise, shrink the rectangle and retry.
   - Manually name them: **Stain 1 = `LFB`**, **Stain 2 = `Cresyl violet`**.
   - Set **Background** using the blank-glass patch (there's an auto-estimate; sanity-check it —
     it should be ~230–250 in all three channels, and roughly neutral).
5. Verify visually: `View → Brightness/Contrast`, then switch the viewer to the **LFB** channel
   (`View → Show channel` or the channel dropdown). White matter should be dark/bright and
   uniform; neuronal somata should nearly vanish. If cell bodies bleed into the LFB channel,
   the vectors are wrong — redo step 4 with a better region.
6. **Copy the resulting stain string.** `Automate → Show script editor`, then run:
   ```groovy
   print getCurrentImageData().getColorDeconvolutionStains().toString()
   ```
   Paste that string into the batch script (`LFB_batch_measure.groovy`, `STAIN_JSON` variable) so
   every image gets identical vectors.

---

## 4. Draw regional annotations

Create the classes first so the export is tidy: in the **Annotations** tab, `+` to add classes —
e.g. `CorpusCallosum`, `Cortex_M1`, `Striatum_WM_bundles`, `AnteriorCommissure`, `Fimbria`,
`CerebellarWM`, `OpticTract`.

Then per section:

1. **Brush/Polygon tool** (`B` / `P`) to outline each region. Wand tool (`W`) is tempting for
   corpus callosum but it follows *intensity* — which is the thing you're measuring — so it will
   bias KO regions to be drawn smaller. **Use geometry-based tools, not intensity-based ones.**
2. Assign the class (select annotation → click class in the list → `Set class`).
3. For structures where area itself is the readout (corpus callosum thickness/area), also record
   `Measure → Add shape features` so you get area and perimeter alongside intensity.
4. Consistency tricks:
   - Draw the same fixed-size ROI everywhere it's defensible. `Objects → Annotations → Specify annotation`
     lets you place an exact-dimension rectangle/ellipse (e.g. a 500 × 200 µm box centered on the CC
     at the midline) — far more reproducible than freehand for a paired WT/KO comparison.
   - Copy annotations between serial sections: select → `Ctrl+C`, open next image → `Ctrl+V`, then
     nudge. Saves enormous time and keeps ROI size constant.
   - Draw on the **cresyl violet channel or the original RGB**, never while viewing the LFB channel —
     otherwise you're picking your ROI based on the outcome variable.

> If ROI drawing becomes the bottleneck (it usually does at n > 20 slides), the alternative is
> **ABBA** (Aligning Big Brains & Atlases) — a free Fiji/QuPath tool that registers mouse sections
> to the Allen CCFv3 and imports every atlas region as a QuPath annotation automatically. Much
> more setup, much more rigorous, and it removes ROI-drawing bias entirely. Worth it if this
> figure is going in the manuscript. Say the word and I'll write that workflow up.

---

## 5. Build the myelin threshold classifier (once, applied to all)

`Classify → Pixel classification → Create thresholder`

| Setting | Suggested value | Notes |
|---|---|---|
| Resolution | ~1–2 µm/px | Full res is unnecessary and slow; keep it constant across all images |
| Channel | **LFB** (the deconvolved channel, not Red/Green/Blue) | The whole point of §3 |
| Prefilter | Gaussian | Suppresses per-pixel noise |
| Smoothing sigma | 1.0–2.0 | Higher = smoother masks; don't go so high you blur out thin fiber bundles |
| Threshold | set on WT tissue, see below | |
| Above threshold | `Myelin` | |
| Below threshold | `Ignore*` | `Ignore*` keeps it out of object creation but still counted in the denominator |

**Setting the threshold honestly:** open 3–4 images spanning your full range of staining
intensity (both genotypes), toggle the classifier overlay, and pick the value where the mask
matches what you'd call myelin *by eye in WT*. Then **freeze it** and never touch it again.
Report the exact value in your methods.

Save as e.g. `LFB_myelin` (`Save` button in the thresholder dialog). It's stored in the project's
`classifiers/pixel_classifiers/` folder.

**Sanity check before committing:** run it on a KO and a WT and eyeball the overlay. If the
threshold lands on a steep part of the intensity histogram, small staining differences will
produce large %area differences — that's the fragile regime. Check robustness by re-running at
threshold ±10% and confirming your WT/KO effect direction and significance don't flip (§8).

---

## 6. Measure and export

**GUI route:** select all annotations (`Ctrl+Alt+A`), then in the pixel classifier dialog press
**Measure**, or `Classify → Pixel classification → …` and add measurements. This writes
`LFB: Myelin area µm^2` and `LFB: Myelin %` onto each annotation.

**Batch route (recommended):** use `LFB_batch_measure.groovy` in this folder.
`Automate → Script editor` → open it → set the three CONFIG variables → `Run → Run for project`
→ select all images.

Then `Measure → Export measurements` → object type **Annotation** → select all images → export
to `LFB_measurements.tsv`.

### Which numbers to report

Report **both**, they fail in different ways:

- **% positive area** (`Myelin %`) — robust, intuitive, but **saturates**: in a densely myelinated
  corpus callosum both WT and KO can be ~100% positive and a real thinning is invisible.
- **Mean LFB optical density** within the ROI — continuous, sensitive to density changes in
  saturated regions, but more sensitive to staining batch. Get this from
  `Analyze → Calculate features → Add intensity features` with **Stain 1 (LFB)** checked. Run it
  once in the GUI, then `Automate → Create script from workflow` to grab the exact parameter JSON
  for the batch script.
- Plus, for corpus callosum specifically: **ROI area / thickness** as an independent structural
  readout.

---

## 7. Statistics

- **n = number of animals, not sections and not ROIs.** Average all sections for a given region
  within an animal → one value per animal per region → test on those. Treating sections as n is
  pseudoreplication and it's the most common reason this kind of figure gets flagged in review.
- Alternative if you want to keep section-level information: **linear mixed model** with
  `animal` as a random effect (`lme4::lmer(pct_myelin ~ genotype + (1|animal)`), and add
  `(1|slide)` if sections were spread across slides.
- Multiple regions → correct across regions (Holm or FDR), or pre-register corpus callosum as
  the primary region and treat the rest as exploratory.
- Report the staining batch and scanner session; if more than one batch was unavoidable, include
  batch as a covariate and say so.

---

## 8. Validation before you trust the numbers

1. **Threshold sensitivity:** re-run at threshold ±10% and ±20%. Effect size should shrink or
   grow smoothly, not reverse.
2. **Stain-vector sensitivity:** re-derive vectors from a *different* reference image, re-run,
   confirm the conclusion holds.
3. **Positive control region:** a structure you don't expect to change (e.g. optic tract, or a
   region unrelated to your OPC/oligodendrocyte hypothesis) should show no genotype difference.
   If everything moves in the same direction, you're measuring staining, not biology.
4. **Blinded re-scoring:** have a second person redraw ROIs on a subset; report the correlation.
5. **Orthogonal modality — do this.** Reviewers rarely accept LFB alone for a myelin claim.
   Back it with MBP or PLP immunostaining (also quantifiable in QuPath by the same route with
   DAB vectors), and ideally FluoroMyelin or EM g-ratio on a subset. Given the BCAT1 project's
   OPC→oligodendrocyte-maturation angle, pairing LFB area with an OLIG2 / CC1 / PDGFRA cell-count
   in the same regions would make the histology self-consistent with the RNA-seq maturation
   phenotype rather than just parallel to it.

---

## 9. Known failure modes

| Symptom | Cause | Fix |
|---|---|---|
| Cell bodies show up in the LFB channel | Stain vectors overlap / poorly fitted | Re-estimate on a region with cleaner separation of the two dyes; do not use the HE-exclusion option |
| Wildly different % between slides of the same genotype | Auto stain-vector or auto-background per image; or staining batch | Hard-code one stain string across all images (§3.6) |
| Everything reads ~100% positive | Threshold too permissive, or true saturation in dense WM | Lower threshold *and* switch primary metric to mean OD |
| Background glass counted as tissue | Background RGB set from tissue, not glass | Re-set background from an empty-glass patch |
| Folds, bubbles, knife marks inflating signal | Artifacts inside ROI | Draw a `Region*`-classed annotation over the artifact, or exclude that section — decide the rule *before* looking at genotype |
| Analysis is agonizingly slow | Classifier resolution set to full res | 1–2 µm/px is plenty for myelin area |

---

## Appendix: the Fiji / ImageJ route

Fiji can do this, and the math is identical — QuPath's deconvolution and Fiji's are both
implementations of Ruifrok & Johnston (2001). Two ways in:

**A. Gabriel Landini's Colour Deconvolution plugin.** `Image → Colour → Colour Deconvolution`,
choose `vectors = [User values]`, enter the same normalized OD vectors you got from QuPath.
Landini's **Colour Deconvolution 2** (available via its own Fiji update site) is the one to use —
it can output **32-bit absorbance** directly, and for a 2-dye stain you right-click at the third
vector prompt to set it to 0,0,0 and have it computed as the cross product.

**B. `LFB_CV_deconvolution.ijm` in this folder** — does the deconvolution arithmetic explicitly
in 32-bit OD with no plugin dependency, loops over a folder, reads ROI sets, and writes a CSV
of mean/median/SD LFB optical density plus % area above a fixed OD threshold. Paste your QuPath
stain vectors into the CONFIG block at the top.

### Three ways the Fiji route goes wrong

1. **The default 8-bit output is transmittance, not absorbance.** 255 = *no* stain. Mean Gray
   Value on that image is inversely and non-linearly related to dye amount, so group differences
   come out compressed and sometimes reversed. Either use Colour Deconvolution 2's 32-bit
   absorbance output, or convert: `32-bit → Min...(1) → Divide by 255 → Log → Multiply by −0.4343`.
2. **Skipping the white/background normalization.** ImageJ's plugin assumes I₀ = 255. If your
   scanner's blank glass sits at 238, every OD is offset — tolerable within one slide, fatal
   across a WT/KO comparison. Measure blank glass once and hard-code I₀.
3. **Splitting RGB and using the blue channel.** Very common, completely wrong here: cresyl
   violet is purple and contributes heavily to blue. That shortcut measures LFB + Nissl. This is
   the entire reason deconvolution is necessary.

### Practical limits vs. QuPath

Fiji has no real whole-slide support — Bio-Formats will open `.svs`/`.ndpi`, but full-resolution
sections will exhaust RAM. The workable pattern is **QuPath for ROI management and region export,
Fiji for the pixel math**: draw your regional annotations in QuPath, export each as a TIFF, point
the macro at that folder.

### Worth doing once

Run the same regions through both QuPath and the macro with identical stain vectors and I₀, and
correlate. They should agree closely. If they don't, one of the two has a background or
transmittance/absorbance problem, and you want to find that out now rather than in review. It also
makes a clean supplementary panel.

---

## Sources

- [Separating stains — QuPath documentation](https://qupath.readthedocs.io/en/stable/docs/tutorials/separating_stains.html)
- [Colour Deconvolution — ImageJ](https://imagej.net/plugins/colour-deconvolution)
- [Colour Deconvolution 2 (Landini) — GitHub](https://github.com/landinig/IJ-Colour_Deconvolution2)
- [Optimizing selection of unitary OD vectors — ImageJ docs](https://imagejdocu.list.lu/plugin/color/colour_deconvolution/optimizing_selection_of_unitary_optical_density_vectors/start)
- [QuPath Preprocessing wiki](https://github.com/qupath/qupath/wiki/Preprocessing)
- [Luxol Fast Blue quantification with Color Deconvolution — image.sc forum](https://forum.image.sc/t/luxol-fast-blue-quantification-with-color-deconvolution/33991)
- [EstimateStainVectors API (QuPath javadoc)](https://qupath.github.io/javadoc/docs/qupath/lib/analysis/algorithms/EstimateStainVectors.html)
- [Luxol fast blue stain — background](https://en.wikipedia.org/wiki/Luxol_fast_blue_stain)
