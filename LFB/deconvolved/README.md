# `deconvolved/` — contents

All 36 images (WT-1,2,3,7,8 · KO-2,3,5,7 × CORTEX/OB/CB/PONS), processed with identical
parameters: QuPath stain vectors, I₀ = 242,242,242, crop 695×780, tissue mask total OD > 0.08.

| Folder | Contents | Type |
|---|---|---|
| `RGB_cropped/` | `<name>_RGB.tif` | 8-bit RGB — **draw your ROIs here** |
| `LFB_OD/` | `<name>_LFB_OD.tif` | 32-bit float, myelin optical density |
| `CV_OD/` | `<name>_CV_OD.tif` | 32-bit float, cresyl violet optical density |
| `RESID_OD/` | `<name>_RESID_OD.tif` | 32-bit float, residual — QC only |
| `QC_panels/` | `<name>_qc.png` | RGB │ LFB │ CV │ residual, side by side |

Top level: `lfb_summary.csv` (per image), `lfb_by_group.csv`, `lfb_stats.csv`, `lfb_plot.png`.

---

## Read this before drawing ROIs

**Everything here is cropped to 695 × 780.** The originals in `image_data/` are 867 × 846 —
the crop removes the burned-in QuPath overlay (thumbnail inset, scale bar, coordinate box).
ROIs drawn on an original file will **not** align with these images. Always open from
`RGB_cropped/`.

**Draw on the RGB or the CV image, never on the LFB image.** Choosing an ROI boundary while
looking at the LFB channel means selecting on your outcome variable, which will manufacture a
genotype difference. Same reason not to use the wand tool — it follows intensity.

**Name ROIs identically across animals.** `measure_ROIs.ijm` groups on the ROI name, so
`corpus_callosum` in one animal and `CC` in another become different structures.

**Prefer fixed-size ROIs where the anatomy allows it.** Draw one, save it, load the same `.zip`
onto the next animal and reposition rather than redraw. Given that between-animal variance was
39% of the total in the whole-field analysis, ROI-size consistency is where you'll recover the
most precision.

**Blind yourself.** Rename files to codes before drawing, and unblind only at the stats step.
With n = 5 vs 4 and a known expected direction, unblinded ROI placement is the weakest link in
the whole pipeline.

---

## Measuring

1. Open a file from `RGB_cropped/`.
2. Polygon (`P`) or brush (`B`) to outline; press `T` to add to ROI Manager; rename each ROI.
3. ROI Manager → More ▸ → Save… → `<basename>.zip`, where basename drops the `_RGB`
   suffix (e.g. `WT-1-PONS.zip`). Collect all `.zip` files in one folder.
4. Run `../measure_ROIs.ijm` (Fiji → `[` script editor → Language: **IJ1 Macro** → Run).
   It applies each ROI to the matching LFB and CV images and writes a CSV with mean/median/SD
   LFB OD, % area above 0.25 OD, and the LFB:CV ratio.

## Sanity checks once you have ROI data

- **Residual channel:** open a few from `RESID_OD/`. Flat noise = vectors are good. Recognizable
  anatomy = re-estimate stain vectors before trusting anything.
- **Compare ROI results to `lfb_summary.csv`.** The whole-field numbers had KO below WT in all
  four fields. If matched ROIs abolish that difference, the whole-field effect was framing, not
  biology — which is exactly what this exercise is designed to find out.
- **The staining-batch question still stands.** ROIs fix field composition, not the fact that
  genotype is perfectly confounded with animal and slide. If WT and KO were stained in separate
  runs, matched ROIs will not rescue it.
