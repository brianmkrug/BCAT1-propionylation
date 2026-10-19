"""Raw ssGSEA vs z-scored colour scales on the Harmony UMAPs, with a
quantitative comparison of how much of the colour range each actually uses.
"""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

W = REPO_ROOT + "/scXenograft/work"
FIG = REPO_ROOT + "/scXenograft/figures"
SETS = [("k27m", "HSJ019 K27M"), ("ko13", "HSJ019 K27M-KO (r1+r3)")]
KEY = ["CHOL_Hallmark_cholesterol_homeostasis", "CHOL_Reactome_cholesterol_biosynthesis",
       "CHOL_WP_mevalonate_pathway", "FA_Hallmark_fatty_acid_metabolism",
       "FA_Reactome_mito_beta_oxidation", "FA_KEGGc_fatty_acid_biosynthesis"]

D, RAW, Z = {}, {}, {}
for w, _ in SETS:
    D[w] = pd.read_csv(f"{W}/hsj019_{w}_annotated.csv.gz", index_col=0)
    RAW[w] = pd.read_csv(f"{W}/hsj019_{w}_ssgsea_raw.csv.gz", index_col=0).loc[D[w].index]
    Z[w] = pd.read_csv(f"{W}/hsj019_{w}_ssgsea_z.csv.gz", index_col=0).loc[D[w].index]

# ---------- quantitative comparison ----------
print("How much of the colour range does each scaling actually use?")
print("(fraction of the full min-max span covered by the central 90% of cells;")
print(" higher = more of the colour bar carries information)\n")
rows = []
for w, ttl in SETS:
    for k in RAW[w].columns:
        r, z = RAW[w][k], Z[w][k]
        span_r = (np.nanpercentile(r, 95) - np.nanpercentile(r, 5)) / (np.nanmax(r) - np.nanmin(r))
        cv = r.std() / r.mean()          # raw scores sit on a large offset
        rows.append(dict(subset=w, signature=k, raw_mean=r.mean(), raw_sd=r.std(),
                         raw_CV=cv, frac_range_used=span_r,
                         z_p5=np.nanpercentile(z, 5), z_p95=np.nanpercentile(z, 95)))
res = pd.DataFrame(rows)
print(res.groupby("subset")[["raw_mean", "raw_sd", "raw_CV", "frac_range_used"]]
      .mean().round(4).to_string())
res.round(4).to_csv(f"{W}/ssgsea_scale_comparison.csv", index=False)

# cross-subset comparability of raw scores
print("\nRaw score offsets differ between the two integrations "
      "(same signature, different subset):")
for k in KEY[:3]:
    print(f"  {k[:44]:<44} K27M {RAW['k27m'][k].mean():8.1f}  "
          f"KO {RAW['ko13'][k].mean():8.1f}  diff {RAW['k27m'][k].mean()-RAW['ko13'][k].mean():+7.1f}")

# ---------- Fig: raw vs z, side by side ----------
for w, ttl in SETS:
    d, r, z = D[w], RAW[w], Z[w]
    fig, axes = plt.subplots(2, len(KEY), figsize=(4*len(KEY), 8.5))
    for j, k in enumerate(KEY):
        for i, (mat, lab, cmap) in enumerate([(r, "raw ssGSEA ES", "viridis"),
                                              (z, "z-scored", "RdBu_r")]):
            v = mat[k].values
            ax = axes[i, j]
            if i == 0:
                lo, hi = np.nanmin(v), np.nanmax(v)   # honest full range, no clipping
            else:
                lo, hi = -2.5, 2.5
            bad = ~np.isfinite(v)
            if bad.any():   # cells with no detected set gene -> grey, not imputed
                ax.scatter(d.UMAP1_harmony.values[bad], d.UMAP2_harmony.values[bad],
                           c="#dddddd", s=3, linewidths=0, rasterized=True)
            o = np.argsort(np.where(bad, -np.inf, v))[bad.sum():]
            sm = ax.scatter(d.UMAP1_harmony.values[o], d.UMAP2_harmony.values[o],
                            c=np.clip(v[o], lo, hi), s=3, cmap=cmap, vmin=lo, vmax=hi,
                            linewidths=0, rasterized=True)
            if i == 0 and bad.any():
                ax.text(.02, .02, f"{100*bad.mean():.1f}% n.d.", transform=ax.transAxes,
                        fontsize=6, color="#777777")
            ax.set_xticks([]); ax.set_yticks([])
            for sp in ax.spines.values(): sp.set_visible(False)
            if i == 0:
                ax.set_title(k.replace("_", " ")[:34], fontsize=8)
            if j == 0:
                ax.set_ylabel(lab, fontsize=9)
            fig.colorbar(sm, ax=ax, fraction=.045)
    fig.suptitle(f"{ttl} — raw ssGSEA (top, full range) vs z-scored (bottom, ±2.5 SD)",
                 fontsize=12)
    fig.tight_layout(); fig.savefig(f"{FIG}/hsj019_{w}_ssgsea_raw_vs_z.png", dpi=160)
    plt.close(fig)

# ---------- Fig: full panel, z-scored ----------
for w, ttl in SETS:
    d, z = D[w], Z[w]
    cols = list(z.columns)
    n = len(cols); ncol = 6; nrow = int(np.ceil(n/ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.1*ncol, 2.9*nrow))
    for ax, k in zip(axes.ravel(), cols):
        v = z[k].values
        bad = ~np.isfinite(v)
        if bad.any():
            ax.scatter(d.UMAP1_harmony.values[bad], d.UMAP2_harmony.values[bad],
                       c="#dddddd", s=2.5, linewidths=0, rasterized=True)
        o = np.argsort(np.where(bad, -np.inf, v))[bad.sum():]
        sm = ax.scatter(d.UMAP1_harmony.values[o], d.UMAP2_harmony.values[o],
                        c=np.clip(v[o], -2.5, 2.5), s=2.5, cmap="RdBu_r",
                        vmin=-2.5, vmax=2.5, linewidths=0, rasterized=True)
        ax.set_title(k.replace("_", " ")[:36], fontsize=7.5)
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values(): sp.set_visible(False)
    for ax in axes.ravel()[n:]: ax.axis("off")
    fig.colorbar(sm, ax=axes.ravel().tolist(), fraction=.012, label="z-scored ssGSEA")
    fig.suptitle(f"{ttl} — metabolic panel ssGSEA (z-scored within subset)", fontsize=13)
    fig.savefig(f"{FIG}/hsj019_{w}_ssgsea_panel_z.png", dpi=160, bbox_inches="tight")
    plt.close(fig)
print("\nfigures written")
