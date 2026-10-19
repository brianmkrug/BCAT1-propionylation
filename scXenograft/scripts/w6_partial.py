"""Partial correlations controlling for detection depth (n_genes), on ranks.
Reports which earlier claims survive depth adjustment."""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import numpy as np, pandas as pd
from scipy.stats import rankdata, spearmanr
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

W = REPO_ROOT + "/scXenograft/work"
FIG = REPO_ROOT + "/scXenograft/figures"

def partial(x, y, z):
    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    rx, ry, rz = (rankdata(v[ok]) for v in (x, y, z))
    def res(a, b):
        A = np.c_[np.ones_like(b), b]
        return a - A @ np.linalg.lstsq(A, a, rcond=None)[0]
    return np.corrcoef(res(rx, rz), res(ry, rz))[0, 1]

rows = []
for w, ttl in [("k27m", "K27M"), ("ko13", "KO r1+r3")]:
    d = pd.read_csv(f"{W}/hsj019_{w}_dpt.csv.gz", index_col=0)
    zz = pd.read_csv(f"{W}/hsj019_{w}_ssgsea_z.csv.gz", index_col=0).loc[d.index]
    zz["CHOL_composite"] = zz[[c for c in zz.columns if c.startswith("CHOL_")]].mean(1)
    zz["FA_composite"] = zz[[c for c in zz.columns if c.startswith("FA_")]].mean(1)
    ng = d.n_genes.values.astype(float)
    pt = d.dpt_pseudotime.values

    print(f"\n=== {ttl}: Spearman with DPT, raw vs partialled on n_genes ===")
    for c in ["RGC", "OPC", "OL", "ASTRO", "PROLIF_sig", "BCAT1", "PDGFRA", "MBP"]:
        if c not in d.columns: continue
        v = d[c].values
        raw = spearmanr(pt, v).statistic
        par = partial(pt, v, ng)
        rows.append(dict(subset=w, feature=c, kind="lineage/marker",
                         rho_raw=raw, rho_partial=par))
        print(f"  {c:<14} raw {raw:+.3f}  partial {par:+.3f}  "
              f"({'survives' if abs(par) > 0.15 else 'collapses'})")
    for c in ["CHOL_composite", "FA_composite",
              "CHOL_Hallmark_cholesterol_homeostasis", "FA_Hallmark_fatty_acid_metabolism",
              "CTX_KEGGc_BCAA_degradation"]:
        v = zz[c].values
        raw = spearmanr(pt[np.isfinite(v)], v[np.isfinite(v)]).statistic
        par = partial(pt, v, ng)
        rows.append(dict(subset=w, feature=c, kind="ssGSEA",
                         rho_raw=raw, rho_partial=par))
        print(f"  {c[:36]:<36} raw {raw:+.3f}  partial {par:+.3f}  "
              f"({'survives' if abs(par) > 0.15 else 'collapses'})")

pd.DataFrame(rows).round(4).to_csv(f"{W}/depth_partial_correlations.csv", index=False)

# ---------- figure: the confound, shown ----------
fig, axes = plt.subplots(2, 3, figsize=(15, 8.5))
for row, (w, ttl) in enumerate([("k27m", "HSJ019 K27M"), ("ko13", "HSJ019 K27M-KO (r1+r3)")]):
    d = pd.read_csv(f"{W}/hsj019_{w}_dpt.csv.gz", index_col=0)
    zz = pd.read_csv(f"{W}/hsj019_{w}_ssgsea_z.csv.gz", index_col=0).loc[d.index]
    rn = pd.read_csv(f"{W}/hsj019_k27m_ssgsea_random_null.csv.gz", index_col=0) \
         if w == "k27m" else None
    comp = zz[[c for c in zz.columns if c.startswith("CHOL_")]].mean(1)

    ax = axes[row, 0]
    ax.scatter(d.n_genes, comp, s=2, alpha=.25, c="#C1272D", linewidths=0, rasterized=True)
    ax.set_xlabel("genes detected"); ax.set_ylabel("CHOL composite (z ssGSEA)")
    ax.set_title(f"{ttl}\nssGSEA vs detection depth "
                 f"(rho {spearmanr(d.n_genes, comp).statistic:+.2f})", fontsize=9)

    ax = axes[row, 1]
    o = np.argsort(d.n_genes.values)
    sm = ax.scatter(d.UMAP1_harmony.values[o], d.UMAP2_harmony.values[o],
                    c=d.n_genes.values[o], s=3, cmap="magma", linewidths=0, rasterized=True)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_visible(False)
    ax.set_title("genes detected on the UMAP", fontsize=9)
    fig.colorbar(sm, ax=ax, fraction=.045)

    ax = axes[row, 2]
    if rn is not None:
        for c, col in zip(rn.columns, ["#444444", "#888888", "#bbbbbb"]):
            v = (rn[c] - rn[c].mean()) / rn[c].std()
            ax.scatter(d.n_genes, v, s=2, alpha=.2, c=col, linewidths=0,
                       label=f"{c} (rho {spearmanr(d.n_genes, rn[c]).statistic:+.2f})",
                       rasterized=True)
        ax.legend(fontsize=6, frameon=False, markerscale=4)
        ax.set_xlabel("genes detected"); ax.set_ylabel("random-set ssGSEA (z)")
        ax.set_title("RANDOM gene sets show the same gradient", fontsize=9)
    else:
        ax.axis("off")
    for a in axes[row, [0, 2]]:
        for sp in ["top", "right"]: a.spines[sp].set_visible(False)
fig.suptitle("ssGSEA on sparse single-cell data tracks detection depth, not pathway activity",
             fontsize=12)
fig.tight_layout(); fig.savefig(f"{FIG}/hsj019_ssgsea_depth_confound.png", dpi=170)
plt.close(fig)
print("\nfigure written")
