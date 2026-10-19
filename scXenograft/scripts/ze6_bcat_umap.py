"""Pooled5 UMAP split K27M | KO, BCAT1 & BCAT2:
(A) log-normalised expression (count metric), common per-gene scale across genotypes
(B) pooled z-score (both genotypes together), common per-gene scale
"""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import numpy as np, pandas as pd, scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

W = REPO_ROOT + "/scXenograft/work"
FIG = REPO_ROOT + "/scXenograft/figures"
GENES = ["BCAT1", "BCAT2"]
emb = pd.read_csv(f"{W}/pooled5_dpt.csv.gz", index_col=0)
A = sc.read_h5ad(f"{W}/hsj019_norm.h5ad")[emb.index]
for g in GENES:
    emb[g] = np.asarray(A[:, g].X.todense()).ravel()
    emb[g+"_z"] = (emb[g] - emb[g].mean()) / emb[g].std()   # pooled z (both genotypes)
gt = [("K27M", "K27M"), ("K27M-KO", "KO (r1+r2+r3)")]
X, Y = emb.UMAP1.values, emb.UMAP2.values

def panel(mode, fname, title):
    fig = plt.figure(figsize=(11, 9))
    gs = GridSpec(2, 3, width_ratios=[1, 1, 0.05], wspace=0.08, hspace=0.14)
    for r, g in enumerate(GENES):
        col = g if mode == "expr" else g+"_z"
        v = emb[col].values
        if mode == "expr":
            vmin, vmax, cmap = 0, np.percentile(v, 99.5), "magma"
        else:
            lim = np.percentile(np.abs(v), 99); vmin, vmax, cmap = -lim, lim, "RdBu_r"
        im = None
        for c, (gg, ttl) in enumerate(gt):
            ax = fig.add_subplot(gs[r, c])
            m = (emb.genotype == gg).values
            o = np.argsort(v[m]) if mode == "expr" else np.argsort(np.abs(v[m]))
            im = ax.scatter(X[m][o], Y[m][o], c=v[m][o], s=5, cmap=cmap,
                            vmin=vmin, vmax=vmax, linewidths=0, rasterized=True)
            ax.set_xticks([]); ax.set_yticks([])
            for sp in ax.spines.values(): sp.set_visible(False)
            if r == 0: ax.set_title(ttl, fontsize=14)
            if c == 0: ax.set_ylabel(g, fontsize=15, fontweight="bold")
        cax = fig.add_subplot(gs[r, 2])
        cb = fig.colorbar(im, cax=cax)
        cb.set_label("log-norm expression" if mode == "expr" else "z-score (pooled)", fontsize=11)
        cb.ax.tick_params(labelsize=10)
    fig.suptitle(title, fontsize=14, y=0.98)
    fig.savefig(f"{FIG}/{fname}", dpi=170, bbox_inches="tight"); plt.close(fig)

panel("expr", "bcat_umap_expression.png",
      "BCAT1 / BCAT2 log-normalised expression on pooled UMAP, split by genotype\n"
      "(common colour scale per gene across K27M and KO)")
panel("z", "bcat_umap_zscore.png",
      "BCAT1 / BCAT2 z-score on pooled UMAP, split by genotype\n"
      "(pooled z per gene — both genotypes on one common scale)")
print("written")
for g in GENES:
    print(f"{g}: K27M mean {emb[emb.genotype=='K27M'][g].mean():.2f}, "
          f"KO mean {emb[emb.genotype=='K27M-KO'][g].mean():.2f}, "
          f"%det K27M {100*(emb[emb.genotype=='K27M'][g]>0).mean():.0f}, "
          f"KO {100*(emb[emb.genotype=='K27M-KO'][g]>0).mean():.0f}")
