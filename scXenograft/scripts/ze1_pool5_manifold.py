"""v2 of pooled_manifold_representation, all 5 libraries (KO r1+r2+r3).
genotype UMAP | pooled5 pseudotime | per-genotype density | per-bin representation."""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

W = REPO_ROOT + "/scXenograft/work"
FIG = REPO_ROOT + "/scXenograft/figures"
GT = {"K27M": "#C1272D", "K27M-KO": "#2B6CB0"}
NB = 20
emb = pd.read_csv(f"{W}/pooled5_dpt.csv.gz", index_col=0)
emb["bin"] = pd.qcut(emb.dpt, NB, labels=False, duplicates="drop")

# per-bin representation (fixed pseudotime, per genotype + per KO library)
rep = emb.groupby("bin").apply(lambda d: pd.Series({
    "dpt_mid": d.dpt.mean(),
    "n_K27M": int((d.genotype == "K27M").sum()),
    "n_KO": int((d.genotype == "K27M-KO").sum())})).reset_index()

def scat(ax, x, y, c, **k):
    ax.scatter(x, y, c=c, s=3, linewidths=0, rasterized=True, **k)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_visible(False)

fig, axes = plt.subplots(1, 4, figsize=(19, 4.6))
for g, sub in emb.groupby("genotype"):
    scat(axes[0], sub.UMAP1, sub.UMAP2, GT[g], label=g)
axes[0].legend(markerscale=4, frameon=False, fontsize=8); axes[0].set_title("genotype", fontsize=10)
o = np.argsort(emb.dpt.values)
sm = axes[1].scatter(emb.UMAP1.values[o], emb.UMAP2.values[o], c=emb.dpt.values[o],
                     s=3, cmap="viridis", linewidths=0, rasterized=True)
axes[1].set_xticks([]); axes[1].set_yticks([])
for sp in axes[1].spines.values(): sp.set_visible(False)
axes[1].set_title("pooled DPT pseudotime", fontsize=10); fig.colorbar(sm, ax=axes[1], fraction=.045)
for g, sub in emb.groupby("genotype"):
    axes[2].hist(sub.dpt, bins=40, range=(0, 1), density=True, histtype="step", lw=2,
                 color=GT[g], label="KO (r1+r2+r3)" if g == "K27M-KO" else g)
axes[2].set_xlabel("pseudotime"); axes[2].set_ylabel("density"); axes[2].legend(frameon=False, fontsize=8)
axes[2].set_title("density per genotype\n(the compositional bias)", fontsize=10)
for sp in ["top", "right"]: axes[2].spines[sp].set_visible(False)
ax = axes[3]
ax.bar(rep.dpt_mid, rep.n_K27M, width=.02, color=GT["K27M"], label="K27M")
ax.bar(rep.dpt_mid, -rep.n_KO, width=.02, color=GT["K27M-KO"], label="KO (r1+r2+r3)")
ax.axhline(0, lw=.6, c="k"); ax.axhline(40, lw=.6, ls=":", c="k"); ax.axhline(-40, lw=.6, ls=":", c="k")
ax.set_xlabel("pseudotime"); ax.set_ylabel("cells per bin"); ax.legend(frameon=False, fontsize=8)
ax.set_title("per-bin representation\n(dotted = n=40 floor)", fontsize=10)
for sp in ["top", "right"]: ax.spines[sp].set_visible(False)
fig.suptitle("HSJ019 pooled K27M + KO manifold, common pseudotime axis — ALL 5 libraries (KO r1+r2+r3)", fontsize=12)
fig.tight_layout(); fig.savefig(f"{FIG}/pooled_manifold_representation_5lib.png", dpi=160, bbox_inches="tight")
print("written; total", len(emb), "cells")
print(emb.groupby(["genotype","sample"],observed=True).size().to_string())
print("\nmin cells/genotype across bins: K27M", int(rep.n_K27M.min()), " KO", int(rep.n_KO.min()))
