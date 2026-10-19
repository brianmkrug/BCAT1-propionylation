"""Pooled5 UMAP coloured by cell-type panel classification (argmax lineage label).
(A) scatter dots by label, (B) hexbin coloured by majority label per hex,
plus soft per-lineage score small-multiples (honest, continuous)."""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

W = REPO_ROOT + "/scXenograft/work"
FIG = REPO_ROOT + "/scXenograft/figures"
CLS = ["RGC","GLIP","OPC","OL","ASTRO","EPEN","NEURON"]
COL = ["#7f7f7f","#9467bd","#e6a010","#f4d03f","#2ca089","#e8746a","#7fbf3f"]
cmap = ListedColormap(COL)

e = pd.read_csv(f"{W}/pooled5_dpt.csv.gz", index_col=0)
e["code"] = e[CLS].values.argmax(1)
e["label"] = [CLS[c] for c in e.code]
x, y = e.UMAP1.values, e.UMAP2.values

fig, axes = plt.subplots(1, 2, figsize=(15, 6.5))
# (A) scatter dots by label
o = np.random.default_rng(0).permutation(len(e))   # shuffle so no type is drawn on top
axes[0].scatter(x[o], y[o], c=[COL[c] for c in e.code.values[o]], s=3.5,
                linewidths=0, rasterized=True)
axes[0].set_title("(A) cells coloured by panel classification (argmax label)", fontsize=11)
# (B) hexbin coloured by majority label
from scipy import stats
def majority(vals): return stats.mode(vals, keepdims=False).mode
hb = axes[1].hexbin(x, y, C=e.code.values, reduce_C_function=majority, gridsize=55,
                    cmap=cmap, vmin=-0.5, vmax=6.5, mincnt=1, linewidths=0.15, edgecolors="white")
axes[1].set_title("(B) hexbin coloured by majority label (de-overplotted)", fontsize=11)
for ax in axes:
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_visible(False)
h = [plt.Line2D([],[],marker="o",ls="",color=COL[i],label=CLS[i]) for i in range(len(CLS))]
fig.legend(handles=h, loc="lower center", ncol=7, frameon=False, fontsize=10, bbox_to_anchor=(0.5,-0.02))
fig.suptitle("HSJ019 pooled UMAP (all 5 libraries) — cell-type panel classification", fontsize=13)
fig.tight_layout(rect=[0,0.03,1,1])
fig.savefig(f"{FIG}/pooled5_umap_celltype.png", dpi=170, bbox_inches="tight"); plt.close(fig)

# soft per-lineage scores (continuous, honest)
fig, axes = plt.subplots(2, 4, figsize=(17, 8))
for ax, c in zip(axes.ravel(), CLS):
    v = e[c].values; lo, hi = np.percentile(v, [2, 98]); od = np.argsort(v)
    sm = ax.scatter(x[od], y[od], c=np.clip(v[od], lo, hi), s=3, cmap="magma",
                    linewidths=0, rasterized=True)
    ax.set_title(f"{c} score", fontsize=10); ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_visible(False)
    fig.colorbar(sm, ax=ax, fraction=.045)
axes.ravel()[-1].axis("off")
fig.suptitle("HSJ019 pooled UMAP — lineage meta-scores (continuous; the graded reality behind the labels)",
             fontsize=12)
fig.tight_layout(); fig.savefig(f"{FIG}/pooled5_umap_lineage_scores.png", dpi=160, bbox_inches="tight")
print("written")
print(e.label.value_counts().to_string())
