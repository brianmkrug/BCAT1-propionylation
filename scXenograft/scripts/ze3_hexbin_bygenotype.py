"""Majority-label hexbin UMAP split by genotype (K27M | KO), shared hex grid."""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from scipy import stats

W = REPO_ROOT + "/scXenograft/work"
FIG = REPO_ROOT + "/scXenograft/figures"
CLS = ["RGC","GLIP","OPC","OL","ASTRO","EPEN","NEURON"]
COL = ["#7f7f7f","#9467bd","#e6a010","#f4d03f","#2ca089","#e8746a","#7fbf3f"]
cmap = ListedColormap(COL)

e = pd.read_csv(f"{W}/pooled5_dpt.csv.gz", index_col=0)
e["code"] = e[CLS].values.argmax(1)
# shared extent + gridsize so hexes align across panels
pad = 0.5
ext = [e.UMAP1.min()-pad, e.UMAP1.max()+pad, e.UMAP2.min()-pad, e.UMAP2.max()+pad]
GS = 48
def majority(v): return stats.mode(v, keepdims=False).mode

fig, axes = plt.subplots(1, 2, figsize=(16, 6.8))
for ax, (g, ttl) in zip(axes, [("K27M","K27M"), ("K27M-KO","KO (r1+r2+r3)")]):
    s = e[e.genotype == g]
    ax.hexbin(s.UMAP1, s.UMAP2, C=s.code.values, reduce_C_function=majority,
              gridsize=GS, extent=ext, cmap=cmap, vmin=-0.5, vmax=6.5,
              mincnt=3, linewidths=0.15, edgecolors="white")
    ax.set_xlim(ext[0], ext[1]); ax.set_ylim(ext[2], ext[3])
    ax.set_title(f"{ttl}  (n={len(s)})", fontsize=12)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_visible(False)
h = [plt.Line2D([],[],marker="h",ls="",color=COL[i],label=CLS[i],ms=10) for i in range(len(CLS))]
fig.legend(handles=h, loc="lower center", ncol=7, frameon=False, fontsize=11, bbox_to_anchor=(0.5,-0.01))
fig.suptitle("HSJ019 pooled UMAP — majority cell-type per hex, split by genotype "
             "(shared hex grid; hexes with <3 cells blank)", fontsize=13)
fig.tight_layout(rect=[0,0.04,1,1])
fig.savefig(f"{FIG}/pooled5_hexbin_celltype_bygenotype.png", dpi=170, bbox_inches="tight")
print("written")
print((pd.crosstab([CLS[c] for c in e.code], e.genotype, normalize="columns")*100).round(1).to_string())
