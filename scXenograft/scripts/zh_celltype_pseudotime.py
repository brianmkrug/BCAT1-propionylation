"""Cell-type association along the pooled pseudotime axis, K27M vs KO side by side,
matched (fixed-width) pseudotime on x. Three representations:
  (A) hard-label stacked-area composition (argmax of lineage meta-scores)
  (B) soft-assignment stacked area (softmax of meta-scores -> expected composition)
  (C) unstacked smoothed line per cell type
Each with a cell-density track so sparse-support regions are visible.
"""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

W = REPO_ROOT + "/scXenograft/work"
FIG = REPO_ROOT + "/scXenograft/figures"
CLS = ["RGC", "GLIP", "OPC", "OL", "ASTRO", "EPEN", "NEURON"]
COL = dict(zip(CLS, ["#7f7f7f", "#9467bd", "#e6a010", "#f4d03f", "#2ca089",
                     "#e8746a", "#7fbf3f"]))
NBIN, SMOOTH = 24, 1.0

emb = pd.read_csv(f"{W}/pooled_dpt.csv.gz", index_col=0)
emb["label"] = emb[CLS].idxmax(1)
# softmax soft assignment per cell (temperature 1 on the z-scored meta-scores)
Z = emb[CLS].values
P = np.exp(Z - Z.max(1, keepdims=True)); P = P / P.sum(1, keepdims=True)
for i, c in enumerate(CLS):
    emb["soft_"+c] = P[:, i]

lo, hi = emb.dpt.min(), emb.dpt.max()
edges = np.linspace(lo, hi, NBIN+1); ctr = (edges[:-1]+edges[1:])/2
emb["fb"] = np.clip(np.digitize(emb.dpt, edges)-1, 0, NBIN-1)

def smooth(y):
    if SMOOTH <= 0: return y
    k = np.array([0.25, 0.5, 0.25])
    yy = y.copy()
    for _ in range(1):
        yy = np.convolve(np.r_[yy[0], yy, yy[-1]], k, "same")[1:-1]
    return yy

def comp(sub, cols):
    frac = np.zeros((len(cols), NBIN)); n = np.zeros(NBIN)
    for b in range(NBIN):
        d = sub[sub.fb == b]; n[b] = len(d)
        if len(d):
            for i, c in enumerate(cols):
                frac[i, b] = d[c].mean() if c.startswith("soft_") else (d.label == c).mean()
    return frac, n

GT = [("K27M", "K27M"), ("K27M-KO", "KO (r1+r3)")]
fig, axes = plt.subplots(3, 2, figsize=(15, 12), sharex=True)
for col, (g, ttl) in enumerate(GT):
    sub = emb[emb.genotype == g]
    # (A) hard stacked area
    frac, n = comp(sub, CLS)
    ax = axes[0, col]
    bottom = np.zeros(NBIN)
    for i, c in enumerate(CLS):
        y = smooth(frac[i])
        ax.fill_between(ctr, bottom, bottom+y, color=COL[c], label=c, lw=0)
        bottom += y
    ax.set_ylim(0, 1); ax.set_title(f"{ttl}\n(A) hard-label composition", fontsize=10)
    axc = ax.twinx(); axc.plot(ctr, n, color="k", lw=1, ls=":"); axc.set_ylim(0, n.max()*1.1)
    axc.set_ylabel("cells/bin", fontsize=7); axc.tick_params(labelsize=6)
    # (B) soft stacked area
    fracs, _ = comp(sub, ["soft_"+c for c in CLS])
    ax = axes[1, col]; bottom = np.zeros(NBIN)
    for i, c in enumerate(CLS):
        y = smooth(fracs[i])
        ax.fill_between(ctr, bottom, bottom+y, color=COL[c], lw=0)
        bottom += y
    ax.set_ylim(0, 1); ax.set_title("(B) soft-assignment composition (softmax)", fontsize=10)
    # (C) unstacked lines (hard fraction)
    ax = axes[2, col]
    for i, c in enumerate(CLS):
        ax.plot(ctr, smooth(frac[i]), color=COL[c], lw=2, label=c)
    ax.set_ylim(0, None); ax.set_title("(C) per-type fraction (lines)", fontsize=10)
    ax.set_xlabel("pooled pseudotime (matched)")
    axes[1, col].set_xlabel("")
for r in range(3):
    axes[r, 0].set_ylabel("fraction of cells")
axes[0, 1].legend(loc="center left", bbox_to_anchor=(1.12, 0.5), frameon=False, fontsize=8)
for ax in axes.ravel():
    for sp in ["top", "right"]: ax.spines[sp].set_visible(False)
fig.suptitle("Cell-type association along pooled pseudotime — K27M vs KO (matched x)\n"
             "dotted = cells per bin (sparse tails: K27M late, KO early)", fontsize=12)
fig.tight_layout(rect=[0, 0, 0.9, 0.96])
fig.savefig(f"{FIG}/celltype_along_pseudotime.png", dpi=170, bbox_inches="tight")
plt.close(fig)

# export per-bin composition tables
for g, ttl in GT:
    sub = emb[emb.genotype == g]
    frac, n = comp(sub, CLS)
    pd.DataFrame(frac.T, index=ctr.round(3), columns=CLS).assign(n_cells=n.astype(int)).round(3)\
      .to_csv(f"{FIG.replace('/figures','')}/celltype_composition_{ttl.split()[0]}.csv")
print("figure written")
print("\nlabel counts by genotype:")
print(pd.crosstab(emb.label, emb.genotype).to_string())
