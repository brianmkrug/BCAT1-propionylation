"""bcaa clustered (5 lib), ECHS1 & HSD17B10 below LDHB, split into 3 groups with
white gaps (after ALDH3A2 and after HSD17B10), with highlighted gene labels."""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.cluster.hierarchy import linkage, leaves_list

OUT = REPO_ROOT + "/scXenograft"; FIG = f"{OUT}/figures"
NB = 20
GREEN = {"BCAT1", "ACAD8", "ECHS1", "HADHA", "HSD17B10", "ACAT1"}
PURPLE = {"PCCA", "PCCB", "MUT", "MCEE"}          # MUT displays as MMUT/MUT
GCOL, PCOL = "#1a7d1a", "#6a1b9a"
GAP = 1   # blank rows between groups

v = pd.read_csv(f"{OUT}/bcaa_gene_heatmap_values_5lib.csv", index_col=0)
kc = [f"K27M_bin{b}" for b in range(NB)]; oc = [f"KO_bin{b}" for b in range(NB)]
feat = np.nan_to_num(np.hstack([v[kc].values, v[oc].values]))
clustered = list(v.index[leaves_list(linkage(feat, method="average", metric="correlation"))])
order = [g for g in clustered if g not in ("ECHS1", "HSD17B10")]
i = order.index("LDHB"); order = order[:i+1] + ["ECHS1", "HSD17B10"] + order[i+1:]

# 3 groups: split after ALDH3A2 and after HSD17B10
b1 = order.index("ALDH3A2") + 1
b2 = order.index("HSD17B10") + 1
groups = [order[:b1], order[b1:b2], order[b2:]]

def stack(colset):
    rows, ypos, ylab, ycol = [], [], [], []
    r = 0
    for gi, grp in enumerate(groups):
        for g in grp:
            rows.append(v.loc[g, colset].values); ypos.append(r)
            ylab.append("MMUT/MUT" if g == "MUT" else g)
            ycol.append(GCOL if g in GREEN else PCOL if g in PURPLE else "k")
            r += 1
        if gi < len(groups)-1:
            for _ in range(GAP): rows.append(np.full(len(colset), np.nan)); r += 1
    return np.array(rows), ypos, ylab, ycol

K, ypos, ylab, ycol = stack(kc)
O, _, _, _ = stack(oc)
lfccol, _, _, _ = stack(["pb_log2FC"])

cmap = plt.cm.RdBu_r.copy(); cmap.set_bad("white")
cmapL = plt.cm.PuOr_r.copy(); cmapL.set_bad("white")
VZ = 2.5
fig = plt.figure(figsize=(14.5, max(11, 0.26*K.shape[0])))
gs = GridSpec(1, 4, width_ratios=[1, 1, 0.06, 0.045], wspace=0.05)
axK, axO = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
for ax, Mx, ttl in [(axK, K, "K27M"), (axO, O, "K27M-KO (r1+r2+r3)")]:
    im = ax.imshow(np.ma.masked_invalid(Mx), aspect="auto", cmap=cmap, vmin=-VZ, vmax=VZ)
    ax.set_xticks([0, NB//2, NB-1]); ax.set_xticklabels(["early", "mid", "late"])
    ax.set_xlabel("pooled pseudotime bin"); ax.set_title(ttl, fontsize=11)
axK.set_yticks(ypos); axK.set_yticklabels(ylab, fontsize=7); axO.set_yticks([])
for t, c in zip(axK.get_yticklabels(), ycol):
    t.set_color(c)
    if c != "k": t.set_fontweight("bold")
axL = fig.add_subplot(gs[0, 2])
axL.imshow(np.ma.masked_invalid(lfccol), aspect="auto", cmap=cmapL, vmin=-1.5, vmax=1.5)
axL.set_xticks([]); axL.set_yticks([]); axL.set_title("log2FC\nK27M/KO", fontsize=6.5)
fig.colorbar(im, cax=fig.add_subplot(gs[0, 3]), label="pooled z (per gene)")
fig.suptitle("BCAA metabolism genes (all 5 libraries), K27M vs KO across pooled pseudotime\n"
             "clustering order (ECHS1/HSD17B10 moved below LDHB); 3 groups; strip = pseudobulk log2FC",
             fontsize=11, y=1.01)
fig.subplots_adjust(left=0.14, top=0.95)
fig.savefig(f"{FIG}/bcaa_gene_heatmap_clustered_5lib_reordered.png", dpi=170, bbox_inches="tight")
print("groups:", [len(g) for g in groups], "| green", sorted(GREEN & set(order)),
      "| purple", sorted(PURPLE & set(order)))
