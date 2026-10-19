"""BCAT1/BCAT2 (5 lib) — K27M stacked on top of K27M-KO, shared x-axis, thin rows."""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

OUT = REPO_ROOT + "/scXenograft"; FIG = f"{OUT}/figures"
NB = 20; GENES = ["BCAT1", "BCAT2"]; VZ = 2.5
v = pd.read_csv(f"{OUT}/bcaa_gene_heatmap_values_5lib.csv", index_col=0)
kc = [f"K27M_bin{b}" for b in range(NB)]; oc = [f"KO_bin{b}" for b in range(NB)]
mats = {"K27M": v.loc[GENES, kc].values, "K27M-KO (r1+r2+r3)": v.loc[GENES, oc].values}
lfc = v.loc[GENES, "pb_log2FC"].values[:, None]

fig = plt.figure(figsize=(9, 2.9))
gs = GridSpec(2, 2, width_ratios=[1, 0.03], height_ratios=[1, 1], hspace=0.28, wspace=0.02)
im = None
for r, (name, M) in enumerate(mats.items()):
    axHM = fig.add_subplot(gs[r, 0]); axL = fig.add_subplot(gs[r, 1])
    im = axHM.imshow(M, aspect="auto", cmap="RdBu_r", vmin=-VZ, vmax=VZ)
    axL.imshow(lfc, aspect="auto", cmap="PuOr_r", vmin=-1.5, vmax=1.5)
    for ax in (axHM, axL):
        for sp in ax.spines.values(): sp.set_visible(True); sp.set_edgecolor("k"); sp.set_linewidth(0.9)
    axHM.set_yticks(range(len(GENES))); axHM.set_yticklabels(GENES, fontsize=12)
    axL.set_xticks([]); axL.set_yticks([])
    axHM.set_title(name, fontsize=13, pad=4)
    bottom = (r == len(mats)-1)
    axHM.set_xticks([0, NB//2, NB-1])
    axHM.set_xticklabels(["early", "mid", "late"] if bottom else [], fontsize=12)
    if bottom: axHM.set_xlabel("pooled pseudotime bin", fontsize=13)
    if r == 0: axL.set_title("log2FC", fontsize=8, pad=4)
# colorbars, right margin
cax1 = fig.add_axes([0.92, 0.55, 0.018, 0.33]); cb1 = fig.colorbar(im, cax=cax1)
cb1.set_label("pooled z", fontsize=11); cb1.ax.tick_params(labelsize=9)
cax2 = fig.add_axes([0.92, 0.12, 0.018, 0.33])
cb2 = fig.colorbar(ScalarMappable(Normalize(-1.5, 1.5), "PuOr_r"), cax=cax2)
cb2.set_label("log2FC (K27M/KO)", fontsize=11); cb2.ax.tick_params(labelsize=9)
fig.suptitle("BCAT1 / BCAT2 across pooled pseudotime (all 5 libraries); pooled z per gene",
             fontsize=12, y=1.02)
fig.subplots_adjust(left=0.10, right=0.88, top=0.86, bottom=0.18)
fig.savefig(f"{FIG}/bcaa_gene_heatmap_bcat_only_5lib_stacked.png", dpi=170, bbox_inches="tight")
print("written; BCAT1 log2FC %.2f, BCAT2 %.2f" % (lfc[0, 0], lfc[1, 0]))
