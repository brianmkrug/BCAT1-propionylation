"""BCAT1 / BCAT2 only, paired K27M | KO pooled-pseudotime heatmap (pooled-z per gene)."""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

OUT = REPO_ROOT + "/scXenograft"; FIG = f"{OUT}/figures"
NB = 20
v = pd.read_csv(f"{OUT}/bcaa_gene_heatmap_values.csv", index_col=0)
genes = ["BCAT1", "BCAT2"]
kc = [f"K27M_bin{b}" for b in range(NB)]; oc = [f"KO_bin{b}" for b in range(NB)]
K = v.loc[genes, kc].values; O = v.loc[genes, oc].values
lfc = v.loc[genes, "pb_log2FC"].values

VZ = 2.5
fig = plt.figure(figsize=(13, 2.2))
gs = GridSpec(1, 4, width_ratios=[1, 1, 0.06, 0.045], wspace=0.06)
axK, axO = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
for ax, Mx, ttl in [(axK, K, "K27M"), (axO, O, "K27M-KO (r1+r3)")]:
    im = ax.imshow(Mx, aspect="auto", cmap="RdBu_r", vmin=-VZ, vmax=VZ)
    ax.set_xticks([0, NB//2, NB-1]); ax.set_xticklabels(["early", "mid", "late"])
    ax.set_xlabel("pooled pseudotime bin"); ax.set_title(ttl, fontsize=11)
axK.set_yticks(range(len(genes))); axK.set_yticklabels(genes, fontsize=11); axO.set_yticks([])
axL = fig.add_subplot(gs[0, 2]); axL.imshow(lfc[:, None], aspect="auto", cmap="PuOr_r", vmin=-1.5, vmax=1.5)
axL.set_xticks([]); axL.set_yticks([]); axL.set_title("log2FC\nK27M/KO", fontsize=6.5)
fig.colorbar(im, cax=fig.add_subplot(gs[0, 3]), label="pooled z (per gene)")
fig.suptitle("BCAT1 / BCAT2 across pooled pseudotime, K27M vs KO (pooled z per gene; strip = pseudobulk log2FC)",
             fontsize=11, y=1.18)
fig.subplots_adjust(left=0.10, top=0.78)
fig.savefig(f"{FIG}/bcaa_gene_heatmap_bcat_only.png", dpi=170, bbox_inches="tight")
print("written; log2FC BCAT1 %.2f  BCAT2 %.2f" % (lfc[0], lfc[1]))
