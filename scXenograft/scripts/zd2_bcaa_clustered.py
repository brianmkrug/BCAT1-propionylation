"""v2 of the BCAA heatmap: rows clustered by similarity (all genes mixed,
no pathway groups, no asterisks). Same pooled-z values as v1."""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.cluster.hierarchy import linkage, leaves_list

OUT = REPO_ROOT + "/scXenograft"; FIG = f"{OUT}/figures"
NB = 20
v = pd.read_csv(f"{OUT}/bcaa_gene_heatmap_values.csv", index_col=0)
kc = [f"K27M_bin{b}" for b in range(NB)]; oc = [f"KO_bin{b}" for b in range(NB)]
K = v[kc].values; O = v[oc].values
lfc = v["pb_log2FC"].values
genes = [g.replace("MUT", "MMUT/MUT") if g == "MUT" else g for g in v.index]

# cluster on the concatenated K27M+KO profile (correlation distance = group by shape)
feat = np.nan_to_num(np.hstack([K, O]))
Z = linkage(feat, method="average", metric="correlation")
order = leaves_list(Z)
K, O, lfc = K[order], O[order], lfc[order]
genes = [genes[i] for i in order]

VZ = 2.5
fig = plt.figure(figsize=(14.5, max(11, 0.26*len(genes))))
gs = GridSpec(1, 4, width_ratios=[1, 1, 0.06, 0.045], wspace=0.05)
axK, axO = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
for ax, Mx, ttl in [(axK, K, "K27M"), (axO, O, "K27M-KO (r1+r3)")]:
    im = ax.imshow(Mx, aspect="auto", cmap="RdBu_r", vmin=-VZ, vmax=VZ)
    ax.set_xticks([0, NB//2, NB-1]); ax.set_xticklabels(["early", "mid", "late"])
    ax.set_xlabel("pooled pseudotime bin"); ax.set_title(ttl, fontsize=11)
axK.set_yticks(range(len(genes))); axK.set_yticklabels(genes, fontsize=7); axO.set_yticks([])
axL = fig.add_subplot(gs[0, 2]); axL.imshow(lfc[:, None], aspect="auto", cmap="PuOr_r", vmin=-1.5, vmax=1.5)
axL.set_xticks([]); axL.set_yticks([]); axL.set_title("log2FC\nK27M/KO", fontsize=6.5)
fig.colorbar(im, cax=fig.add_subplot(gs[0, 3]), label="pooled z (per gene)")
fig.suptitle("v2 — BCAA metabolism genes clustered by similarity, K27M vs KO across pooled pseudotime\n"
             "pooled z per gene (shared scale); rows = hierarchical clustering (correlation distance); "
             "right strip = pseudobulk log2FC", fontsize=11)
fig.subplots_adjust(left=0.14, top=0.95)
fig.savefig(f"{FIG}/bcaa_gene_heatmap_clustered.png", dpi=170, bbox_inches="tight")
print("clustered order (top->bottom):")
print(", ".join(genes))
