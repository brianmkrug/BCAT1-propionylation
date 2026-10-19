"""Trimmed grouped paired heatmap (all 5 libs): cholesterol biosynthesis,
citrate cycle (TCA), sphingolipid metabolism only."""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
W = REPO_ROOT + "/scXenograft/work"
FIG = REPO_ROOT + "/scXenograft/figures"
NB = 20
d = pd.read_csv(f"{W}/pooled5_dpt.csv.gz", index_col=0)
d["bin"] = pd.qcut(d.dpt, NB, labels=False, duplicates="drop")
D = pd.read_csv(f"{W}/pooled5_metab_resid.csv.gz", index_col=0).loc[d.index]

ROWS = [("CHOL_Reactome_cholesterol_biosynthesis", "Cholesterol biosynthesis (Reactome)"),
        ("KEGG_Citrate cycle (TCA cycle)", "Citrate cycle (TCA)"),
        ("KEGG_Sphingolipid metabolism", "Sphingolipid metabolism"),
        ("KEGG_Oxidative phosphorylation", "Oxidative phosphorylation")]
paths = [p for p, _ in ROWS]; labels = [l for _, l in ROWS]
row_mean = np.array([np.nanmean(D[p].values) for p in paths])[:, None]
def mat(gt):
    return np.vstack([[np.nanmean(D[p].values[(d.genotype==gt).values & (d.bin==b).values])
                       for b in range(NB)] for p in paths])
K = mat("K27M") - row_mean; O = mat("K27M-KO") - row_mean
vmax = float(np.nanpercentile(np.abs(np.r_[K.ravel(), O.ravel()]), 98))
nK = [int(((d.genotype=="K27M")&(d.bin==b)).sum()) for b in range(NB)]
nO = [int(((d.genotype=="K27M-KO")&(d.bin==b)).sum()) for b in range(NB)]

fig = plt.figure(figsize=(15.5, 3.4))
gs = GridSpec(1, 3, width_ratios=[1, 1, 0.045], wspace=0.06)
for j, (Mx, ttl) in enumerate([(K, "K27M"), (O, "K27M-KO (r1+r2+r3)")]):
    ax = fig.add_subplot(gs[0, j])
    im = ax.imshow(Mx, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    ax.set_xticks([0, NB//2, NB-1]); ax.set_xticklabels(["early", "mid", "late"])
    ax.set_xlabel("pooled pseudotime bin")
    if j == 0: ax.set_yticks(range(len(paths))); ax.set_yticklabels(labels, fontsize=9.5)
    else: ax.set_yticks([])
    ax.set_title(ttl, fontsize=11)
fig.colorbar(im, cax=fig.add_subplot(gs[0, 2]), label="control-bin score, row-centred\n(shared scale; ±98th pct)")
fig.suptitle("Metabolic pathways across pooled pseudotime, K27M vs KO (all 5 libraries, KO r1+r2+r3)\n"
             f"shared colour scale; bins {min(nK+nO)}-{max(nK+nO)} cells/genotype", fontsize=11, y=1.13)
fig.subplots_adjust(left=0.17, top=0.80)
fig.savefig(f"{FIG}/pooled_paired_heatmaps_grouped_v2_subset.png", dpi=170, bbox_inches="tight")
print("written")
