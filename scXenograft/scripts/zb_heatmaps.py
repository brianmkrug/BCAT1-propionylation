"""Paired K27M | KO heatmaps: pathways (rows) x pooled pseudotime bins (columns),
common scale. Values = pooled depth-residualised control-bin score, averaged per
(genotype, bin), then each row centred on its pooled (both-genotype) mean so the
between-genotype difference is preserved and baseline offsets don't dominate.
"""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

W = REPO_ROOT + "/scXenograft/work"
OUT = REPO_ROOT + "/scXenograft"; FIG = f"{OUT}/figures"
NB = 20
PATH = ["CHOL_Reactome_cholesterol_biosynthesis", "KEGG_Citrate cycle (TCA cycle)",
        "KEGG_Oxidative phosphorylation", "CHOL_KEGGc_cholesterol_metabolism",
        "KEGG_Sphingolipid metabolism"]
LAB = ["Cholesterol biosynthesis\n(Reactome)", "Citrate cycle (TCA)",
       "Oxidative phosphorylation", "Cholesterol metabolism\n(uptake/transport, KEGG)",
       "Sphingolipid metabolism"]

emb = pd.read_csv(f"{W}/pooled_dpt.csv.gz", index_col=0)
emb["bin"] = pd.qcut(emb.dpt, NB, labels=False, duplicates="drop")
D = pd.read_csv(f"{W}/pooled_metab_resid.csv.gz", index_col=0).loc[emb.index]

# per-(genotype,bin) mean; centre each row on the pooled per-pathway cell mean
def mat(gt):
    m = np.vstack([[np.nanmean(D[p].values[(emb.genotype == gt).values & (emb.bin == b).values])
                    for b in range(NB)] for p in PATH])
    return m
row_mean = np.array([np.nanmean(D[p].values) for p in PATH])[:, None]
K = mat("K27M") - row_mean
O = mat("K27M-KO") - row_mean

vmax = np.nanmax(np.abs(np.concatenate([K, O])))      # shared symmetric scale
nK = [int(((emb.genotype == "K27M") & (emb.bin == b)).sum()) for b in range(NB)]
nO = [int(((emb.genotype == "K27M-KO") & (emb.bin == b)).sum()) for b in range(NB)]

fig = plt.figure(figsize=(15, 4.6))
gs = GridSpec(1, 3, width_ratios=[1, 1, 0.05], wspace=0.08)
for j, (M, ttl, nn) in enumerate([(K, "K27M", nK), (O, "K27M-KO (r1+r3)", nO)]):
    ax = fig.add_subplot(gs[0, j])
    im = ax.imshow(M, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    ax.set_xticks([0, NB//2, NB-1]); ax.set_xticklabels(["early", "mid", "late"])
    ax.set_xlabel("pooled pseudotime bin")
    if j == 0:
        ax.set_yticks(range(len(PATH))); ax.set_yticklabels(LAB, fontsize=8.5)
    else:
        ax.set_yticks([])
    ax.set_title(ttl, fontsize=11)
    for b in range(NB):                     # mark sparse (<80-cell) columns
        if nn[b] < 80:
            ax.add_patch(plt.Rectangle((b-.5, -.5), 1, len(PATH), fill=False,
                         edgecolor="k", lw=0, hatch="///", alpha=.25))
cax = fig.add_subplot(gs[0, 2])
fig.colorbar(im, cax=cax, label="control-bin score,\nrow-centred (common scale)")
fig.suptitle("Metabolic pathways across pooled pseudotime, K27M vs KO (shared colour scale)\n"
             "hatched columns = <80 cells of that genotype (noisier)", fontsize=11)
fig.tight_layout(); fig.savefig(f"{FIG}/pooled_paired_heatmaps.png", dpi=170,
                                bbox_inches="tight")
plt.close(fig)

# also export the underlying matrices
for M, g in [(K, "K27M"), (O, "KO")]:
    pd.DataFrame(M, index=PATH, columns=[f"bin{b}" for b in range(NB)]).round(4).to_csv(
        f"{OUT}/pooled_heatmap_{g}_rowcentred.csv")
print("figure written; shared |max| =", round(float(vmax), 4))
print("per-bin cells  K27M:", nK)
print("per-bin cells  KO:  ", nO)
