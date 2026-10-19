"""pooled_manifold_representation_5lib v2 — pseudotime panel as hexbin (mean dpt
per hex) for presentation clarity. Other panels unchanged."""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
W = REPO_ROOT + "/scXenograft/work"
FIG = REPO_ROOT + "/scXenograft/figures"
GT = {"K27M": "#C1272D", "K27M-KO": "#2B6CB0"}
NB = 20
emb = pd.read_csv(f"{W}/pooled5_dpt.csv.gz", index_col=0)
emb["bin"] = pd.qcut(emb.dpt, NB, labels=False, duplicates="drop")
rep = emb.groupby("bin").apply(lambda d: pd.Series({
    "dpt_mid": d.dpt.mean(),
    "n_K27M": int((d.genotype == "K27M").sum()),
    "n_KO": int((d.genotype == "K27M-KO").sum())})).reset_index()
pad = 0.5
ext = [emb.UMAP1.min()-pad, emb.UMAP1.max()+pad, emb.UMAP2.min()-pad, emb.UMAP2.max()+pad]

def scat(ax, x, y, c, **k):
    ax.scatter(x, y, c=c, s=4, linewidths=0, rasterized=True, **k)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_visible(False)

fig, axes = plt.subplots(1, 4, figsize=(19, 4.8))
# genotype (scatter)
for g, sub in emb.groupby("genotype"):
    scat(axes[0], sub.UMAP1, sub.UMAP2, GT[g], label=g)
axes[0].legend(markerscale=4, frameon=False, fontsize=10); axes[0].set_title("genotype", fontsize=12)
# pseudotime HEXBIN (mean dpt per hex)
hb = axes[1].hexbin(emb.UMAP1, emb.UMAP2, C=emb.dpt.values, reduce_C_function=np.mean,
                    gridsize=45, extent=ext, cmap="viridis", mincnt=1,
                    linewidths=0.2, edgecolors="white")
axes[1].set_xlim(ext[0], ext[1]); axes[1].set_ylim(ext[2], ext[3])
axes[1].set_xticks([]); axes[1].set_yticks([])
for sp in axes[1].spines.values(): sp.set_visible(False)
axes[1].set_title("pooled DPT pseudotime (hexbin, mean/hex)", fontsize=12)
cb = fig.colorbar(hb, ax=axes[1], fraction=.045); cb.set_label("pseudotime", fontsize=11)
# density per genotype
for g, sub in emb.groupby("genotype"):
    axes[2].hist(sub.dpt, bins=40, range=(0, 1), density=True, histtype="step", lw=2,
                 color=GT[g], label="KO (r1+r2+r3)" if g == "K27M-KO" else g)
axes[2].set_xlabel("pseudotime", fontsize=12); axes[2].set_ylabel("density", fontsize=12)
axes[2].legend(frameon=False, fontsize=10); axes[2].set_title("density per genotype", fontsize=12)
for sp in ["top", "right"]: axes[2].spines[sp].set_visible(False)
# per-bin representation
ax = axes[3]
ax.bar(rep.dpt_mid, rep.n_K27M, width=.02, color=GT["K27M"], label="K27M")
ax.bar(rep.dpt_mid, -rep.n_KO, width=.02, color=GT["K27M-KO"], label="KO (r1+r2+r3)")
ax.axhline(0, lw=.6, c="k"); ax.axhline(40, lw=.6, ls=":", c="k"); ax.axhline(-40, lw=.6, ls=":", c="k")
ax.set_xlabel("pseudotime", fontsize=12); ax.set_ylabel("cells per bin", fontsize=12)
ax.legend(frameon=False, fontsize=10); ax.set_title("per-bin representation", fontsize=12)
for sp in ["top", "right"]: ax.spines[sp].set_visible(False)
fig.suptitle("HSJ019 pooled K27M + KO manifold — ALL 5 libraries (KO r1+r2+r3); pseudotime as hexbin",
             fontsize=13)
fig.tight_layout(); fig.savefig(f"{FIG}/pooled_manifold_representation_5lib_hex.png", dpi=170, bbox_inches="tight")
print("written")
