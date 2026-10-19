"""Stage 8: diagnostic figures - UMAP by genotype/replicate/label, meta-score maps,
replicate-resolved composition. Usage: python3 s8_figs.py bt245|hsj019"""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import sys, json
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

W = REPO_ROOT + "/scXenograft/work"
FIG = REPO_ROOT + "/scXenograft/figures"
import os; os.makedirs(FIG, exist_ok=True)
tag = sys.argv[1]
d = pd.read_csv(f"{W}/{tag}_annotated.csv.gz", index_col=0)
h = pd.read_csv(f"{W}/{tag}_embed_harmony.csv.gz", index_col=0).loc[d.index]
CLS = ["RGC", "GLIP", "OPC", "OL", "ASTRO", "EPEN", "NEURON", "NON_NEURAL"]
GT = {"K27M": "#C1272D", "K27M-KO": "#2B6CB0"}

def scat(ax, x, y, c, s=3, **kw):
    ax.scatter(x, y, c=c, s=s, linewidths=0, rasterized=True, **kw)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_visible(False)

# --- Fig 1: uncorrected vs harmony, by genotype and by replicate ---
fig, axes = plt.subplots(2, 3, figsize=(13.5, 9))
for row, (emb, ttl) in enumerate([(d, "uncorrected merge"), (h, "Harmony on sample")]):
    ax = axes[row, 0]
    for g, sub in emb.groupby("genotype"):
        scat(ax, sub.UMAP1, sub.UMAP2, GT[g], label=g)
    ax.set_title(f"{ttl}\ngenotype", fontsize=10)
    if row == 0: ax.legend(markerscale=4, frameon=False, fontsize=8)
    ax = axes[row, 1]
    cols = plt.cm.tab10(np.linspace(0, 1, emb["sample"].nunique()))
    for (s, sub), col in zip(emb.groupby("sample"), cols):
        scat(ax, sub.UMAP1, sub.UMAP2, [col], label=s.split("_", 1)[1])
    ax.set_title("replicate / library", fontsize=10)
    if row == 0: ax.legend(markerscale=4, frameon=False, fontsize=6)
    ax = axes[row, 2]
    lb = d.label
    cols = dict(zip(sorted(lb.unique()), plt.cm.Set2(np.linspace(0, 1, lb.nunique()))))
    for l, sub in emb.groupby(lb):
        scat(ax, sub.UMAP1, sub.UMAP2, [cols[l]], label=l)
    ax.set_title("atlas label (from uncorrected clustering)", fontsize=10)
    if row == 0: ax.legend(markerscale=4, frameon=False, fontsize=7)
fig.suptitle(f"{tag.upper()} — embedding diagnostics", fontsize=12)
fig.tight_layout(); fig.savefig(f"{FIG}/{tag}_umap_diagnostics.png", dpi=170); plt.close(fig)

# --- Fig 2: meta-score maps on uncorrected UMAP ---
fig, axes = plt.subplots(2, 4, figsize=(16, 7.5))
for ax, c in zip(axes.ravel(), CLS):
    v = d[c].values
    lo, hi = np.percentile(v, [2, 98])
    o = np.argsort(v)
    sm = ax.scatter(d.UMAP1.values[o], d.UMAP2.values[o], c=np.clip(v[o], lo, hi),
                    s=3, cmap="RdYlBu_r", linewidths=0, rasterized=True)
    ax.set_title(c, fontsize=10); ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_visible(False)
    fig.colorbar(sm, ax=ax, fraction=0.045)
fig.suptitle(f"{tag.upper()} — Jessa 2022 atlas meta-scores (uncorrected UMAP)", fontsize=12)
fig.tight_layout(); fig.savefig(f"{FIG}/{tag}_metascore_umaps.png", dpi=170); plt.close(fig)

# --- Fig 3: replicate-resolved composition + meta-score distributions ---
comp = pd.crosstab(d["sample"], d.label, normalize="index").mul(100)
fig, axes = plt.subplots(1, 2, figsize=(14, 4.6),
                         gridspec_kw={"width_ratios": [1.1, 1.5]})
lbls = sorted(d.label.unique())
cols = plt.cm.Set2(np.linspace(0, 1, len(lbls)))
bot = np.zeros(len(comp))
for l, c in zip(lbls, cols):
    if l in comp:
        axes[0].bar(range(len(comp)), comp[l], bottom=bot, color=c, label=l, width=.7)
        bot += comp[l].values
axes[0].set_xticks(range(len(comp)))
axes[0].set_xticklabels([i.split("_", 1)[1] for i in comp.index], rotation=45, ha="right", fontsize=8)
axes[0].set_ylabel("% of cells"); axes[0].legend(fontsize=7, frameon=False,
                                                 bbox_to_anchor=(1, 1), loc="upper left")
axes[0].set_title("composition per library", fontsize=10)

pos, ticks = 0, []
for c in CLS[:6]:
    for s, sub in d.groupby("sample"):
        g = sub.genotype.iloc[0]
        axes[1].violinplot(sub[c].values, positions=[pos], widths=.8, showextrema=False)
        axes[1].scatter([pos], [sub[c].median()], color=GT[g], s=14, zorder=3)
        pos += 1
    ticks.append(pos - d["sample"].nunique()/2 - .5); pos += 1.2
axes[1].set_xticks(ticks); axes[1].set_xticklabels(CLS[:6], fontsize=9)
axes[1].set_ylabel("meta-score"); axes[1].axhline(0, lw=.6, color="k")
axes[1].set_title("meta-score per library (dot = median; red K27M / blue KO)", fontsize=10)
for sp in ["top", "right"]: axes[1].spines[sp].set_visible(False)
fig.suptitle(f"{tag.upper()} — replicate-resolved composition and lineage scores", fontsize=12)
fig.tight_layout(); fig.savefig(f"{FIG}/{tag}_composition_by_replicate.png", dpi=170); plt.close(fig)

comp.to_csv(f"{FIG}/{tag}_composition_pct.csv")
d.groupby("sample")[CLS + ["PROLIF_sig"]].median().round(3).to_csv(f"{FIG}/{tag}_metascore_medians.csv")
print(tag, "figures written")
print(comp.round(1).to_string())
