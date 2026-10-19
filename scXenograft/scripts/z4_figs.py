import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

W = REPO_ROOT + "/scXenograft/work"
OUT = REPO_ROOT + "/scXenograft"; FIG = f"{OUT}/figures"
GT = {"K27M": "#C1272D", "K27M-KO": "#2B6CB0"}
emb = pd.read_csv(f"{W}/pooled_dpt.csv.gz", index_col=0)
NB = 20
emb["bin"] = pd.qcut(emb.dpt, NB, labels=False, duplicates="drop")
rep = pd.read_csv(f"{W}/pooled_bin_representation.csv", index_col=0)
R = pd.read_csv(f"{OUT}/pooled_pathway_genotype_binconditional.csv", index_col=0)
D = pd.read_csv(f"{W}/pooled_metab_resid.csv.gz", index_col=0).loc[emb.index]

def scat(ax, x, y, c, s=3, **k):
    ax.scatter(x, y, c=c, s=s, linewidths=0, rasterized=True, **k)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_visible(False)

# ---- Fig 1: pooled manifold + representation ----
fig, axes = plt.subplots(1, 4, figsize=(19, 4.6))
for g, sub in emb.groupby("genotype"):
    scat(axes[0], sub.UMAP1, sub.UMAP2, GT[g], label=g)
axes[0].legend(markerscale=4, frameon=False, fontsize=8); axes[0].set_title("genotype", fontsize=10)
o = np.argsort(emb.dpt.values)
sm = axes[1].scatter(emb.UMAP1.values[o], emb.UMAP2.values[o], c=emb.dpt.values[o],
                     s=3, cmap="viridis", linewidths=0, rasterized=True)
axes[1].set_xticks([]); axes[1].set_yticks([])
for sp in axes[1].spines.values(): sp.set_visible(False)
axes[1].set_title("pooled DPT pseudotime", fontsize=10); fig.colorbar(sm, ax=axes[1], fraction=.045)
for g, sub in emb.groupby("genotype"):
    axes[2].hist(sub.dpt, bins=40, range=(0, 1), density=True, histtype="step",
                 lw=2, color=GT[g], label=g)
axes[2].set_xlabel("pseudotime"); axes[2].set_ylabel("density"); axes[2].legend(frameon=False, fontsize=8)
axes[2].set_title("density per genotype\n(the compositional bias)", fontsize=10)
for sp in ["top", "right"]: axes[2].spines[sp].set_visible(False)
ax = axes[3]
ax.bar(rep.dpt_mid, rep.n_K27M, width=.02, color=GT["K27M"], label="K27M")
ax.bar(rep.dpt_mid, -rep.n_KO, width=.02, color=GT["K27M-KO"], label="KO (r1+r3)")
ax.axhline(0, lw=.6, c="k"); ax.axhline(40, lw=.6, ls=":", c="k"); ax.axhline(-40, lw=.6, ls=":", c="k")
ax.set_xlabel("pseudotime"); ax.set_ylabel("cells per bin"); ax.legend(frameon=False, fontsize=8)
ax.set_title("per-bin representation\n(dotted = n=40 floor; both met everywhere)", fontsize=10)
for sp in ["top", "right"]: ax.spines[sp].set_visible(False)
fig.suptitle("HSJ019 pooled K27M + KO manifold with a common pseudotime axis", fontsize=12)
fig.tight_layout(); fig.savefig(f"{FIG}/pooled_manifold_representation.png", dpi=160); plt.close(fig)

# ---- Fig 2: naive vs bin-adjusted delta ----
fig, ax = plt.subplots(figsize=(7.2, 6.8))
ax.scatter(R.naive_delta, R.binadj_delta, s=22, c="#888888")
lim = max(R.naive_delta.abs().max(), R.binadj_delta.abs().max())*1.1
ax.plot([-lim, lim], [-lim, lim], ls="--", lw=.8, c="k")
ax.axhline(0, lw=.5, c="k"); ax.axvline(0, lw=.5, c="k")
lab = set(R.binadj_delta.nlargest(6).index) | set(R.binadj_delta.nsmallest(6).index)
for p in lab:
    ax.annotate(p.replace("KEGG_","").replace("HM_","HM ").replace("CHOL_","").replace("CTX_","")[:30],
                (R.naive_delta[p], R.binadj_delta[p]), fontsize=6.5, xytext=(3,2), textcoords="offset points")
ax.set_xlabel("naive genotype delta (marginal, ignores pseudotime)")
ax.set_ylabel("bin-adjusted genotype delta (conditioned on pseudotime)")
ax.set_title("Conditioning on pseudotime barely changes the genotype comparison\n"
             "(Spearman 0.97, 3/108 flip) -> composition is not the confound here", fontsize=10)
for sp in ["top", "right"]: ax.spines[sp].set_visible(False)
fig.tight_layout(); fig.savefig(f"{FIG}/pooled_naive_vs_binadjusted.png", dpi=160); plt.close(fig)

# ---- Fig 3: genotype profiles along shared pseudotime for key pathways ----
KEY = ["CHOL_Reactome_cholesterol_biosynthesis","KEGG_Terpenoid backbone biosynthesis",
       "HM_MTORC1_SIGNALING","KEGG_Citrate cycle (TCA cycle)",
       "KEGG_Oxidative phosphorylation","CHOL_KEGGc_cholesterol_metabolism",
       "FA_KEGGc_unsaturated_FA_biosynthesis","KEGG_Sphingolipid metabolism"]
mids = rep.set_index(rep.index).dpt_mid
fig, axes = plt.subplots(2, 4, figsize=(18, 8))
for ax, c in zip(axes.ravel(), KEY):
    v = D[c].values
    for g in ["K27M", "K27M-KO"]:
        prof = [np.nanmean(v[(emb.genotype == g).values & (emb.bin == b).values])
                for b in range(NB)]
        ax.plot(mids.reindex(range(NB)).values, prof, lw=2.2, color=GT[g], label=g)
    bd = R.binadj_delta[c]
    ax.set_title(c.replace("KEGG_","").replace("HM_","HM ").replace("CHOL_","").replace("CTX_","")[:34]
                 + f"\nbinadj Δ(K27M-KO)={bd:+.3f}", fontsize=8.5)
    ax.set_xlabel("shared pseudotime"); ax.axhline(0, lw=.5, c="k")
    for sp in ["top", "right"]: ax.spines[sp].set_visible(False)
axes.ravel()[0].legend(frameon=False, fontsize=8)
fig.suptitle("Metabolic pathways along the SHARED pseudotime axis, by genotype "
             "(depth-robust, common scale)", fontsize=12)
fig.tight_layout(); fig.savefig(f"{FIG}/pooled_pathway_profiles_by_genotype.png", dpi=160); plt.close(fig)
print("figures written")
print("\nBiosynthesis vs uptake split:")
for c in ["CHOL_Reactome_cholesterol_biosynthesis","KEGG_Terpenoid backbone biosynthesis",
          "CHOL_KEGGc_cholesterol_metabolism"]:
    print(f"  {c[:44]:<44} binadj Δ {R.binadj_delta[c]:+.3f}  "
          f"({'K27M' if R.binadj_delta[c]>0 else 'KO'} higher)")
