"""v2 of pooled_pathway_profiles_by_genotype (8 key pathways along shared axis),
using all 5 libraries (KO r1+r2+r3)."""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
W = REPO_ROOT + "/scXenograft/work"
FIG = REPO_ROOT + "/scXenograft/figures"
GT = {"K27M": "#C1272D", "K27M-KO": "#2B6CB0"}
NB = 20
KEY = ["CHOL_Reactome_cholesterol_biosynthesis", "KEGG_Terpenoid backbone biosynthesis",
       "HM_MTORC1_SIGNALING", "KEGG_Citrate cycle (TCA cycle)",
       "KEGG_Oxidative phosphorylation", "CHOL_KEGGc_cholesterol_metabolism",
       "FA_KEGGc_unsaturated_FA_biosynthesis", "KEGG_Sphingolipid metabolism"]

d = pd.read_csv(f"{W}/pooled5_dpt.csv.gz", index_col=0)
D = pd.read_csv(f"{W}/pooled5_metab_resid.csv.gz", index_col=0).loc[d.index]
d["bin"] = pd.qcut(d.dpt, NB, labels=False, duplicates="drop")
mids = d.groupby("bin").dpt.mean().reindex(range(NB)).values
isK = (d.genotype == "K27M").values

fig, axes = plt.subplots(2, 4, figsize=(18, 8))
for ax, c in zip(axes.ravel(), KEY):
    v = D[c].values
    kb = np.array([np.nanmean(v[isK & (d.bin == b).values]) for b in range(NB)])
    ob = np.array([np.nanmean(v[~isK & (d.bin == b).values]) for b in range(NB)])
    ax.plot(mids, kb, lw=2.2, color=GT["K27M"], label="K27M")
    ax.plot(mids, ob, lw=2.2, color=GT["K27M-KO"], label="KO (r1+r2+r3)")
    binadj = np.nanmean(kb - ob)
    ax.set_title(c.replace("KEGG_", "").replace("HM_", "HM ").replace("CHOL_", "")
                 .replace("CTX_", "")[:34] + f"\nbinadj Δ(K27M-KO)={binadj:+.3f}", fontsize=8.5)
    ax.set_xlabel("shared pseudotime"); ax.axhline(0, lw=.5, c="k")
    for sp in ["top", "right"]: ax.spines[sp].set_visible(False)
axes.ravel()[0].legend(frameon=False, fontsize=8)
fig.suptitle("v2 — Metabolic pathways along the SHARED pseudotime axis, by genotype "
             "(all 5 libraries, KO r1+r2+r3; depth-robust, common scale)", fontsize=12)
fig.tight_layout(); fig.savefig(f"{FIG}/pooled_pathway_profiles_by_genotype_v2.png", dpi=160)
print("profiles v2 written")
