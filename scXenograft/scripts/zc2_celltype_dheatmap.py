"""Full metabolic panel (curated KEGG-78 + Reactome/WP + Hallmark) x cell type
effect-size heatmap, K27M vs KO. Same pipeline as the violins: depth-robust
control-bin scores, per-(pathway,cell type) Cohen's d, clean library-separation
flag (all 5 libraries), low-representation flag. Cells are cross-hatched where
low-rep; a dot marks clean separation."""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from itertools import combinations

W = REPO_ROOT + "/scXenograft/work"
OUT = REPO_ROOT + "/scXenograft"; FIG = f"{OUT}/figures"
CLS = ["RGC","GLIP","OPC","OL","ASTRO","EPEN","NEURON"]

d = pd.read_csv(f"{W}/pooled5_dpt.csv.gz", index_col=0)
d["label"] = d[CLS].idxmax(1); d["lib"] = d["sample"]
d["geno"] = np.where(d.genotype == "K27M", "K27M", "KO")
D = pd.read_csv(f"{W}/pooled5_metab_resid.csv.gz", index_col=0).loc[d.index]
paths = list(D.columns)

def prettify(p):
    return (p.replace("KEGG_","").replace("HM_","Hallmark ").replace("CHOL_","chol:")
             .replace("FA_","FA:").replace("CTX_KEGGc_","").replace("CTX_Hallmark_","Hallmark ")
             .replace("_"," "))[:46]

rows = []
for p in paths:
    v = D[p].values
    for ct in CLS:
        m = (d.label == ct).values
        k = v[m & (d.geno=="K27M").values]; o = v[m & (d.geno=="KO").values]
        nk, no = len(k), len(o)
        sub = d[m]
        klib = sub[sub.geno=="K27M"].groupby("lib").size(); olib = sub[sub.geno=="KO"].groupby("lib").size()
        minlib = int(min(klib.min() if len(klib) else 0, olib.min() if len(olib) else 0))
        psd = np.sqrt((np.var(k)+np.var(o))/2) if nk>1 and no>1 else np.nan
        cohen = (np.mean(k)-np.mean(o))/psd if psd else np.nan
        lm = pd.Series(v[m]).groupby(sub.lib.values).mean(); lg = sub.groupby("lib").geno.first()
        km = lm[lg=="K27M"].values; om = lm[lg=="KO"].values
        sep = (len(km) and len(om) and (km.min()>om.max() or om.min()>km.max()))
        low = (min(nk,no)<50) or (minlib<15)
        rows.append(dict(pathway=p, cell_type=ct, cohen_d=cohen, clean=bool(sep), low_rep=low))
R = pd.DataFrame(rows)
R.to_csv(f"{OUT}/celltype_allpathway_stats.csv", index=False)

Dmat = R.pivot(index="pathway", columns="cell_type", values="cohen_d")[CLS]
CLEAN = R.pivot(index="pathway", columns="cell_type", values="clean")[CLS]
LOW = R.pivot(index="pathway", columns="cell_type", values="low_rep")[CLS]
# order pathways by mean d across adequately-repped cell types
score = Dmat.where(~LOW).mean(1)
Dmat = Dmat.loc[score.sort_values(ascending=False).index]
CLEAN = CLEAN.loc[Dmat.index]; LOW = LOW.loc[Dmat.index]

n = len(Dmat)
fig, ax = plt.subplots(figsize=(9, max(14, 0.23*n)))
vmax = 1.2
im = ax.imshow(Dmat.values, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
ax.set_xticks(range(len(CLS))); ax.set_xticklabels(CLS, fontsize=9)
ax.set_yticks(range(n)); ax.set_yticklabels([prettify(p) for p in Dmat.index], fontsize=5.6)
for i in range(n):
    for j in range(len(CLS)):
        if bool(LOW.values[i, j]):
            ax.add_patch(plt.Rectangle((j-0.5, i-0.5), 1, 1, fill=False, hatch="xxx",
                         edgecolor="#888", lw=0, alpha=0.9))
        elif bool(CLEAN.values[i, j]):
            ax.plot(j, i, marker=".", color="k", ms=2.5)
ax.set_title("Metabolic pathway × cell-type effect size (Cohen's d, K27M − KO)\n"
             "all 5 libraries; depth-robust control-bin. dot = all 5 libraries separate; "
             "hatch = low representation", fontsize=10)
fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02, label="Cohen's d  (red K27M>KO, blue KO>K27M)")
fig.tight_layout()
fig.savefig(f"{FIG}/celltype_allpathway_dheatmap.png", dpi=170, bbox_inches="tight")
plt.close(fig)

print(f"{n} pathways x {len(CLS)} cell types")
print("\nTop 12 K27M>KO (by mean d over well-repped cell types):")
print(score.sort_values(ascending=False).head(12).round(2).to_string())
print("\nTop 12 KO>K27M:")
print(score.sort_values().head(12).round(2).to_string())
# consistency summary
cons = R[~R.low_rep].groupby("pathway").agg(
    mean_d=("cohen_d","mean"), n_clean=("clean","sum"),
    n_ct=("clean","size")).sort_values("mean_d")
cons["frac_clean"] = (cons.n_clean/cons.n_ct).round(2)
cons.round(3).to_csv(f"{OUT}/celltype_allpathway_consistency.csv")
print("\nmost consistent K27M>KO (frac cell types with clean separation, mean_d>0):")
print(cons[cons.mean_d>0].sort_values("frac_clean",ascending=False).head(10)[["mean_d","frac_clean"]].to_string())
print("\nmost consistent KO>K27M:")
print(cons[cons.mean_d<0].sort_values("frac_clean",ascending=False).head(10)[["mean_d","frac_clean"]].to_string())
