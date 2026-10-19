"""BCAA heatmap (5 lib) — ACADSB force-included, clustered, ECHS1/HSD17B10 moved
below LDHB, 3 groups rendered as SEPARATE outlined objects with white gaps.
Green/purple highlighted labels. Pooled z per gene; log2FC strip."""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import numpy as np, pandas as pd, scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.cluster.hierarchy import linkage, leaves_list

B = REPO_ROOT + ""; W = f"{B}/scXenograft/work"
OUT = f"{B}/scXenograft"; FIG = f"{OUT}/figures"
NB = 20; DET = 12; FORCE = set(); ALIAS = {"MMUT": "MUT"}   # ACADSB dropped (sparse detection)
GREEN = {"BCAT1","ACAD8","ECHS1","HADHA","HSD17B10","ACAT1"}
PURPLE = {"PCCA","PCCB","MUT","MCEE"}
GCOL, PCOL = "#1a7d1a", "#6a1b9a"
GENES = ["BCAT1","BCAT2","BCKDHA","BCKDHB","DBT","DLD","IVD","ACAD8","ACADSB","ACADM","ACADS",
 "MCCC1","MCCC2","AUH","HMGCL","HMGCLL1","HMGCS1","OXCT1","OXCT2","ACAT1","ACAT2",
 "ECHS1","HIBCH","HIBADH","ALDH6A1","ABAT","HADH","HADHA","HADHB","HSD17B10","ACAA1","ACAA2",
 "EHHADH","AACS","ALDH2","ALDH3A2","ALDH7A1","ALDH9A1","ALDH1B1","AOX1","IL4I1","ACSF3","AGXT2",
 "PCCA","PCCB","MCEE","MMUT","SUCLG1","SUCLG2","SUCLA2","ACSS1","ACSS2","ACSS3","ACACA","ACACB",
 "MLYCD","ACOX1","ACOX3","ECHDC1","LDHA","LDHB","LDHC"]

emb = pd.read_csv(f"{W}/pooled5_dpt.csv.gz", index_col=0)
emb["bin"] = pd.qcut(emb.dpt, NB, labels=False, duplicates="drop")
A = sc.read_h5ad(f"{W}/hsj019_norm.h5ad")[emb.index]
de = pd.read_csv(f"{OUT}/sc_pseudobulk_deseq2_K27M_vs_KO.csv", index_col=0)
vset = set(A.var_names)
def ref(g): return ALIAS.get(g, g)
def det(g):
    r = ref(g); return 100*(np.asarray(A[:, r].X.todense()).ravel() > 0).mean() if r in vset else 0
ordered, seen = [], set()
for g in GENES:
    r = ref(g)
    if r in vset and r not in seen and (det(g) >= DET or g in FORCE):
        ordered.append(r); seen.add(r)
E = pd.DataFrame(np.asarray(A[:, ordered].X.todense()), index=emb.index, columns=ordered)
def zmat(gt):
    return np.vstack([[E[g].values[(emb.genotype==gt).values & (emb.bin==b).values].mean()
                       for b in range(NB)] for g in ordered])
Kr, Or = zmat("K27M"), zmat("K27M-KO")
mu = np.array([np.r_[Kr[i],Or[i]].mean() for i in range(len(ordered))])[:,None]
sd = np.array([np.r_[Kr[i],Or[i]].std() for i in range(len(ordered))])[:,None]
Kz, Oz = (Kr-mu)/sd, (Or-mu)/sd
val = pd.DataFrame(np.hstack([Kz, Oz]), index=ordered,
                   columns=[f"K27M_{b}" for b in range(NB)]+[f"KO_{b}" for b in range(NB)])
val["lfc"] = [de.log2FoldChange.get(g, np.nan) for g in ordered]

clustered = list(val.index[leaves_list(linkage(np.nan_to_num(val.iloc[:, :2*NB].values),
                                                method="average", metric="correlation"))])
order = [g for g in clustered if g not in ("ECHS1", "HSD17B10")]
i = order.index("LDHB"); order = order[:i+1] + ["ECHS1", "HSD17B10"] + order[i+1:]
b1 = order.index("ALDH3A2") + 1; b2 = order.index("HSD17B10") + 1
groups = [order[:b1], order[b1:b2], order[b2:]]

VZ = 2.5
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
def rows(gl, cols): return val.loc[gl, cols].values
fig = plt.figure(figsize=(14.5*0.70, max(11, 0.28*len(order))))   # 30% narrower
gs = GridSpec(3, 3, width_ratios=[1, 1, 0.06],
              height_ratios=[len(g) for g in groups], hspace=0.06, wspace=0.05)
im = None
for gi, grp in enumerate(groups):
    Kc = rows(grp, [f"K27M_{b}" for b in range(NB)])
    Oc = rows(grp, [f"KO_{b}" for b in range(NB)])
    Lc = val.loc[grp, "lfc"].values[:, None]
    axK = fig.add_subplot(gs[gi, 0]); axO = fig.add_subplot(gs[gi, 1]); axL = fig.add_subplot(gs[gi, 2])
    im = axK.imshow(Kc, aspect="auto", cmap="RdBu_r", vmin=-VZ, vmax=VZ)
    axO.imshow(Oc, aspect="auto", cmap="RdBu_r", vmin=-VZ, vmax=VZ)
    axL.imshow(Lc, aspect="auto", cmap="PuOr_r", vmin=-1.5, vmax=1.5)
    for ax in (axK, axO, axL):                       # thin black outline on each object
        for sp in ax.spines.values(): sp.set_visible(True); sp.set_edgecolor("k"); sp.set_linewidth(0.9)
    axK.set_yticks(range(len(grp)))
    axK.set_yticklabels(["MMUT/MUT" if g == "MUT" else g for g in grp], fontsize=11)
    for t, g in zip(axK.get_yticklabels(), grp):
        if g in GREEN: t.set_color(GCOL); t.set_fontweight("bold")
        elif g in PURPLE: t.set_color(PCOL); t.set_fontweight("bold")
    axO.set_yticks([]); axL.set_xticks([]); axL.set_yticks([])
    bottom = (gi == len(groups)-1)
    for ax in (axK, axO):
        ax.set_xticks([0, NB//2, NB-1])
        ax.set_xticklabels(["early", "mid", "late"] if bottom else [], fontsize=13)
        if bottom: ax.set_xlabel("pooled pseudotime bin", fontsize=15)
    if gi == 0:
        axK.set_title("K27M", fontsize=16); axO.set_title("K27M-KO (r1+r2+r3)", fontsize=16)
# two equal-height colorbars in the right margin: pooled z (top), log2FC (below)
CBH, CBW, CBL = 0.15, 0.020, 0.90
cax1 = fig.add_axes([CBL, 0.55, CBW, CBH])
cb1 = fig.colorbar(im, cax=cax1); cb1.set_label("pooled z (per gene)", fontsize=14)
cb1.ax.tick_params(labelsize=12)
cax2 = fig.add_axes([CBL, 0.35, CBW, CBH])
sm = ScalarMappable(norm=Normalize(-1.5, 1.5), cmap="PuOr_r")
cb2 = fig.colorbar(sm, cax=cax2); cb2.set_label("log2FC (K27M/KO)", fontsize=14)
cb2.ax.tick_params(labelsize=12)
fig.suptitle("BCAA metabolism genes (all 5 libraries), K27M vs KO across pooled pseudotime\n"
             "clustering order (ECHS1/HSD17B10 moved below LDHB); 3 separate groups; strip = pseudobulk log2FC",
             fontsize=14, y=0.995)
fig.subplots_adjust(left=0.14, right=0.87, top=0.95)
fig.savefig(f"{FIG}/bcaa_gene_heatmap_clustered_5lib_reordered.png", dpi=170, bbox_inches="tight")
print(f"{len(order)} genes; ACADSB in group",
      [gi for gi, g in enumerate(groups) if "ACADSB" in g], "; group sizes", [len(g) for g in groups])
