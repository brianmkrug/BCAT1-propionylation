"""v2 of pooled_paired_heatmaps_grouped, all 5 libraries (KO r1+r2+r3)."""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import numpy as np, pandas as pd, scipy.sparse as sp, scanpy as sc
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
W = REPO_ROOT + "/scXenograft/work"
FIG = REPO_ROOT + "/scXenograft/figures"
NB = 20
d = pd.read_csv(f"{W}/pooled5_dpt.csv.gz", index_col=0)
d["bin"] = pd.qcut(d.dpt, NB, labels=False, duplicates="drop")
D = pd.read_csv(f"{W}/pooled5_metab_resid.csv.gz", index_col=0).loc[d.index]

# proliferation via same control-bin pipeline
S = """MCM5 PCNA TYMS FEN1 MCM2 MCM4 RRM1 UNG GINS2 MCM6 CDCA7 DTL PRIM1 UHRF1 MLF1IP
HELLS RFC2 RPA2 NASP RAD51AP1 GMNN WDR76 SLBP CCNE2 UBR7 POLD3 MSH2 ATAD2 RAD51 RRM2
CDC45 CDC6 EXO1 TIPIN DSCC1 BLM CASP8AP2 USP1 CLSPN POLA1 CHAF1B BRIP1 E2F8""".split()
G2M = """HMGB2 CDK1 NUSAP1 UBE2C BIRC5 TPX2 TOP2A NDC80 CKS2 NUF2 CKS1B MKI67 TMPO
CENPF TACC3 FAM64A SMC4 CCNB2 CKAP2L CKAP2 AURKB BUB1 KIF11 ANP32E TUBB4B GTSE1 KIF20B
HJURP CDCA3 HN1 CDC20 TTK CDC25C KIF2C RANGAP1 NCAPD2 DLGAP5 CDCA2 CDCA8 ECT2 KIF23
HMMR AURKA PSRC1 ANLN LBR CKAP5 CENPE CTCF NEK2 G2E3 GAS2L3 CBX5 CENPA""".split()
A = sc.read_h5ad(f"{W}/hsj019_norm.h5ad")[d.index]
X = A.X.tocsc() if sp.issparse(A.X) else sp.csc_matrix(A.X)
genes = pd.Index(A.var_names); NG = X.shape[1]
avg = np.asarray(X.mean(0)).ravel(); order = np.argsort(avg)
bg = np.zeros(NG, int); bg[order] = np.minimum((np.arange(NG)/NG*25).astype(int), 24)
members = {b: np.where(bg == b)[0] for b in range(25)}
rng = np.random.default_rng(0)
idx = np.unique(genes.get_indexer(S+G2M)); idx = idx[idx >= 0]
setv = np.asarray(X[:, idx].mean(1)).ravel()
ctrl = []
for b, cnt in zip(*np.unique(bg[idx], return_counts=True)):
    pool = np.setdiff1d(members[b], idx); take = min(len(pool), cnt*50)
    if take: ctrl.append(rng.choice(pool, take, replace=False))
ctrl = np.concatenate(ctrl)
prolif = setv - np.asarray(X[:, ctrl].mean(1)).ravel()
ng = d.n_genes.values.astype(float); M = np.c_[np.ones_like(ng), ng, ng**2]
D["PROLIFERATION"] = prolif - M @ np.linalg.lstsq(M, prolif, rcond=None)[0]

GROUPS = [
 ("K27M-enriched", [("CHOL_Reactome_cholesterol_biosynthesis","Cholesterol biosynthesis (Reactome)"),
                    ("KEGG_Citrate cycle (TCA cycle)","Citrate cycle (TCA)")]),
 ("KO-enriched", [("KEGG_Oxidative phosphorylation","Oxidative phosphorylation"),
                  ("CHOL_KEGGc_cholesterol_metabolism","Cholesterol metabolism (uptake/transport)"),
                  ("KEGG_Sphingolipid metabolism","Sphingolipid metabolism")]),
 ("Common", [("PROLIFERATION","Proliferation (S+G2M)"),
             ("KEGG_Purine metabolism","Purine metabolism"),
             ("KEGG_One carbon pool by folate","One carbon pool by folate")]),
]
paths = [p for _, rows in GROUPS for p, _ in rows]
labels = [l for _, rows in GROUPS for _, l in rows]
gsizes = [(n, len(r)) for n, r in GROUPS]
row_mean = np.array([np.nanmean(D[p].values) for p in paths])[:, None]
def mat(gt):
    return np.vstack([[np.nanmean(D[p].values[(d.genotype==gt).values & (d.bin==b).values])
                       for b in range(NB)] for p in paths])
K = mat("K27M") - row_mean; O = mat("K27M-KO") - row_mean
vmax = float(np.nanpercentile(np.abs(np.r_[K.ravel(), O.ravel()]), 98))
nK = [int(((d.genotype=="K27M")&(d.bin==b)).sum()) for b in range(NB)]
nO = [int(((d.genotype=="K27M-KO")&(d.bin==b)).sum()) for b in range(NB)]
GCOL = {"K27M-enriched":"#C1272D","KO-enriched":"#2B6CB0","Common":"#555555"}

fig = plt.figure(figsize=(15.5, 5.6))
gs = GridSpec(1, 3, width_ratios=[1,1,0.045], wspace=0.06)
for j,(Mx,ttl) in enumerate([(K,"K27M"),(O,"K27M-KO (r1+r2+r3)")]):
    ax = fig.add_subplot(gs[0,j])
    im = ax.imshow(Mx, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    ax.set_xticks([0,NB//2,NB-1]); ax.set_xticklabels(["early","mid","late"]); ax.set_xlabel("pooled pseudotime bin")
    if j==0: ax.set_yticks(range(len(paths))); ax.set_yticklabels(labels, fontsize=8.5)
    else: ax.set_yticks([])
    ax.set_title(ttl, fontsize=11)
    off=0
    for _,sz in gsizes:
        if off>0: ax.axhline(off-0.5, color="k", lw=1.6)
        off+=sz
fig.colorbar(im, cax=fig.add_subplot(gs[0,2]), label="control-bin score, row-centred\n(shared scale; ±98th pct)")
off=0
for name,sz in gsizes:
    yc=(off+off+sz-1)/2
    fig.text(0.012, 1-(0.13+0.74*(yc+0.5)/len(paths)), name, rotation=90, va="center", ha="center",
             fontsize=10, fontweight="bold", color=GCOL[name]); off+=sz
fig.suptitle("v2 — Metabolic pathways across pooled pseudotime, K27M vs KO (all 5 libraries, KO r1+r2+r3)\n"
             f"shared colour scale; bins {min(nK+nO)}-{max(nK+nO)} cells/genotype", fontsize=11)
fig.subplots_adjust(left=0.16)
fig.savefig(f"{FIG}/pooled_paired_heatmaps_grouped_v2.png", dpi=170, bbox_inches="tight")
print("grouped heatmap v2 written; bin cell range K27M", min(nK), max(nK), "KO", min(nO), max(nO))
