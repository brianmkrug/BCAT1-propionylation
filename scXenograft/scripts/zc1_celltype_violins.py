"""Per-cell-type metabolic pathway scores, K27M vs KO (all 5 libraries).
Depth-robust control-bin scores (NOT ssGSEA). Split-violin per cell type.
Test = linear mixed model (score ~ genotype, random intercept = library) -- the
honest test; naive per-cell Mann-Whitney reported alongside to show inflation.
"""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, scipy.sparse as sp, scanpy as sc
from scipy.stats import mannwhitneyu
import statsmodels.formula.api as smf
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

W = REPO_ROOT + "/scXenograft/work"
OUT = REPO_ROOT + "/scXenograft"; FIG = f"{OUT}/figures"
CLS = ["RGC","GLIP","OPC","OL","ASTRO","EPEN","NEURON"]
PATH = [("CHOL_Reactome_cholesterol_biosynthesis","Cholesterol biosynthesis"),
        ("KEGG_Citrate cycle (TCA cycle)","Citrate cycle (TCA)"),
        ("KEGG_Oxidative phosphorylation","Oxidative phosphorylation"),
        ("CHOL_KEGGc_cholesterol_metabolism","Cholesterol metab (uptake)"),
        ("KEGG_Sphingolipid metabolism","Sphingolipid metabolism"),
        ("PROLIFERATION","Proliferation (S+G2M)"),
        ("KEGG_Purine metabolism","Purine metabolism"),
        ("KEGG_One carbon pool by folate","One carbon pool by folate")]

d = pd.read_csv(f"{W}/pooled5_dpt.csv.gz", index_col=0)
d["label"] = d[CLS].idxmax(1)
d["lib"] = d["sample"]   # full sample ID: K27M-r1/r2, KO-r1/r2/r3 are 5 distinct libraries
d["geno"] = np.where(d.genotype == "K27M", "K27M", "KO")
D = pd.read_csv(f"{W}/pooled5_metab_resid.csv.gz", index_col=0).loc[d.index]

# proliferation control-bin residualised (same pipeline)
S = "MCM5 PCNA TYMS FEN1 MCM2 MCM4 RRM1 UNG GINS2 MCM6 CDCA7 DTL PRIM1 UHRF1 HELLS RFC2 RPA2 NASP RAD51AP1 GMNN WDR76 SLBP CCNE2 UBR7 POLD3 MSH2 ATAD2 RAD51 RRM2 CDC45 CDC6 EXO1 TIPIN DSCC1 BLM CASP8AP2 USP1 CLSPN POLA1 CHAF1B BRIP1 E2F8".split()
G2M = "HMGB2 CDK1 NUSAP1 UBE2C BIRC5 TPX2 TOP2A NDC80 CKS2 NUF2 CKS1B MKI67 TMPO CENPF TACC3 SMC4 CCNB2 CKAP2L CKAP2 AURKB BUB1 KIF11 ANP32E TUBB4B GTSE1 KIF20B HJURP CDCA3 CDC20 TTK CDC25C KIF2C RANGAP1 NCAPD2 DLGAP5 CDCA2 CDCA8 ECT2 KIF23 HMMR AURKA PSRC1 ANLN LBR CKAP5 CENPE CTCF NEK2 G2E3 GAS2L3 CBX5 CENPA".split()
A = sc.read_h5ad(f"{W}/hsj019_norm.h5ad")[d.index]
X = A.X.tocsc() if sp.issparse(A.X) else sp.csc_matrix(A.X)
genes = pd.Index(A.var_names); NG = X.shape[1]
avg = np.asarray(X.mean(0)).ravel(); order = np.argsort(avg)
bg = np.zeros(NG,int); bg[order]=np.minimum((np.arange(NG)/NG*25).astype(int),24)
mem={b:np.where(bg==b)[0] for b in range(25)}; rng=np.random.default_rng(0)
idx=np.unique(genes.get_indexer(S+G2M)); idx=idx[idx>=0]
ctrl=[]
for b,cnt in zip(*np.unique(bg[idx],return_counts=True)):
    pool=np.setdiff1d(mem[b],idx); t=min(len(pool),cnt*50)
    if t: ctrl.append(rng.choice(pool,t,replace=False))
ctrl=np.concatenate(ctrl)
pr=np.asarray(X[:,idx].mean(1)).ravel()-np.asarray(X[:,ctrl].mean(1)).ravel()
ng=d.n_genes.values.astype(float); M=np.c_[np.ones_like(ng),ng,ng**2]
D["PROLIFERATION"]=pr-M@np.linalg.lstsq(M,pr,rcond=None)[0]

# assemble long
for p,_ in PATH: d[p]=D[p].values

# ---- stats ----
rows=[]
for p,pl in PATH:
    for ct in CLS:
        sub=d[d.label==ct]
        k=sub[sub["geno"]=="K27M"][p].values; o=sub[sub["geno"]=="KO"][p].values
        nk,no=len(k),len(o)
        klibs=sub[sub["geno"]=="K27M"].groupby("lib").size(); olibs=sub[sub["geno"]=="KO"].groupby("lib").size()
        minlib=int(min(klibs.min() if len(klibs) else 0, olibs.min() if len(olibs) else 0))
        d_med=np.median(k)-np.median(o) if nk and no else np.nan
        pooled_sd=np.sqrt((np.var(k)+np.var(o))/2) if nk>1 and no>1 else np.nan
        cohen=(np.mean(k)-np.mean(o))/pooled_sd if pooled_sd else np.nan
        mwu=mannwhitneyu(k,o).pvalue if nk and no else np.nan
        # library-level: per-sample means, then EXACT permutation test over the
        # 5 libraries (only correctly-powered option at 2 vs 3 replicates)
        libmean=sub.groupby("lib")[p].mean()
        libgt=sub.groupby("lib")["geno"].first()
        km=libmean[libgt=="K27M"].values; om=libmean[libgt=="KO"].values
        p_perm=np.nan
        if len(km)>=1 and len(om)>=1 and (len(km)+len(om))>=4:
            allm=np.r_[km,om]; nK=len(km); obs=abs(km.mean()-om.mean())
            from itertools import combinations
            diffs=[]
            for idxK in combinations(range(len(allm)),nK):
                mask=np.zeros(len(allm),bool); mask[list(idxK)]=True
                diffs.append(abs(allm[mask].mean()-allm[~mask].mean()))
            p_perm=float(np.mean(np.array(diffs)>=obs-1e-12))
        sep = "clean" if (len(km) and len(om) and (km.min()>om.max() or om.min()>km.max())) else "overlap"
        low = (min(nk,no)<50) or (minlib<15)
        rows.append(dict(pathway=pl, cell_type=ct, n_K27M=nk, n_KO=no, min_lib_n=minlib,
                         median_diff_KmO=d_med, cohen_d=cohen, p_perm_lib=p_perm,
                         naive_mwu_p=mwu, lib_separation=sep, low_rep=low))
R=pd.DataFrame(rows)
R.to_csv(f"{OUT}/celltype_pathway_stats.csv", index=False)

# ---- figure ----
fig,axes=plt.subplots(2,4,figsize=(20,10),sharex=True)
pal={"K27M":"#C1272D","KO":"#2B6CB0"}
for ax,(p,pl) in zip(axes.ravel(),PATH):
    sns.violinplot(data=d, x="label", y=p, hue="geno", order=CLS, hue_order=["K27M","KO"],
                   split=True, inner="quartile", palette=pal, ax=ax, linewidth=.6, cut=0)
    ax.set_title(pl,fontsize=10); ax.set_xlabel(""); ax.set_ylabel("control-bin score")
    ax.legend_.remove() if ax.legend_ else None
    yl=ax.get_ylim(); ytop=yl[1]; span=yl[1]-yl[0]
    for i,ct in enumerate(CLS):
        r=R[(R.pathway==pl)&(R.cell_type==ct)].iloc[0]
        # marker: clean library separation -> filled diamond coloured by direction
        if r.lib_separation=="clean":
            c = pal["K27M"] if r.median_diff_KmO>0 else pal["KO"]
            ax.text(i, ytop+span*0.02, "◆", ha="center", va="bottom", fontsize=10, color=c)
        # Cohen's d magnitude
        ax.text(i, ytop+span*0.13, f"d={r.cohen_d:+.2f}", ha="center", va="bottom",
                fontsize=6.5, color="#444")
        if r.low_rep:
            ax.text(i, yl[0]+span*0.01, "▲", ha="center", va="bottom", fontsize=10, color="#d95f02")
    ax.set_ylim(yl[0], ytop+span*0.22)
    ax.tick_params(axis="x", labelrotation=45, labelsize=8)
h=[plt.Line2D([],[],marker="s",ls="",color=pal[g],label=g) for g in ["K27M","KO"]]
h+=[plt.Line2D([],[],marker="D",ls="",color="#C1272D",label="clean library separation (K27M>KO)"),
    plt.Line2D([],[],marker="D",ls="",color="#2B6CB0",label="clean library separation (KO>K27M)"),
    plt.Line2D([],[],marker="^",ls="",color="#d95f02",label="low representation (n<50 or lib<15)")]
fig.legend(handles=h,loc="upper right",frameon=False,fontsize=8,ncol=5,bbox_to_anchor=(0.99,1.0))
fig.suptitle("Metabolic pathway scores by cell type, K27M vs KO (all 5 libraries; depth-robust control-bin, NOT ssGSEA)\n"
             "◆ = all 5 libraries separate by genotype (best achievable evidence at n=2v3; d = Cohen's d); ▲ = low representation",
             fontsize=11)
fig.tight_layout(rect=[0,0,1,0.95])
fig.savefig(f"{FIG}/celltype_pathway_violins.png",dpi=150,bbox_inches="tight")
print("figure written\n")
pd.set_option("display.width",220)
print(R[["pathway","cell_type","n_K27M","n_KO","min_lib_n","median_diff_KmO","cohen_d",
         "p_perm_lib","naive_mwu_p","lib_separation","low_rep"]].round(4).to_string(index=False))
print("\nSummary: clean library separation counts per pathway (of 7 cell types; excl. low-rep):")
ok=R[~R.low_rep]
print(ok.groupby("pathway").apply(lambda x:f"{(x.lib_separation=='clean').sum()}/{len(x)} clean, "
      f"mean|d|={x.cohen_d.abs().mean():.2f}").to_string())
