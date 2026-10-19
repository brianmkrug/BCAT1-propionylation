"""BCAA gene matrix on ALL 5 libraries (pooled5), then (1) clustered heatmap and
(2) BCAT1/BCAT2-only heatmap — the 5-library versions."""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import numpy as np, pandas as pd, scanpy as sc
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.cluster.hierarchy import linkage, leaves_list

B = REPO_ROOT + ""; W = f"{B}/scXenograft/work"
OUT = f"{B}/scXenograft"; FIG = f"{OUT}/figures"
NB = 20; DET = 12
ALIAS = {"MMUT": "MUT"}
GROUPS = [  # same curated union as the 4-lib figure
 ("Transamination", ["BCAT1","BCAT2"]),
 ("BCKDH complex", ["BCKDHA","BCKDHB","DBT","DLD"]),
 ("Acyl-CoA dehydrogenases", ["IVD","ACAD8","ACADSB","ACADM","ACADS"]),
 ("Leucine ketogenic", ["MCCC1","MCCC2","AUH","HMGCL","HMGCLL1","HMGCS1","OXCT1","OXCT2","ACAT1","ACAT2"]),
 ("Val/Ile", ["ECHS1","HIBCH","HIBADH","ALDH6A1","ABAT","HADH","HADHA","HADHB","HSD17B10",
              "ACAA1","ACAA2","EHHADH","AACS","ALDH2","ALDH3A2","ALDH7A1","ALDH9A1","ALDH1B1","AOX1","IL4I1","ACSF3","AGXT2"]),
 ("Propionyl-CoA node", ["PCCA","PCCB","MCEE","MMUT"]),
 ("Propanoate/other", ["SUCLG1","SUCLG2","SUCLA2","ACSS1","ACSS2","ACSS3","ACACA","ACACB","MLYCD","ACOX1","ACOX3","ECHDC1","LDHA","LDHB","LDHC"]),
]
emb = pd.read_csv(f"{W}/pooled5_dpt.csv.gz", index_col=0)
emb["bin"] = pd.qcut(emb.dpt, NB, labels=False, duplicates="drop")
A = sc.read_h5ad(f"{W}/hsj019_norm.h5ad")[emb.index]
de = pd.read_csv(f"{OUT}/sc_pseudobulk_deseq2_K27M_vs_KO.csv", index_col=0)
vset = set(A.var_names)
def ref(g): return ALIAS.get(g, g)
def det(g):
    r = ref(g); return 100*(np.asarray(A[:, r].X.todense()).ravel() > 0).mean() if r in vset else 0
ordered, labels, seen = [], [], set()
for _, gg in GROUPS:
    for g in gg:
        r = ref(g)
        if r in vset and det(g) >= DET and r not in seen:
            ordered.append(r); seen.add(r)
            labels.append("MMUT/MUT" if g == "MMUT" else g)
E = pd.DataFrame(np.asarray(A[:, ordered].X.todense()), index=emb.index, columns=ordered)
def zmat(gt):
    return np.vstack([[E[g].values[(emb.genotype==gt).values & (emb.bin==b).values].mean()
                       for b in range(NB)] for g in ordered])
Kr, Or = zmat("K27M"), zmat("K27M-KO")
mu = np.array([np.r_[Kr[i],Or[i]].mean() for i in range(len(ordered))])[:,None]
sd = np.array([np.r_[Kr[i],Or[i]].std() for i in range(len(ordered))])[:,None]
K, O = (Kr-mu)/sd, (Or-mu)/sd
lfc = np.array([de.log2FoldChange.get(g, np.nan) for g in ordered])
vals = pd.DataFrame(K, index=ordered, columns=[f"K27M_bin{b}" for b in range(NB)]).join(
       pd.DataFrame(O, index=ordered, columns=[f"KO_bin{b}" for b in range(NB)]))
vals["pb_log2FC"] = lfc
vals.round(3).to_csv(f"{OUT}/bcaa_gene_heatmap_values_5lib.csv")

def paired(Ksub, Osub, gl, lsub, lfcsub, title, fname, h):
    VZ = 2.5
    fig = plt.figure(figsize=(14.5, h))
    gs = GridSpec(1, 4, width_ratios=[1,1,0.06,0.045], wspace=0.05)
    axK, axO = fig.add_subplot(gs[0,0]), fig.add_subplot(gs[0,1])
    for ax, Mx, ttl in [(axK,Ksub,"K27M"),(axO,Osub,"K27M-KO (r1+r2+r3)")]:
        im = ax.imshow(Mx, aspect="auto", cmap="RdBu_r", vmin=-VZ, vmax=VZ)
        ax.set_xticks([0,NB//2,NB-1]); ax.set_xticklabels(["early","mid","late"]); ax.set_xlabel("pooled pseudotime bin")
        ax.set_title(ttl, fontsize=11)
    axK.set_yticks(range(len(gl))); axK.set_yticklabels(lsub, fontsize=7 if len(gl)>5 else 11); axO.set_yticks([])
    axL = fig.add_subplot(gs[0,2]); axL.imshow(lfcsub[:,None], aspect="auto", cmap="PuOr_r", vmin=-1.5, vmax=1.5)
    axL.set_xticks([]); axL.set_yticks([]); axL.set_title("log2FC\nK27M/KO", fontsize=6.5)
    fig.colorbar(im, cax=fig.add_subplot(gs[0,3]), label="pooled z (per gene)")
    fig.suptitle(title, fontsize=11, y=1.02 if len(gl)>5 else 1.18)
    fig.subplots_adjust(left=0.14 if len(gl)>5 else 0.10, top=0.94 if len(gl)>5 else 0.78)
    fig.savefig(f"{FIG}/{fname}", dpi=170, bbox_inches="tight"); plt.close(fig)

# (1) clustered, all genes mixed
feat = np.nan_to_num(np.hstack([K, O]))
order = leaves_list(linkage(feat, method="average", metric="correlation"))
gl = [labels[i] for i in order]
paired(K[order], O[order], gl, gl, lfc[order],
       "v2 (all 5 libraries) — BCAA metabolism genes clustered by similarity, K27M vs KO\n"
       "pooled z per gene; rows = hierarchical clustering (correlation distance); strip = pseudobulk log2FC",
       "bcaa_gene_heatmap_clustered_5lib.png", max(11, 0.26*len(gl)))

# (2) BCAT1 / BCAT2 only
bi = [ordered.index("BCAT1"), ordered.index("BCAT2")]
paired(K[bi], O[bi], ["BCAT1","BCAT2"], ["BCAT1","BCAT2"], lfc[bi],
       "BCAT1 / BCAT2 across pooled pseudotime, K27M vs KO (all 5 libraries; pooled z per gene; strip = pseudobulk log2FC)",
       "bcaa_gene_heatmap_bcat_only_5lib.png", 2.2)
print(f"{len(ordered)} genes; BCAT1 log2FC {lfc[bi[0]]:+.2f}, BCAT2 {lfc[bi[1]]:+.2f}")
print("clustered order:", ", ".join(gl))
