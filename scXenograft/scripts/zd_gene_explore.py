"""Per-gene stats for the two cholesterol pathways on the pooled cells:
detection rate, mean expression, matched-pseudotime genotype difference (z),
and sc-pseudobulk log2FC. Used to pick genes and identify leading edge."""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import numpy as np, pandas as pd, scipy.sparse as sp, scanpy as sc

B = REPO_ROOT + ""
W = f"{B}/scXenograft/work"; OUT = f"{B}/scXenograft"
NB = 20

def read_set(name):
    for ln in open(f"{W}/metabolic_panel.gmt"):
        f = ln.rstrip("\n").split("\t")
        if f[0] == name: return [x for x in f[2:] if x]
BIOSYN = read_set("CHOL_Reactome_cholesterol_biosynthesis")
UPTAKE = read_set("CHOL_KEGGc_cholesterol_metabolism")
print(f"Reactome cholesterol biosynthesis: {len(BIOSYN)} genes")
print(f"KEGG cholesterol metabolism: {len(UPTAKE)} genes")

emb = pd.read_csv(f"{W}/pooled_dpt.csv.gz", index_col=0)
emb["bin"] = pd.qcut(emb.dpt, NB, labels=False, duplicates="drop")
A = sc.read_h5ad(f"{W}/hsj019_norm.h5ad")[emb.index]
de = pd.read_csv(f"{OUT}/sc_pseudobulk_deseq2_K27M_vs_KO.csv", index_col=0)

allg = [g for g in dict.fromkeys(BIOSYN + UPTAKE) if g in A.var_names]
E = pd.DataFrame(np.asarray(A[:, allg].X.todense()), index=emb.index, columns=allg)
isK = (emb.genotype == "K27M").values

rows = []
for g in allg:
    v = E[g].values
    det = 100*(v > 0).mean()
    # matched-pseudotime genotype diff on pooled z of bin-means
    kb = np.array([v[isK & (emb.bin == b).values].mean() for b in range(NB)])
    ob = np.array([v[~isK & (emb.bin == b).values].mean() for b in range(NB)])
    pooled = np.r_[kb, ob]; mu, sd = pooled.mean(), pooled.std()
    zk, zo = (kb-mu)/sd, (ob-mu)/sd
    binadj_z = np.nanmean(zk - zo)
    rows.append(dict(gene=g, in_biosyn=g in BIOSYN, in_uptake=g in UPTAKE,
                     det_pct=det, mean_expr=v.mean(),
                     binadj_z_KmO=binadj_z,
                     pb_log2FC=de.log2FoldChange.get(g, np.nan),
                     pb_padj=de.padj.get(g, np.nan)))
R = pd.DataFrame(rows).set_index("gene")
R.to_csv(f"{OUT}/cholesterol_gene_stats.csv")
pd.set_option("display.width", 200)

print("\n=== BIOSYNTHESIS genes (detection-sorted) ===")
b = R[R.in_biosyn].sort_values("det_pct", ascending=False)
print(b[["det_pct","mean_expr","binadj_z_KmO","pb_log2FC","pb_padj"]].round(3).to_string())
print("\n=== UPTAKE/TRANSPORT genes (detection-sorted) ===")
u = R[R.in_uptake].sort_values("det_pct", ascending=False)
print(u[["det_pct","mean_expr","binadj_z_KmO","pb_log2FC","pb_padj"]].round(3).to_string())

print("\nLeading edge for K27M>KO (biosynthesis, binadj_z>0, det>=10):")
le = R[(R.det_pct>=10)].sort_values("binadj_z_KmO", ascending=False)
print(le.head(15)[["det_pct","binadj_z_KmO","pb_log2FC","in_biosyn","in_uptake"]].round(3).to_string())
print("\nLeading edge for KO>K27M (binadj_z<0, det>=10):")
print(le.tail(15)[["det_pct","binadj_z_KmO","pb_log2FC","in_biosyn","in_uptake"]].round(3).to_string())
