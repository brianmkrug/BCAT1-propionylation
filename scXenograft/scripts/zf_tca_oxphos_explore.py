"""Per-gene stats for TCA cycle and OxPhos on pooled cells, classified by
step/complex. Same metric as the cholesterol analysis."""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import re
import numpy as np, pandas as pd, scanpy as sc

B = REPO_ROOT + ""
W = f"{B}/scXenograft/work"; OUT = f"{B}/scXenograft"
NB = 20

def read_set(name):
    for ln in open(f"{B}/xHSJ019_K27M_vs_KO_analysis/KEGG_metabolism_expanded.gmt"):
        f = ln.rstrip("\n").split("\t")
        if f[0] == name: return [x for x in f[2:] if x]
    return []
TCA = read_set("Citrate cycle (TCA cycle)")
OXP = read_set("Oxidative phosphorylation")
print(f"TCA set: {len(TCA)} genes | OxPhos set: {len(OXP)} genes")

def oxphos_complex(g):
    if g.startswith("NDUF") or g in ("MT-ND1","MT-ND2","MT-ND3","MT-ND4","MT-ND4L",
                                      "MT-ND5","MT-ND6"): return "CI (NADH dh)"
    if g.startswith("SDH"): return "CII (SDH)"
    if g.startswith("UQCR") or g in ("CYC1","MT-CYB"): return "CIII (bc1)"
    if g.startswith("COX") or g.startswith("COA") or g.startswith("MT-CO") or g=="NDUFA4":
        return "CIV (COX)"
    if g.startswith("ATP5") or g in ("MT-ATP6","MT-ATP8"): return "CV (ATP synthase)"
    if g.startswith("ATP6V") or g.startswith("ATP6AP") or g=="TCIRG1": return "V-ATPase (lysosomal)"
    return "other"

TCA_STEP = {
 "PDH feed-in": ["PDHA1","PDHB","DLAT","DLD","PDHX","PC"],
 "Citrate synthase": ["CS"],
 "Aconitase": ["ACO1","ACO2"],
 "Isocitrate dh": ["IDH1","IDH2","IDH3A","IDH3B","IDH3G"],
 "a-KG dh": ["OGDH","OGDHL","DLST"],
 "Succinyl-CoA synth": ["SUCLA2","SUCLG1","SUCLG2"],
 "Succinate dh (CII)": ["SDHA","SDHB","SDHC","SDHD"],
 "Fumarase": ["FH"],
 "Malate dh": ["MDH1","MDH2","MDH1B"],
 "Anaplerosis/cataplerosis": ["PCK1","PCK2","ACLY"],
}
def tca_step(g):
    for k, v in TCA_STEP.items():
        if g in v: return k
    return "other"

emb = pd.read_csv(f"{W}/pooled_dpt.csv.gz", index_col=0)
emb["bin"] = pd.qcut(emb.dpt, NB, labels=False, duplicates="drop")
A = sc.read_h5ad(f"{W}/hsj019_norm.h5ad")[emb.index]
de = pd.read_csv(f"{OUT}/sc_pseudobulk_deseq2_K27M_vs_KO.csv", index_col=0)
isK = (emb.genotype == "K27M").values

def stats(genes, kind):
    genes = [g for g in genes if g in A.var_names]
    E = pd.DataFrame(np.asarray(A[:, genes].X.todense()), index=emb.index, columns=genes)
    rows = []
    for g in genes:
        v = E[g].values
        kb = np.array([v[isK & (emb.bin==b).values].mean() for b in range(NB)])
        ob = np.array([v[~isK & (emb.bin==b).values].mean() for b in range(NB)])
        pooled = np.r_[kb, ob]; mu, sd = pooled.mean(), pooled.std()
        binadj = np.nanmean((kb-mu)/sd - (ob-mu)/sd) if sd>0 else np.nan
        grp = oxphos_complex(g) if kind=="oxp" else tca_step(g)
        rows.append(dict(gene=g, group=grp, det_pct=100*(v>0).mean(),
                         mean_expr=v.mean(), binadj_z_KmO=binadj,
                         pb_log2FC=de.log2FoldChange.get(g,np.nan),
                         pb_padj=de.padj.get(g,np.nan)))
    return pd.DataFrame(rows).set_index("gene")

st_t = stats(TCA, "tca"); st_t.to_csv(f"{OUT}/tca_gene_stats.csv")
st_o = stats(OXP, "oxp"); st_o.to_csv(f"{OUT}/oxphos_gene_stats.csv")
pd.set_option("display.width", 200)
print("\n=== TCA by step (det>=12) ===")
print(st_t[st_t.det_pct>=12].sort_values(["group","binadj_z_KmO"])
      [["group","det_pct","binadj_z_KmO","pb_log2FC","pb_padj"]].round(3).to_string())
print("\n=== OxPhos: mean genotype lean per complex (det>=12) ===")
o = st_o[st_o.det_pct>=12]
print(o.groupby("group").agg(n=("binadj_z_KmO","size"),
      mean_z_KmO=("binadj_z_KmO","mean"), mean_log2FC=("pb_log2FC","mean")).round(3).to_string())
print("\nOxPhos genes count by complex (det>=12):", o.group.value_counts().to_dict())
print("\nMT-encoded genes present:", [g for g in OXP if g.startswith("MT-") and g in A.var_names])
