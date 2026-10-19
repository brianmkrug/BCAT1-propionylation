"""Prerank GSEA on the sc pseudobulk Wald statistic, matching the bulk run
(Hallmark + curated KEGG metabolism, 1000 perms, size 5-500, seed 42)."""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import numpy as np, pandas as pd, gseapy as gp

B = REPO_ROOT + ""
OUT = f"{B}/scXenograft"
res = pd.read_csv(f"{OUT}/sc_pseudobulk_deseq2_K27M_vs_KO.csv", index_col=0)
rnk = res.stat.dropna().sort_values(ascending=False)
rnk.to_csv(f"{OUT}/sc_pseudobulk_ranked_genes.rnk", sep="\t", header=False)
print(f"ranked genes: {len(rnk)}")

for tag, gmt in [("hallmark", f"{B}/signatures/h.all.v7.5.symbols.gmt"),
                 ("kegg_metab", f"{B}/xHSJ019_K27M_vs_KO_analysis/KEGG_metabolism_expanded.gmt")]:
    pre = gp.prerank(rnk=rnk.to_frame(), gene_sets=gmt, permutation_num=1000,
                     min_size=5, max_size=500, seed=42, threads=4,
                     outdir=None, verbose=False)
    r = pre.res2d.copy()
    for c in ["NES", "FDR q-val", "NOM p-val"]:
        r[c] = pd.to_numeric(r[c], errors="coerce")
    r = r.sort_values("NES", ascending=False)
    r.to_csv(f"{OUT}/sc_pseudobulk_gsea_{tag}.csv", index=False)
    sig = r[r["FDR q-val"] < 0.25]
    print(f"\n=== {tag}: {len(sig)} sets at FDR<0.25 of {len(r)} ===")
    cols = ["Term", "NES", "NOM p-val", "FDR q-val"]
    print("-- top 12 enriched in K27M --")
    print(sig.head(12)[cols].to_string(index=False))
    print("-- top 6 enriched in KO --")
    print(sig.tail(6)[cols].to_string(index=False))

# direct comparison to the bulk GSEA
bh = pd.read_csv(f"{B}/xHSJ019_K27M_vs_KO_analysis/gsea_hallmark_results.csv")
sh = pd.read_csv(f"{OUT}/sc_pseudobulk_gsea_hallmark.csv")
tcol = "Term" if "Term" in bh.columns else bh.columns[0]
m = sh.set_index("Term").NES.to_frame("NES_sc").join(
    bh.set_index(tcol).NES.to_frame("NES_bulk"), how="inner")
from scipy.stats import spearmanr, pearsonr
print(f"\n=== Hallmark NES, sc pseudobulk vs bulk ({len(m)} sets) ===")
print(f"  Pearson r {pearsonr(m.NES_sc, m.NES_bulk)[0]:+.3f}, "
      f"Spearman {spearmanr(m.NES_sc, m.NES_bulk).statistic:+.3f}, "
      f"same sign {100*(np.sign(m.NES_sc)==np.sign(m.NES_bulk)).mean():.0f}%")
key = ["HALLMARK_E2F_TARGETS","HALLMARK_G2M_CHECKPOINT","HALLMARK_MYC_TARGETS_V1",
       "HALLMARK_MTORC1_SIGNALING","HALLMARK_CHOLESTEROL_HOMEOSTASIS",
       "HALLMARK_GLYCOLYSIS","HALLMARK_FATTY_ACID_METABOLISM","HALLMARK_HYPOXIA",
       "HALLMARK_TNFA_SIGNALING_VIA_NFKB","HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]
print(m.reindex([k for k in key if k in m.index]).round(2).to_string())
m.round(3).to_csv(f"{OUT}/sc_vs_bulk_hallmark_NES.csv")
