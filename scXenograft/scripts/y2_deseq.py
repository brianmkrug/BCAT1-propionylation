"""PyDESeq2 on the library-level pseudobulk, restricted to the gene symbols
shared between the hg19 single-cell annotation and the hg38 bulk annotation.

Matches the bulk pipeline: design ~condition, contrast K27M vs KO,
prefilter total count >= 10. Positive log2FC = higher in K27M.
"""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import numpy as np, pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

B = REPO_ROOT + ""
OUT = f"{B}/scXenograft"

sc_pb = pd.read_csv(f"{OUT}/sc_pseudobulk_counts_by_library.csv", index_col=0)
bulk = pd.read_csv(f"{B}/xHSJ019_K27M_vs_KO_analysis/gene_counts_raw.csv", index_col=0)
bulk_de = pd.read_csv(f"{B}/xHSJ019_K27M_vs_KO_analysis/deseq2_K27M_vs_KO.csv", index_col=0)

# ---- gene universe intersection ----
sc_g, bk_g = set(sc_pb.index), set(bulk.index)
shared = sorted(sc_g & bk_g)
print("Gene universe")
print(f"  hg19 single-cell (post gene-filter): {len(sc_g):>6}")
print(f"  hg38 bulk (all quantified):          {len(bk_g):>6}")
print(f"  shared symbols:                      {len(shared):>6}")
print(f"  sc-only: {len(sc_g-bk_g)}   bulk-only: {len(bk_g-sc_g)}")
tested = set(bulk_de.index)
print(f"  shared AND tested in the bulk DE table: {len(set(shared) & tested)}")
pd.Series(sorted(sc_g & bk_g)).to_csv(f"{OUT}/shared_gene_universe_hg19_hg38.csv",
                                      index=False, header=["gene"])

X = sc_pb.loc[shared]
X = X.loc[X.sum(1) >= 10]                      # same prefilter as the bulk run
print(f"\ngenes tested after prefilter (total count >= 10): {X.shape[0]}")

meta = pd.DataFrame(
    {"condition": ["K27M" if c.startswith("K27M") else "KO" for c in X.columns]},
    index=X.columns)
print(meta.to_string())

dds = DeseqDataSet(counts=X.T, metadata=meta, design="~condition", refit_cooks=True)
dds.deseq2()
st = DeseqStats(dds, contrast=["condition", "K27M", "KO"])
st.summary()
res = st.results_df.sort_values("pvalue")
res.index.name = "gene"
res.to_csv(f"{OUT}/sc_pseudobulk_deseq2_K27M_vs_KO.csv")

sfkey = "size_factors" if "size_factors" in dds.obsm else ("size_factors" if "size_factors" in dds.obs else None)
sf = pd.Series(dds.obs["size_factors"].values if sfkey else np.nan, index=X.columns)
print("\nsize factors:"); print(sf.round(3).to_string())
sig = res.padj < 0.05
print(f"\npadj<0.05: {int(sig.sum())}")
print(f"padj<0.05 & |log2FC|>1: up in K27M {int((sig & (res.log2FoldChange>1)).sum())}, "
      f"up in KO {int((sig & (res.log2FoldChange<-1)).sum())}")
print("\nTop 15 by p-value:")
print(res.head(15)[["baseMean", "log2FoldChange", "stat", "padj"]].round(3).to_string())
print("\nBCAT1:")
print(res.loc[["BCAT1"]][["baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]].round(4).to_string()
      if "BCAT1" in res.index else "  not tested")
