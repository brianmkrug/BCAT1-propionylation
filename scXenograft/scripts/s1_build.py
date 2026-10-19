"""Stage 1: per-sample human-cell extraction + methods QC -> one h5ad per sample.
Usage: python3 s1_build.py <sample_prefix> [...]
Keeps hg19 features only (prefix stripped), raw integer counts in .X.
"""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import gzip, os, sys
import numpy as np, pandas as pd, scipy.io as sio, scipy.sparse as sp
import anndata as ad

D = REPO_ROOT + "/scXenograft/datasets"
W = REPO_ROOT + "/scXenograft/work"
os.makedirs(W, exist_ok=True)

def build(s):
    feat = pd.read_csv(f"{D}/{s}.features.tsv.gz", sep="\t", header=None,
                       names=["id", "symbol", "type"])
    bc = pd.read_csv(f"{D}/{s}.barcodes.tsv.gz", sep="\t", header=None)[0].values
    with gzip.open(f"{D}/{s}.matrix.mtx.gz") as fh:
        M = sio.mmread(fh).tocsr()                      # genes x cells
    hg = feat.id.str.startswith("hg19_").values
    hg_umi = np.asarray(M[hg].sum(0)).ravel()
    mm_umi = np.asarray(M[~hg].sum(0)).ravel()
    tot = np.maximum(hg_umi + mm_umi, 1)
    frac = hg_umi / tot

    A = ad.AnnData(sp.csr_matrix(M[hg].T.astype(np.float32)))   # cells x hg19 genes
    A.obs_names = [f"{s}_{b}" for b in bc]
    hf = feat[hg].copy()
    A.var_names = pd.Index(hf.symbol.str.replace("hg19_", "", regex=False)).astype(str)
    A.var["ensembl"] = hf.id.str.replace("hg19_", "", regex=False).values
    A.var_names_make_unique()

    A.obs["sample"] = s
    A.obs["gsm"] = s.split("_")[0]
    A.obs["line"] = "BT245" if "_bt245_" in s else "HSJ019"
    A.obs["genotype"] = "K27M" if "_k27m_" in s else "K27M-KO"
    A.obs["replicate"] = s.rsplit("_", 1)[-1]
    A.obs["hg19_umi"] = hg_umi
    A.obs["mm10_umi"] = mm_umi
    A.obs["frac_hg19"] = frac
    A.obs["mm10_frac"] = 1 - frac

    mito = A.var_names.str.startswith("MT-")
    A.obs["n_genes"] = np.asarray((A.X > 0).sum(1)).ravel()
    A.obs["n_umi"] = np.asarray(A.X.sum(1)).ravel()
    A.obs["pct_mt"] = 100 * np.asarray(A[:, mito].X.sum(1)).ravel() / np.maximum(A.obs.n_umi, 1)

    n0 = A.n_obs
    A = A[A.obs.frac_hg19 > 0.75].copy()                # human cells
    n_hum = A.n_obs
    A = A[(A.obs.pct_mt <= 50) & (A.obs.n_genes >= 500) & (A.obs.n_umi >= 500)].copy()
    for col in ["n_genes", "n_umi"]:                    # +/- 2 SD, per sample, sequential
        m, sd = A.obs[col].mean(), A.obs[col].std()
        A = A[(A.obs[col] > m - 2*sd) & (A.obs[col] < m + 2*sd)].copy()

    A.write_h5ad(f"{W}/{s}.h5ad")
    print(f"{s}: {n0} barcodes -> {n_hum} human -> {A.n_obs} post-QC | "
          f"{A.n_vars} hg19 genes", flush=True)

for s in sys.argv[1:]:
    build(s)
