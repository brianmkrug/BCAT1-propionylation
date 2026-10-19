"""HSJ019 split by genotype: independent HVG + PCA within each subset.
Harmony is applied within genotype only, so the batch variable (library) is no
longer nested inside the variable of interest.
Usage: python3 t1_split_pca.py k27m|ko
"""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import sys
import numpy as np, pandas as pd, scanpy as sc

W = REPO_ROOT + "/scXenograft/work"
which = sys.argv[1]
gt = "K27M" if which == "k27m" else "K27M-KO"

A = sc.read_h5ad(f"{W}/hsj019_norm.h5ad")
A = A[A.obs.genotype == gt].copy()
sc.pp.filter_genes(A, min_cells=3)
print(f"{which}: {A.shape}", A.obs["sample"].value_counts().to_dict(), flush=True)

sc.pp.highly_variable_genes(A, n_top_genes=2000, flavor="seurat", batch_key="sample")
A.obs.to_csv(f"{W}/hsj019_{which}_obs.csv.gz")
H = A[:, A.var.highly_variable].copy()
del A
sc.pp.scale(H, max_value=10)
sc.tl.pca(H, n_comps=50, svd_solver="randomized", random_state=0)
np.savez_compressed(f"{W}/hsj019_{which}_pca.npz",
                    X_pca=H.obsm["X_pca"].astype(np.float32),
                    variance_ratio=H.uns["pca"]["variance_ratio"])
vr = H.uns["pca"]["variance_ratio"]
print(f"{which}: PCA {H.shape}; PC1-10 var {vr[:10].round(3)}; cum50 {vr.sum():.3f}")
