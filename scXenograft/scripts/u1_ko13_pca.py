"""HSJ019 K27M-KO, r2 excluded (user call: r2 is an outlier w.r.t. expected cell
states). HVG + PCA on r1 + r3 only.
"""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import numpy as np, pandas as pd, scanpy as sc
W = REPO_ROOT + "/scXenograft/work"

A = sc.read_h5ad(f"{W}/hsj019_norm.h5ad")
A = A[(A.obs.genotype == "K27M-KO") & (A.obs.replicate != "r2")].copy()
sc.pp.filter_genes(A, min_cells=3)
print("ko_r1r3:", A.shape, A.obs["sample"].value_counts().to_dict(), flush=True)

sc.pp.highly_variable_genes(A, n_top_genes=2000, flavor="seurat", batch_key="sample")
A.obs.to_csv(f"{W}/hsj019_ko13_obs.csv.gz")
H = A[:, A.var.highly_variable].copy()
del A
sc.pp.scale(H, max_value=10)
sc.tl.pca(H, n_comps=50, svd_solver="randomized", random_state=0)
np.savez_compressed(f"{W}/hsj019_ko13_pca.npz",
                    X_pca=H.obsm["X_pca"].astype(np.float32),
                    variance_ratio=H.uns["pca"]["variance_ratio"])
vr = H.uns["pca"]["variance_ratio"]
print(f"ko13: PCA {H.shape}; PC1-10 var {vr[:10].round(3)}; cum50 {vr.sum():.3f}")
