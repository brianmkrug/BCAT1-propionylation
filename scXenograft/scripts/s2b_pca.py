"""Stage 2b: scale HVG subset, PCA. Saves compact companion .npz (not the big object).
Usage: python3 s2b_pca.py bt245|hsj019"""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import sys
import numpy as np, scanpy as sc
W = REPO_ROOT + "/scXenograft/work"
tag = sys.argv[1]
A = sc.read_h5ad(f"{W}/{tag}_norm.h5ad")
hv = A.var.highly_variable.values
H = A[:, hv].copy()
del A
sc.pp.scale(H, max_value=10)
sc.tl.pca(H, n_comps=50, svd_solver="randomized", random_state=0)
np.savez_compressed(f"{W}/{tag}_pca.npz",
                    X_pca=H.obsm["X_pca"].astype(np.float32),
                    obs_names=np.array(H.obs_names, dtype=object),
                    variance_ratio=H.uns["pca"]["variance_ratio"])
vr = H.uns["pca"]["variance_ratio"]
print(f"{tag}: {H.shape} PCA done; PC1-10 var = {vr[:10].round(3)}; cum50={vr.sum():.3f}")
