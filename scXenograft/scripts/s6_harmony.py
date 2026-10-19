"""Stage 6: Harmony correction on `sample` (NOT genotype), then re-embed.
Usage: python3 s6_harmony.py bt245|hsj019
"""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import sys
import numpy as np, pandas as pd, scanpy as sc, anndata as ad, harmonypy

W = REPO_ROOT + "/scXenograft/work"
tag = sys.argv[1]
z = np.load(f"{W}/{tag}_pca.npz", allow_pickle=True)
obs = pd.read_csv(f"{W}/{tag}_obs.csv.gz", index_col=0)
pcs = z["X_pca"][:, :30]

ho = harmonypy.run_harmony(pcs, obs, ["sample"], max_iter_harmony=20)
H = np.asarray(ho.Z_corr)
if H.shape[0] != len(obs):
    H = H.T
H = H.astype(np.float32)
assert H.shape == (len(obs), pcs.shape[1]), H.shape
np.savez_compressed(f"{W}/{tag}_harmony.npz", X_harmony=H)

A = ad.AnnData(np.zeros((len(obs), 1), dtype=np.float32), obs=obs)
A.obsm["X_harmony"] = H
sc.pp.neighbors(A, n_neighbors=15, use_rep="X_harmony", random_state=0)
sc.tl.umap(A, random_state=0)
sc.tl.leiden(A, resolution=0.8, key_added="leiden", flavor="igraph",
             n_iterations=2, directed=False, random_state=0)
out = A.obs.copy()
out["UMAP1"], out["UMAP2"] = A.obsm["X_umap"][:, 0], A.obsm["X_umap"][:, 1]
out.to_csv(f"{W}/{tag}_embed_harmony.csv.gz")
print(tag, "clusters:", out.leiden.nunique())
print("\ncluster x sample (% of cluster):")
print((pd.crosstab(out.leiden, out["sample"], normalize="index")*100).round(1).to_string())
