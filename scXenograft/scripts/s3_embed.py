"""Stage 3: neighbors + UMAP + Leiden on the PCA space. Operates on a lightweight
AnnData (obs + X_pca only) so it never touches the big expression matrix.
Usage: python3 s3_embed.py bt245|hsj019
"""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import sys
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
W = REPO_ROOT + "/scXenograft/work"
tag = sys.argv[1]
z = np.load(f"{W}/{tag}_pca.npz", allow_pickle=True)
obs = pd.read_csv(f"{W}/{tag}_obs.csv.gz", index_col=0)
A = ad.AnnData(np.zeros((len(obs), 1), dtype=np.float32), obs=obs)
A.obsm["X_pca"] = z["X_pca"]
sc.pp.neighbors(A, n_neighbors=15, n_pcs=30, use_rep="X_pca", random_state=0)
sc.tl.umap(A, random_state=0)
sc.tl.leiden(A, resolution=1.0, key_added="leiden", flavor="igraph",
             n_iterations=2, directed=False, random_state=0)
out = A.obs.copy()
out["UMAP1"] = A.obsm["X_umap"][:, 0]
out["UMAP2"] = A.obsm["X_umap"][:, 1]
out.to_csv(f"{W}/{tag}_embed.csv.gz")
print(tag, "clusters:", out.leiden.nunique())
ct = pd.crosstab(out.leiden, out["sample"])
print(ct.to_string())
print("\nper-cluster genotype fraction:")
print((pd.crosstab(out.leiden, out.genotype, normalize="index")*100).round(1).to_string())
