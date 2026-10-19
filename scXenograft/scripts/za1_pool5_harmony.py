"""Pooled manifold with ALL 5 HSJ019 libraries (KO r2 included). Harmony on sample."""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import numpy as np, pandas as pd, scanpy as sc, anndata as ad, harmonypy
W = REPO_ROOT + "/scXenograft/work"

A = sc.read_h5ad(f"{W}/hsj019_norm.h5ad")   # all 16100 cells = 5 libraries
sc.pp.filter_genes(A, min_cells=3)
print("pooled5 cells:", A.n_obs)
print(A.obs.groupby(["genotype", "sample"], observed=True).size().to_string())

sc.pp.highly_variable_genes(A, n_top_genes=2000, flavor="seurat", batch_key="sample")
H = A[:, A.var.highly_variable].copy()
sc.pp.scale(H, max_value=10)
sc.tl.pca(H, n_comps=50, svd_solver="randomized", random_state=0)
ho = harmonypy.run_harmony(H.obsm["X_pca"][:, :30], A.obs, ["sample"], max_iter_harmony=20)
Z = np.asarray(ho.Z_corr); Z = Z.T if Z.shape[0] != A.n_obs else Z
Z = Z.astype(np.float32)

E = ad.AnnData(np.zeros((A.n_obs, 1), np.float32), obs=A.obs.copy())
E.obsm["X_harmony"] = Z
sc.pp.neighbors(E, n_neighbors=15, use_rep="X_harmony", random_state=0)
sc.tl.umap(E, random_state=0)
sc.tl.leiden(E, resolution=1.0, key_added="leiden", flavor="igraph",
             n_iterations=2, directed=False, random_state=0)
out = E.obs.copy()
out["UMAP1"], out["UMAP2"] = E.obsm["X_umap"][:, 0], E.obsm["X_umap"][:, 1]
out.to_csv(f"{W}/pooled5_embed.csv.gz")
np.savez_compressed(f"{W}/pooled5_harmony.npz", X_harmony=Z,
                    obs_names=np.array(A.obs_names, dtype=object))
# library mixing among KO libs specifically (does r2 integrate now?)
ko = out[out.genotype == "K27M-KO"]
print("\nKO cluster x library (% of cluster), clusters with >=50 KO cells:")
ct = pd.crosstab(ko.leiden, ko["sample"], normalize="index").mul(100)
big = ko.leiden.value_counts(); ct = ct.loc[big[big >= 50].index.sort_values()]
ct.columns = [c.split("_", 2)[-1] for c in ct.columns]
print(ct.round(0).to_string())
