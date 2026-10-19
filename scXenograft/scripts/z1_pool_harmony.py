"""Pooled K27M + KO(r1+r3) manifold for a COMMON pseudotime axis.

Harmony is applied on `sample` (4 libraries). Genotype is nested in library, so
this also removes genotype from the *embedding* -- which is intended: the
embedding only defines a genotype-neutral lineage coordinate. The genotype
comparison is done later in pathway-score space, stratified by pseudotime bin.
"""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import numpy as np, pandas as pd, scanpy as sc, anndata as ad, harmonypy

W = REPO_ROOT + "/scXenograft/work"

A = sc.read_h5ad(f"{W}/hsj019_norm.h5ad")
keep = ((A.obs.genotype == "K27M") |
        ((A.obs.genotype == "K27M-KO") & (A.obs.replicate != "r2")))
A = A[keep].copy()
sc.pp.filter_genes(A, min_cells=3)
print("pooled cells:", A.n_obs)
print(A.obs.groupby(["genotype", "sample"], observed=True).size().to_string())

sc.pp.highly_variable_genes(A, n_top_genes=2000, flavor="seurat", batch_key="sample")
hv = A.var.highly_variable.values
H = A[:, hv].copy()
sc.pp.scale(H, max_value=10)
sc.tl.pca(H, n_comps=50, svd_solver="randomized", random_state=0)
pcs = H.obsm["X_pca"][:, :30]

ho = harmonypy.run_harmony(pcs, A.obs, ["sample"], max_iter_harmony=20)
Z = np.asarray(ho.Z_corr)
if Z.shape[0] != A.n_obs: Z = Z.T
Z = Z.astype(np.float32)

E = ad.AnnData(np.zeros((A.n_obs, 1), np.float32), obs=A.obs.copy())
E.obsm["X_harmony"] = Z
sc.pp.neighbors(E, n_neighbors=15, use_rep="X_harmony", random_state=0)
sc.tl.umap(E, random_state=0)
sc.tl.leiden(E, resolution=1.0, key_added="leiden", flavor="igraph",
             n_iterations=2, directed=False, random_state=0)

out = E.obs.copy()
out["UMAP1"], out["UMAP2"] = E.obsm["X_umap"][:, 0], E.obsm["X_umap"][:, 1]
out.to_csv(f"{W}/pooled_embed.csv.gz")
np.savez_compressed(f"{W}/pooled_harmony.npz", X_harmony=Z,
                    obs_names=np.array(A.obs_names, dtype=object))

# mixing diagnostic (genotype, not library): expected genotype fraction
exp = (out.genotype == "K27M").mean()
ct = out.groupby("leiden", observed=True).apply(
    lambda d: (d.genotype == "K27M").mean())
print(f"\noverall K27M fraction {exp:.2f}")
print("per-cluster K27M fraction (0=pure KO, 1=pure K27M):")
print(ct.round(2).to_string())
print(f"\nclusters that are >90% one genotype: "
      f"{int(((ct>0.9)|(ct<0.1)).sum())}/{len(ct)}")
