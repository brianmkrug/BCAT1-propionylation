"""Harmony on `sample` WITHIN a genotype subset, then UMAP + Leiden.
Here library is a pure technical variable (not nested in genotype), so this
correction is legitimate. Usage: python3 t2_harmony.py k27m|ko
"""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import sys
import numpy as np, pandas as pd, scanpy as sc, anndata as ad, harmonypy

W = REPO_ROOT + "/scXenograft/work"
which = sys.argv[1]
obs = pd.read_csv(f"{W}/hsj019_{which}_obs.csv.gz", index_col=0)
pcs = np.load(f"{W}/hsj019_{which}_pca.npz")["X_pca"][:, :30]

ho = harmonypy.run_harmony(pcs, obs, ["sample"], max_iter_harmony=20)
H = np.asarray(ho.Z_corr)
if H.shape[0] != len(obs):
    H = H.T
H = H.astype(np.float32)
np.savez_compressed(f"{W}/hsj019_{which}_harmony.npz", X_harmony=H)

A = ad.AnnData(np.zeros((len(obs), 1), dtype=np.float32), obs=obs)
A.obsm["X_pca"], A.obsm["X_harmony"] = pcs, H
res = {}
for rep, key in [("X_pca", "uncorr"), ("X_harmony", "harmony")]:
    sc.pp.neighbors(A, n_neighbors=15, use_rep=rep, random_state=0)
    sc.tl.umap(A, random_state=0)
    sc.tl.leiden(A, resolution=0.8, key_added=f"leiden_{key}", flavor="igraph",
                 n_iterations=2, directed=False, random_state=0)
    res[f"UMAP1_{key}"] = A.obsm["X_umap"][:, 0]
    res[f"UMAP2_{key}"] = A.obsm["X_umap"][:, 1]

out = A.obs.copy()
for k, v in res.items():
    out[k] = v
out.to_csv(f"{W}/hsj019_{which}_embed.csv.gz")
print(f"{which}: {out.shape[0]} cells, "
      f"{out.leiden_harmony.nunique()} harmony clusters "
      f"({out.leiden_uncorr.nunique()} uncorrected)")
for key in ["uncorr", "harmony"]:
    ct = pd.crosstab(out[f"leiden_{key}"], out["sample"], normalize="index").mul(100)
    ct.columns = [c.split("_", 2)[-1] for c in ct.columns]
    print(f"\n--- {key}: cluster x library (% of cluster) ---")
    print(ct.round(1).to_string())
