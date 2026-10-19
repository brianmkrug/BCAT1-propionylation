"""Stage 7: assemble the final merged .h5ad per cell line with all annotations.
Usage: python3 s7_finalize.py bt245|hsj019
"""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import sys
import numpy as np, pandas as pd, scanpy as sc

W = REPO_ROOT + "/scXenograft/work"
OUT = REPO_ROOT + "/scXenograft"
tag = sys.argv[1]

A = sc.read_h5ad(f"{W}/{tag}_norm.h5ad")
ann = pd.read_csv(f"{W}/{tag}_annotated.csv.gz", index_col=0).loc[A.obs_names]
sig = pd.read_csv(f"{W}/{tag}_sigscores.csv.gz", index_col=0).loc[A.obs_names]
har = pd.read_csv(f"{W}/{tag}_embed_harmony.csv.gz", index_col=0).loc[A.obs_names]
pca = np.load(f"{W}/{tag}_pca.npz", allow_pickle=True)["X_pca"]
hz = np.load(f"{W}/{tag}_harmony.npz")["X_harmony"]

for c in ["leiden", "label", "RGC", "GLIP", "OPC", "OL", "ASTRO", "EPEN",
          "NEURON", "NON_NEURAL", "PROLIF_sig"]:
    A.obs[c] = ann[c].values
A.obs["leiden"] = A.obs["leiden"].astype(str).astype("category")
A.obs["label"] = A.obs["label"].astype("category")
A.obs["leiden_harmony"] = har["leiden"].astype(str).astype("category").values

A.obsm["X_pca"] = pca
A.obsm["X_harmony"] = hz
A.obsm["X_umap"] = ann[["UMAP1", "UMAP2"]].values
A.obsm["X_umap_harmony"] = har[["UMAP1", "UMAP2"]].values
A.obsm["jessa_signatures"] = sig.values
A.uns["jessa_signature_names"] = list(sig.columns)
A.uns["provenance"] = (
    "CellRanger 3.1.0 hg19+mm10 filtered matrices; human cells = >75% hg19 UMIs; "
    "QC per Jessa/Krug methods (<=50% MT, >=500 genes & UMIs, +/-2SD genes & UMIs per sample); "
    "CP10K + natural log; 2000 HVGs (seurat, batch_key=sample); 50 PCs. "
    "X = log-normalised; layers['counts'] = raw integer counts. "
    "X_umap/leiden = uncorrected merge (matches published methods). "
    "X_umap_harmony/leiden_harmony = Harmony on 'sample' -- SEE CAVEAT: samples are "
    "nested within genotype, so this correction also removes the genotype effect and "
    "should NOT be used for genotype comparisons."
)
A.write_h5ad(f"{OUT}/{tag}_merged.h5ad")
print(f"wrote {tag}_merged.h5ad", A.shape)
print(A.obs.groupby(["genotype", "replicate"], observed=True).size().to_string())
