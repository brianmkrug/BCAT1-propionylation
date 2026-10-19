"""DPT on the 5-library pooled manifold, rooted at max-RGC; representation diagnostic."""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
W = REPO_ROOT + "/scXenograft/work"
CLS = ["RGC","GLIP","OPC","OL","ASTRO","EPEN","NEURON"]

emb = pd.read_csv(f"{W}/pooled5_embed.csv.gz", index_col=0)
Z = np.load(f"{W}/pooled5_harmony.npz", allow_pickle=True)["X_harmony"]
meta = pd.read_csv(f"{W}/hsj019_metascores.csv.gz", index_col=0).loc[emb.index]
A = ad.AnnData(np.zeros((len(emb),1), np.float32), obs=emb.copy())
A.obsm["X_harmony"] = Z
for c in CLS + ["PROLIF_sig"]: A.obs[c] = meta[c].values
sc.pp.neighbors(A, n_neighbors=15, use_rep="X_harmony", random_state=0)
sc.tl.diffmap(A, n_comps=15)
cm = A.obs.groupby("leiden", observed=True)[CLS].mean(); clab = cm.idxmax(1)
A.obs["cluster_label"] = [clab[c] for c in A.obs.leiden]
root = A.obs.loc[A.obs.cluster_label=="RGC","RGC"].idxmax()
A.uns["iroot"] = int(np.where(A.obs_names==root)[0][0])
sc.tl.dpt(A, n_branchings=0)
emb["dpt"] = A.obs.dpt_pseudotime.values
for c in CLS + ["PROLIF_sig"]: emb[c] = meta[c].values
emb.to_csv(f"{W}/pooled5_dpt.csv.gz")
print(f"root {root} (RGC {A.obs.loc[root,'RGC']:.2f})")
print("\npseudotime per genotype:")
print(emb.groupby("genotype", observed=True).dpt.describe()[["mean","50%","count"]].round(3).to_string())
print("\npseudotime per KO library (does r2 sit differently?):")
print(emb[emb.genotype=="K27M-KO"].groupby("sample").dpt.describe()[["mean","50%","count"]].round(3).to_string())
