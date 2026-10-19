"""DPT on the pooled Harmony manifold, rooted at the max-RGC cell.
Then quantify per-bin genotype representation and define the overlap zone
where the bin-conditional comparison is supported."""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
from scipy.stats import spearmanr

W = REPO_ROOT + "/scXenograft/work"
CLS = ["RGC", "GLIP", "OPC", "OL", "ASTRO", "EPEN", "NEURON"]

emb = pd.read_csv(f"{W}/pooled_embed.csv.gz", index_col=0)
Z = np.load(f"{W}/pooled_harmony.npz", allow_pickle=True)["X_harmony"]
meta = pd.read_csv(f"{W}/hsj019_metascores.csv.gz", index_col=0).loc[emb.index]

A = ad.AnnData(np.zeros((len(emb), 1), np.float32), obs=emb.copy())
A.obsm["X_harmony"] = Z
for c in CLS + ["PROLIF_sig"]:
    A.obs[c] = meta[c].values
sc.pp.neighbors(A, n_neighbors=15, use_rep="X_harmony", random_state=0)
sc.tl.diffmap(A, n_comps=15)

# label clusters, root in an RGC cluster
cm = A.obs.groupby("leiden", observed=True)[CLS].mean()
clab = cm.idxmax(1)
A.obs["cluster_label"] = [clab[c] for c in A.obs.leiden]
rgc_cells = A.obs.cluster_label == "RGC"
root = A.obs.loc[rgc_cells, "RGC"].idxmax()
A.uns["iroot"] = int(np.where(A.obs_names == root)[0][0])
sc.tl.dpt(A, n_branchings=0)

emb["dpt"] = A.obs.dpt_pseudotime.values
for c in CLS + ["PROLIF_sig"]:
    emb[c] = meta[c].values
emb.to_csv(f"{W}/pooled_dpt.csv.gz")
print(f"pooled DPT rooted at {root} (cluster {clab[emb.loc[root,'leiden']]}, "
      f"RGC {A.obs.loc[root,'RGC']:.2f})")

# pseudotime distribution per genotype
print("\npseudotime per genotype:")
print(emb.groupby("genotype", observed=True).dpt.describe()[
    ["mean", "25%", "50%", "75%", "count"]].round(3).to_string())

# per-bin representation -> overlap zone
NB = 20
emb["bin"] = pd.qcut(emb.dpt, NB, labels=False, duplicates="drop")
rep = emb.groupby("bin", observed=True).apply(
    lambda d: pd.Series({
        "n_K27M": int((d.genotype == "K27M").sum()),
        "n_KO": int((d.genotype == "K27M-KO").sum()),
        "n_KO_r1": int((d["sample"].str.endswith("ko_r1")).sum()),
        "n_KO_r3": int((d["sample"].str.endswith("ko_r3")).sum()),
        "dpt_mid": d.dpt.mean()}))
rep["min_genotype_n"] = rep[["n_K27M", "n_KO"]].min(1)
# overlap zone: >=40 cells of each genotype in the bin
rep["overlap"] = rep.min_genotype_n >= 40
rep.to_csv(f"{W}/pooled_bin_representation.csv")
print("\nper-bin genotype representation (bin = pseudotime quantile):")
print(rep[["dpt_mid", "n_K27M", "n_KO", "n_KO_r1", "n_KO_r3", "overlap"]]
      .round(3).to_string())
ov = rep[rep.overlap]
print(f"\noverlap zone: bins {list(ov.index)}  "
      f"(pseudotime {ov.dpt_mid.min():.2f}-{ov.dpt_mid.max():.2f}); "
      f"{int(ov.n_K27M.sum())} K27M + {int(ov.n_KO.sum())} KO cells")
