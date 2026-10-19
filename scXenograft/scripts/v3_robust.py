"""Robustness checks on the DPT orderings.

1. Root sensitivity: re-root at 20 random cells drawn from the top 1% of RGC
   score and correlate the resulting pseudotimes with the reported one.
2. Replicate independence (KO): fit DPT on r3 alone and correlate with the
   joint r1+r3 ordering on the shared cells.
3. Connectivity of the diffusion graph (does the pseudotime gap reflect a
   disconnected component?).
"""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import sys
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
from scipy.stats import spearmanr
import scipy.sparse as sp

W = REPO_ROOT + "/scXenograft/work"
rng = np.random.default_rng(0)

def build(which, mask=None):
    e = pd.read_csv(f"{W}/hsj019_{which}_annotated.csv.gz", index_col=0)
    H = np.load(f"{W}/hsj019_{which}_harmony.npz")["X_harmony"]
    if mask is not None:
        H, e = H[mask.values], e[mask.values]
    A = ad.AnnData(np.zeros((len(e), 1), dtype=np.float32), obs=e.copy())
    A.obsm["X_harmony"] = H
    sc.pp.neighbors(A, n_neighbors=15, use_rep="X_harmony", random_state=0)
    sc.tl.diffmap(A, n_comps=15)
    return A

def dpt_from(A, root_name):
    A.uns["iroot"] = int(np.where(A.obs_names == root_name)[0][0])
    sc.tl.dpt(A, n_branchings=0)
    return A.obs.dpt_pseudotime.copy()

for which in ["k27m", "ko13"]:
    ref = pd.read_csv(f"{W}/hsj019_{which}_dpt.csv.gz", index_col=0)
    A = build(which)
    # graph connectivity
    ncomp, lab = sp.csgraph.connected_components(A.obsp["connectivities"],
                                                 directed=False)
    sizes = np.bincount(lab)
    print(f"\n=== {which} ===")
    print(f"kNN graph components: {ncomp} (sizes {sorted(sizes)[::-1][:5]})")

    # root sensitivity
    thr = np.percentile(A.obs.RGC, 99)
    pool = A.obs_names[A.obs.RGC >= thr]
    rhos = []
    for rn in rng.choice(pool, size=min(20, len(pool)), replace=False):
        pt = dpt_from(A, rn)
        ok = np.isfinite(pt) & np.isfinite(ref.dpt_pseudotime.values)
        rhos.append(spearmanr(pt[ok], ref.dpt_pseudotime.values[ok]).statistic)
    rhos = np.array(rhos)
    print(f"root sensitivity: Spearman vs reported ordering over 20 alternative "
          f"roots in top-1% RGC = {rhos.mean():.3f} "
          f"(min {rhos.min():.3f}, max {rhos.max():.3f})")

# replicate independence for KO: r3 alone
e = pd.read_csv(f"{W}/hsj019_ko13_annotated.csv.gz", index_col=0)
mask = e["sample"].str.endswith("ko_r3")
A3 = build("ko13", mask=mask)
rgc = A3.obs.cluster_label.astype(str).str.endswith(":RGC") if "cluster_label" in A3.obs else None
root3 = A3.obs.RGC.idxmax()
pt3 = dpt_from(A3, root3)
ref = pd.read_csv(f"{W}/hsj019_ko13_dpt.csv.gz", index_col=0)
shared = pt3.index
r = spearmanr(pt3.values, ref.loc[shared, "dpt_pseudotime"].values)
print(f"\nKO r3 alone vs joint r1+r3 ordering (n={len(shared)}): "
      f"Spearman rho {r.statistic:+.3f}, p {r.pvalue:.2e}")

# correlation of r3-only pseudotime with lineage scores
print("r3-only pseudotime vs meta-scores:")
for c in ["RGC", "OPC", "OL", "ASTRO", "PROLIF_sig", "BCAT1"]:
    rr = spearmanr(pt3.values, A3.obs[c].values)
    print(f"  {c:<12} rho {rr.statistic:+.3f}")
pd.DataFrame({"dpt_r3only": pt3}).to_csv(f"{W}/hsj019_ko_r3only_dpt.csv.gz")
