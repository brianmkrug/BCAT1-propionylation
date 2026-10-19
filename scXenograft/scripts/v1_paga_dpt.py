"""PAGA connectivity check, then DPT rooted at the max-RGC cell.

Run per genotype subset on the within-genotype Harmony space:
  k27m  = HSJ019 K27M r1+r2
  ko13  = HSJ019 K27M-KO r1+r3 (r2 excluded)

PAGA is computed and reported BEFORE DPT so the RGC->OL and RGC->ASTRO paths can
be checked for connectivity rather than assumed.
Usage: python3 v1_paga_dpt.py k27m|ko13
"""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import sys
import numpy as np, pandas as pd, scanpy as sc, anndata as ad

W = REPO_ROOT + "/scXenograft/work"
which = sys.argv[1]
CLS = ["RGC", "GLIP", "OPC", "OL", "ASTRO", "EPEN", "NEURON"]

e = pd.read_csv(f"{W}/hsj019_{which}_annotated.csv.gz", index_col=0)
H = np.load(f"{W}/hsj019_{which}_harmony.npz")["X_harmony"]

A = ad.AnnData(np.zeros((len(e), 1), dtype=np.float32), obs=e.copy())
A.obsm["X_harmony"] = H
A.obs["leiden"] = A.obs["leiden_harmony"].astype(str).astype("category")
sc.pp.neighbors(A, n_neighbors=15, use_rep="X_harmony", random_state=0)

# --- label each cluster by its dominant lineage meta-score ---
cm = A.obs.groupby("leiden", observed=True)[CLS].mean()
clab = cm.idxmax(1)
A.obs["cluster_label"] = [f"{c}:{clab[c]}" for c in A.obs.leiden]
A.obs["cluster_label"] = A.obs["cluster_label"].astype("category")

# --- PAGA on the labelled clusters ---
sc.tl.paga(A, groups="leiden")
conn = pd.DataFrame(A.uns["paga"]["connectivities"].todense(),
                    index=cm.index, columns=cm.index)
conn.index = [f"{c}:{clab[c]}" for c in conn.index]
conn.columns = conn.index

rgc = [c for c in conn.index if c.endswith(":RGC")]
ol = [c for c in conn.index if c.endswith(":OL")]
astro = [c for c in conn.index if c.endswith(":ASTRO")]
print(f"=== {which}: {A.n_obs} cells, {len(cm)} clusters ===")
print("cluster labels:", dict(zip(cm.index, clab.values)))
print(f"\nRGC clusters {rgc}\nOL clusters {ol}\nASTRO clusters {astro}")

def path_strength(src, dst, thr=0.05):
    """Is dst reachable from src through the PAGA graph above threshold?"""
    adj = conn.values > thr
    idx = {n: i for i, n in enumerate(conn.index)}
    seen, front, steps = set(idx[s] for s in src), set(idx[s] for s in src), 0
    while front and steps < 20:
        nxt = set()
        for f in front:
            nxt |= {j for j in np.where(adj[f])[0] if j not in seen}
        seen |= nxt; front = nxt; steps += 1
        if any(idx[d] in seen for d in dst):
            return steps
    return None

for name, dst in [("OL", ol), ("ASTRO", astro)]:
    if not rgc or not dst:
        print(f"  RGC -> {name}: MISSING endpoint cluster"); continue
    s = path_strength(rgc, dst)
    direct = conn.loc[rgc, dst].values.max()
    print(f"  RGC -> {name}: connected in {s} PAGA step(s)" if s else
          f"  RGC -> {name}: NOT CONNECTED above 0.05",
          f"| strongest direct edge {direct:.3f}")

# --- diffusion map + DPT rooted at max-RGC cell inside an RGC cluster ---
sc.tl.diffmap(A, n_comps=15)
rgc_cells = A.obs.cluster_label.astype(str).str.endswith(":RGC")
if rgc_cells.sum() == 0:
    raise SystemExit("no RGC-labelled cluster; not rooting arbitrarily")
root_name = A.obs.loc[rgc_cells, "RGC"].idxmax()
A.uns["iroot"] = int(np.where(A.obs_names == root_name)[0][0])
sc.tl.dpt(A, n_branchings=0)
print(f"\nroot cell: {root_name} (RGC score {A.obs.loc[root_name,'RGC']:.2f}, "
      f"cluster {A.obs.loc[root_name,'cluster_label']})")

out = A.obs.copy()
out["dpt_pseudotime"] = A.obs.dpt_pseudotime
for i in range(3):
    out[f"DC{i+1}"] = A.obsm["X_diffmap"][:, i+1]
out.to_csv(f"{W}/hsj019_{which}_dpt.csv.gz")
conn.round(3).to_csv(f"{W}/hsj019_{which}_paga_connectivities.csv")
np.savez_compressed(f"{W}/hsj019_{which}_paga_pos.npz",
                    conn=conn.values, names=np.array(conn.index, dtype=object))

print("\nmean pseudotime by lineage label:")
print(out.groupby("label").dpt_pseudotime.agg(["mean", "count"]).round(3).sort_values("mean").to_string())
print("\nSpearman of pseudotime vs meta-scores:")
from scipy.stats import spearmanr
for c in CLS + ["PROLIF_sig", "BCAT1"]:
    if c in out.columns:
        r, p = spearmanr(out.dpt_pseudotime, out[c])
        print(f"  {c:<12} rho {r:+.3f}  p {p:.2e}")
