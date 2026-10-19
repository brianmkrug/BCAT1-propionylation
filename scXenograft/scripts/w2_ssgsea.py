"""Single-sample GSEA (Barbie et al. 2009) per cell, vectorised.

The standard implementation walks the ranked gene list accumulating
  ES = sum_j [ P_w(j) - P_no(j) ]
which is O(genes) per cell per set. Because both P_w and P_no are step
functions, the sum has a closed form in the positions of the set genes only:

  ES = sum_i w_i (G - p_i) / sum_i w_i
       - [ G(G+1)/2 - sum_i (G - p_i) ] / (G - k)

with p_i the 0-based position of set gene i in the descending-ranked list,
w_i = x_i^alpha (alpha = 0.25), G genes total, k genes in the set. That is
O(k) per cell per set, so one argsort per cell chunk serves the whole panel.
`--validate` checks this against a literal random-walk implementation.

Usage: python3 w2_ssgsea.py k27m|ko13 [--validate]
"""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import sys
import numpy as np, pandas as pd, scanpy as sc

W = REPO_ROOT + "/scXenograft/work"
ALPHA = 0.25
CHUNK = 1500

def read_gmt(p):
    d = {}
    for ln in open(p):
        f = ln.rstrip("\n").split("\t")
        d[f[0]] = [x for x in f[2:] if x]
    return d

def es_closed(Xc, ranks, idx):
    """Xc (c,G) dense lognorm; ranks (c,G) 0-based descending position; idx set indices."""
    G = Xc.shape[1]; k = len(idx)
    p = ranks[:, idx].astype(np.float64)          # (c,k)
    w = np.power(np.abs(Xc[:, idx]), ALPHA, dtype=np.float64)
    tail = G - p
    Wsum = w.sum(1)
    Wsum[Wsum == 0] = np.nan
    pos = (w * tail).sum(1) / Wsum
    neg = (G * (G + 1) / 2 - tail.sum(1)) / (G - k)
    return pos - neg

def es_walk(x, idx):
    """Literal random walk for one cell — reference implementation."""
    G = len(x); k = len(idx)
    order = np.argsort(-x, kind="stable")
    ind = np.zeros(G, bool); ind[idx] = True
    ind_o = ind[order]
    w = np.power(np.abs(x[order]), ALPHA) * ind_o
    Pw = np.cumsum(w) / w.sum()
    Pno = np.cumsum(~ind_o) / (G - k)
    return (Pw - Pno).sum()

which = sys.argv[1]
validate = "--validate" in sys.argv
sets = read_gmt(f"{W}/metabolic_panel.gmt")

e = pd.read_csv(f"{W}/hsj019_{which}_annotated.csv.gz", index_col=0)
A = sc.read_h5ad(f"{W}/hsj019_norm.h5ad")[e.index]
genes = pd.Index(A.var_names)
sidx = {}
for k, v in sets.items():
    ii = genes.get_indexer(v)
    ii = np.unique(ii[ii >= 0])
    if len(ii) >= 5:
        sidx[k] = ii
print(f"{which}: {A.n_obs} cells x {A.n_vars} genes; "
      f"{len(sidx)}/{len(sets)} sets with >=5 genes detected")
print("gene recovery:", {k: f"{len(v)}/{len(sets[k])}" for k, v in list(sidx.items())[:4]}, "...")

res = np.zeros((A.n_obs, len(sidx)), dtype=np.float64)
names = list(sidx)
for s in range(0, A.n_obs, CHUNK):
    Xc = np.asarray(A.X[s:s+CHUNK].todense(), dtype=np.float32)
    order = np.argsort(-Xc, axis=1, kind="stable")
    ranks = np.empty_like(order)
    np.put_along_axis(ranks, order,
                      np.broadcast_to(np.arange(Xc.shape[1]), Xc.shape), axis=1)
    for j, k in enumerate(names):
        res[s:s+CHUNK, j] = es_closed(Xc, ranks, sidx[k])
    if validate and s == 0:
        print("\nvalidation (closed form vs literal random walk):")
        for k in names[:3]:
            for c in [0, 7, 123]:
                a = es_closed(Xc[c:c+1], ranks[c:c+1], sidx[k])[0]
                b = es_walk(Xc[c].astype(np.float64), sidx[k])
                print(f"  {k[:38]:<38} cell {c:<4} closed {a:+.6f}  walk {b:+.6f}"
                      f"  diff {abs(a-b):.2e}")
    print(f"  cells {s}-{min(s+CHUNK, A.n_obs)} done", flush=True)

df = pd.DataFrame(res, index=e.index, columns=names)
df.to_csv(f"{W}/hsj019_{which}_ssgsea_raw.csv.gz")
# standard ssGSEA normalisation: divide by the range of ES across all cells
norm = df / (df.max().max() - df.min().min())
norm.to_csv(f"{W}/hsj019_{which}_ssgsea_norm.csv.gz")
z = (df - df.mean()) / df.std()
z.to_csv(f"{W}/hsj019_{which}_ssgsea_z.csv.gz")
print(f"\nwrote raw / range-normalised / z-scored scores for {which}")
print(df.describe().T[["mean", "std", "min", "max"]].round(3).to_string())
