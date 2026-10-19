"""Score the full metabolic pathway compendium with depth-robust control-bin
scoring, then residualise on detection depth.

Sets: the curated human KEGG-78 list + Hallmark metabolic sets.
Also scores size-matched random gene sets (10 per size stratum) as a null.

Usage: python3 x1_score_all.py k27m|ko13
"""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import sys
import numpy as np, pandas as pd, scipy.sparse as sp, scanpy as sc

B = REPO_ROOT + ""
W = f"{B}/scXenograft/work"
N_CTRL, N_BINS, N_RAND = 50, 25, 10
which = sys.argv[1]
rng = np.random.default_rng(0)

def read_gmt(p, prefix=""):
    d = {}
    for ln in open(p):
        f = ln.rstrip("\n").split("\t")
        g = [x for x in f[2:] if x]
        if g: d[prefix + f[0]] = g
    return d

kegg = read_gmt(f"{B}/xHSJ019_K27M_vs_KO_analysis/KEGG_metabolism_expanded.gmt", "KEGG_")
hall = read_gmt(f"{B}/signatures/h.all.v7.5.symbols.gmt")
HALL_METAB = ["HALLMARK_CHOLESTEROL_HOMEOSTASIS", "HALLMARK_FATTY_ACID_METABOLISM",
              "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_GLYCOLYSIS",
              "HALLMARK_ADIPOGENESIS", "HALLMARK_PEROXISOME", "HALLMARK_MTORC1_SIGNALING",
              "HALLMARK_HYPOXIA", "HALLMARK_BILE_ACID_METABOLISM",
              "HALLMARK_XENOBIOTIC_METABOLISM", "HALLMARK_HEME_METABOLISM"]
sets = dict(kegg)
for h in HALL_METAB:
    if h in hall: sets["HM_" + h.replace("HALLMARK_", "")] = hall[h]
print(f"{len(sets)} metabolic pathways loaded ({len(kegg)} curated KEGG, "
      f"{len(sets)-len(kegg)} Hallmark)")

e = pd.read_csv(f"{W}/hsj019_{which}_dpt.csv.gz", index_col=0)
A = sc.read_h5ad(f"{W}/hsj019_norm.h5ad")[e.index]
X = A.X.tocsc() if sp.issparse(A.X) else sp.csc_matrix(A.X)
genes = pd.Index(A.var_names); NG = X.shape[1]

avg = np.asarray(X.mean(0)).ravel()
order = np.argsort(avg)
bins = np.zeros(NG, dtype=int)
bins[order] = np.minimum((np.arange(NG)/NG*N_BINS).astype(int), N_BINS-1)
members = {b: np.where(bins == b)[0] for b in range(N_BINS)}

def build(idx_map):
    names, rs, cs, vs, rc, cc_, vc, meta = [], [], [], [], [], [], [], []
    for nm, idx in idx_map.items():
        if len(idx) < 8:      # 8-gene floor: below this, single-cell dropout dominates
            continue
        i = len(names); names.append(nm); k = len(idx)
        meta.append((nm, k))
        rs += [i]*k; cs += list(idx); vs += [1.0/k]*k
        ctrl = []
        for b, cnt in zip(*np.unique(bins[idx], return_counts=True)):
            pool = np.setdiff1d(members[b], idx)
            take = min(len(pool), cnt*N_CTRL)
            if take: ctrl.append(rng.choice(pool, take, replace=False))
        ctrl = np.concatenate(ctrl); kc = len(ctrl)
        rc += [i]*kc; cc_ += list(ctrl); vc += [1.0/kc]*kc
    S = sp.csr_matrix((vs, (rs, cs)), shape=(len(names), NG))
    C = sp.csr_matrix((vc, (rc, cc_)), shape=(len(names), NG))
    sco = np.asarray((X @ S.T).todense()) - np.asarray((X @ C.T).todense())
    return pd.DataFrame(sco, index=e.index, columns=names), pd.DataFrame(meta, columns=["pathway", "n_genes_detected"])

real_idx = {}
for nm, gl in sets.items():
    ii = genes.get_indexer(gl); ii = np.unique(ii[ii >= 0])
    if len(ii) >= 8: real_idx[nm] = ii
df, meta = build(real_idx)

# size-matched random null
szs = sorted(set(len(v) for v in real_idx.values()))
rand_idx = {}
for sz in szs[::max(1, len(szs)//12)]:
    for r in range(N_RAND):
        rand_idx[f"RND_{sz}_{r}"] = rng.choice(NG, sz, replace=False)
rdf, _ = build(rand_idx)

# depth residualisation: remove quadratic n_genes trend from every score
ng = e.n_genes.values.astype(float)
M = np.c_[np.ones_like(ng), ng, ng**2]
def resid(D):
    beta, *_ = np.linalg.lstsq(M, D.values, rcond=None)
    return pd.DataFrame(D.values - M @ beta, index=D.index, columns=D.columns)

df.to_csv(f"{W}/hsj019_{which}_allmetab_score.csv.gz")
resid(df).to_csv(f"{W}/hsj019_{which}_allmetab_resid.csv.gz")
resid(rdf).to_csv(f"{W}/hsj019_{which}_allmetab_null_resid.csv.gz")
meta.to_csv(f"{W}/hsj019_{which}_allmetab_sizes.csv", index=False)
print(f"{which}: scored {df.shape[1]} pathways + {rdf.shape[1]} random sets "
      f"over {df.shape[0]} cells")
from scipy.stats import spearmanr
r = [abs(spearmanr(ng, df[c]).statistic) for c in df.columns]
rr = [abs(spearmanr(ng, resid(df)[c]).statistic) for c in df.columns]
print(f"  |rho| with n_genes: raw scores mean {np.mean(r):.3f} max {np.max(r):.3f}; "
      f"after residualisation mean {np.mean(rr):.3f} max {np.max(rr):.3f}")
