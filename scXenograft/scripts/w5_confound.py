"""Is the single-cell ssGSEA gradient real biology or detection depth?

All panel signatures track each other almost exactly, which is the signature of
a shared technical driver. Tests:
  1. correlation of each score with n_genes / n_umi / pct_mt
  2. correlation between unrelated signatures (should be low if biological)
  3. the same, after regressing out n_genes
  4. random gene sets of matched size as a null
"""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import numpy as np, pandas as pd
from scipy.stats import spearmanr
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

W = REPO_ROOT + "/scXenograft/work"
FIG = REPO_ROOT + "/scXenograft/figures"
SETS = [("k27m", "HSJ019 K27M"), ("ko13", "HSJ019 K27M-KO (r1+r3)")]

print("1) Correlation of ssGSEA scores with technical covariates (Spearman)\n")
store = {}
for w, ttl in SETS:
    d = pd.read_csv(f"{W}/hsj019_{w}_dpt.csv.gz", index_col=0)
    z = pd.read_csv(f"{W}/hsj019_{w}_ssgsea_z.csv.gz", index_col=0).loc[d.index]
    store[w] = (d, z)
    rows = []
    for k in z.columns:
        v = z[k].values; ok = np.isfinite(v)
        rows.append(dict(signature=k,
                         n_genes=spearmanr(d.n_genes.values[ok], v[ok]).statistic,
                         n_umi=spearmanr(d.n_umi.values[ok], v[ok]).statistic,
                         pct_mt=spearmanr(d.pct_mt.values[ok], v[ok]).statistic,
                         dpt=spearmanr(d.dpt_pseudotime.values[ok], v[ok]).statistic))
    r = pd.DataFrame(rows).set_index("signature")
    print(f"--- {ttl} ---")
    print(r.round(3).to_string())
    print(f"  MEAN |rho| with n_genes: {r.n_genes.abs().mean():.3f}"
          f"   with dpt: {r.dpt.abs().mean():.3f}\n")
    r.round(4).to_csv(f"{W}/ssgsea_{w}_technical_correlations.csv")

print("\n2) Pairwise correlation between signatures "
      "(if biology, unrelated pathways should NOT track)")
for w, ttl in SETS:
    d, z = store[w]
    cc = z.corr(method="spearman")
    iu = np.triu_indices_from(cc, k=1)
    print(f"  {ttl}: median pairwise rho {np.nanmedian(cc.values[iu]):.3f}, "
          f"CHOL_Hallmark vs CTX_glycolysis "
          f"{cc.loc['CHOL_Hallmark_cholesterol_homeostasis','CTX_KEGGc_glycolysis']:.3f}, "
          f"vs OxPhos {cc.loc['CHOL_Hallmark_cholesterol_homeostasis','CTX_KEGGc_oxidative_phosphorylation']:.3f}")

print("\n3) After regressing n_genes out of every score")
resid_store = {}
for w, ttl in SETS:
    d, z = store[w]
    x = d.n_genes.values.astype(float)
    X = np.c_[np.ones_like(x), x, x**2]
    R = {}
    for k in z.columns:
        v = z[k].values; ok = np.isfinite(v)
        beta, *_ = np.linalg.lstsq(X[ok], v[ok], rcond=None)
        r_ = np.full_like(v, np.nan); r_[ok] = v[ok] - X[ok] @ beta
        R[k] = r_
    R = pd.DataFrame(R, index=z.index)
    resid_store[w] = R
    R.to_csv(f"{W}/hsj019_{w}_ssgsea_z_resid.csv.gz")
    cc = R.corr(method="spearman")
    iu = np.triu_indices_from(cc, k=1)
    dpt_r = [spearmanr(d.dpt_pseudotime.values[np.isfinite(R[k])],
                       R[k].values[np.isfinite(R[k])]).statistic for k in R.columns]
    print(f"  {ttl}: median pairwise rho {np.nanmedian(cc.values[iu]):+.3f} "
          f"(was {np.nanmedian(z.corr(method='spearman').values[iu]):+.3f}); "
          f"mean |rho| with DPT {np.mean(np.abs(dpt_r)):.3f}")

print("\n4) Random gene sets of matched size as a null (K27M)")
d, z = store["k27m"]
import scanpy as sc
A = sc.read_h5ad(f"{W}/hsj019_norm.h5ad")[d.index]
rng = np.random.default_rng(0)
# reuse ssGSEA closed form on random sets
ALPHA = 0.25
sizes = [25, 74, 158]
out = {}
CH = 1500
for sz in sizes:
    out[f"RANDOM_{sz}"] = np.zeros(A.n_obs)
ridx = {f"RANDOM_{sz}": rng.choice(A.n_vars, sz, replace=False) for sz in sizes}
for s in range(0, A.n_obs, CH):
    Xc = np.asarray(A.X[s:s+CH].todense(), dtype=np.float32)
    order = np.argsort(-Xc, axis=1, kind="stable")
    ranks = np.empty_like(order)
    np.put_along_axis(ranks, order,
                      np.broadcast_to(np.arange(Xc.shape[1]), Xc.shape), axis=1)
    G = Xc.shape[1]
    for k, idx in ridx.items():
        p = ranks[:, idx].astype(np.float64)
        wgt = np.power(np.abs(Xc[:, idx]), ALPHA, dtype=np.float64)
        tail = G - p; Ws = wgt.sum(1); Ws[Ws == 0] = np.nan
        out[k][s:s+CH] = (wgt*tail).sum(1)/Ws - (G*(G+1)/2 - tail.sum(1))/(G-len(idx))
rnd = pd.DataFrame(out, index=d.index)
for k in rnd.columns:
    v = rnd[k].values; ok = np.isfinite(v)
    print(f"  {k}: rho with n_genes {spearmanr(d.n_genes.values[ok], v[ok]).statistic:+.3f}, "
          f"with DPT {spearmanr(d.dpt_pseudotime.values[ok], v[ok]).statistic:+.3f}")
rnd.to_csv(f"{W}/hsj019_k27m_ssgsea_random_null.csv.gz")
