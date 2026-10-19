"""Which metabolic pathways vary most along the trajectory?

Metrics per pathway, on depth-residualised control-bin scores:
  eta2   fraction of score variance explained by pseudotime bin (captures
         monotonic and non-monotonic trends alike)
  amp    max-min of the binned profile, in SD units of that pathway's score
  rho    Spearman with pseudotime (direction)
  libr   correlation between the binned profiles of the two libraries
         (replicate consistency: a pathway that varies for real should trace
         the same shape in each library independently)
  z_null eta2 expressed as a z-score against 120 size-matched random gene sets

Usage: python3 x2_variability.py
"""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import numpy as np, pandas as pd
from scipy.stats import spearmanr

W = REPO_ROOT + "/scXenograft/work"
NBIN = 20

def profile(v, pt_bin, nb=NBIN, minn=15):
    out = np.full(nb, np.nan)
    for b in range(nb):
        m = pt_bin == b
        if m.sum() >= minn: out[b] = np.nanmean(v[m])
    return out

def eta2(v, pt_bin, nb=NBIN):
    ok = np.isfinite(v)
    v, g = v[ok], pt_bin[ok]
    gm = v.mean()
    ss_t = ((v - gm)**2).sum()
    ss_b = 0.0
    for b in range(nb):
        m = g == b
        if m.sum() > 1: ss_b += m.sum() * (v[m].mean() - gm)**2
    return ss_b / ss_t if ss_t > 0 else np.nan

RES = {}
for which, ttl in [("k27m", "K27M"), ("ko13", "KO r1+r3")]:
    d = pd.read_csv(f"{W}/hsj019_{which}_dpt.csv.gz", index_col=0)
    S = pd.read_csv(f"{W}/hsj019_{which}_allmetab_resid.csv.gz", index_col=0).loc[d.index]
    N = pd.read_csv(f"{W}/hsj019_{which}_allmetab_null_resid.csv.gz", index_col=0).loc[d.index]
    sizes = pd.read_csv(f"{W}/hsj019_{which}_allmetab_sizes.csv").set_index("pathway")

    pt = d.dpt_pseudotime.values
    pt_bin = pd.qcut(pt, NBIN, labels=False, duplicates="drop")
    libs = sorted(d["sample"].unique())
    lm = {l: (d["sample"] == l).values for l in libs}

    null_e = np.array([eta2(N[c].values, pt_bin) for c in N.columns])
    mu, sd = null_e.mean(), null_e.std()

    rows = []
    for c in S.columns:
        v = S[c].values
        sdv = np.nanstd(v)
        prof = profile(v, pt_bin)
        e = eta2(v, pt_bin)
        p1 = profile(v[lm[libs[0]]], pt_bin[lm[libs[0]]])
        p2 = profile(v[lm[libs[1]]], pt_bin[lm[libs[1]]])
        ok = np.isfinite(p1) & np.isfinite(p2)
        lr = np.corrcoef(p1[ok], p2[ok])[0, 1] if ok.sum() >= 6 else np.nan
        rows.append(dict(pathway=c, n_genes=int(sizes.n_genes_detected.get(c, np.nan)),
                         eta2=e, z_null=(e - mu)/sd,
                         amp=(np.nanmax(prof)-np.nanmin(prof))/sdv,
                         rho=spearmanr(pt[np.isfinite(v)], v[np.isfinite(v)]).statistic,
                         lib_consistency=lr))
    r = pd.DataFrame(rows).set_index("pathway").sort_values("eta2", ascending=False)
    RES[which] = r
    r.round(4).to_csv(f"{W}/hsj019_{which}_pathway_variability.csv")
    print(f"\n=== {ttl}: null eta2 = {mu:.4f} +/- {sd:.4f} "
          f"(120 size-matched random sets) ===")
    print(f"pathways above null z>3 AND library-consistent (r>0.5): "
          f"{((r.z_null > 3) & (r.lib_consistency > 0.5)).sum()} of {len(r)}\n")
    top = r[(r.z_null > 3) & (r.lib_consistency > 0.5)].head(20)
    print(top[["n_genes", "eta2", "z_null", "amp", "rho", "lib_consistency"]]
          .round(3).to_string())

# ---- genotype comparison ----
k, o = RES["k27m"], RES["ko13"]
comp = k[["eta2", "rho", "z_null", "lib_consistency"]].join(
    o[["eta2", "rho", "z_null", "lib_consistency"]], lsuffix="_K27M", rsuffix="_KO")
comp["delta_eta2"] = comp.eta2_KO - comp.eta2_K27M
comp = comp.sort_values("delta_eta2")
comp.round(4).to_csv(f"{W}/pathway_variability_genotype_comparison.csv")
print("\n\n=== Pathways whose trajectory-variation differs most between genotypes ===")
print("\n-- much more variable in K27M --")
print(comp.head(10)[["eta2_K27M", "eta2_KO", "rho_K27M", "rho_KO", "delta_eta2"]].round(3).to_string())
print("\n-- much more variable in KO --")
print(comp.tail(10)[["eta2_K27M", "eta2_KO", "rho_K27M", "rho_KO", "delta_eta2"]].round(3).to_string())
