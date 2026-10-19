"""Is the top of the ranking just the cell cycle?

The strongest K27M hits (one-carbon, pyrimidine, purine, folate) are nucleotide
metabolism, and pseudotime tracks proliferation almost perfectly there
(PROLIF_sig rho -0.86). Two controls:

  A. partial correlation of each pathway with pseudotime, controlling S+G2M
  B. recompute eta2 within G1 cells only (cycling cells removed entirely)
"""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import numpy as np, pandas as pd
from scipy.stats import spearmanr, rankdata

W = REPO_ROOT + "/scXenograft/work"
NBIN = 20

def eta2(v, g, nb):
    ok = np.isfinite(v) & np.isfinite(g)
    v, g = v[ok], g[ok]
    gm = v.mean(); ss_t = ((v-gm)**2).sum(); ss_b = 0.0
    for b in range(nb):
        m = g == b
        if m.sum() > 1: ss_b += m.sum()*(v[m].mean()-gm)**2
    return ss_b/ss_t if ss_t > 0 else np.nan

def partial(x, y, Z):
    ok = np.isfinite(x) & np.isfinite(y) & np.all(np.isfinite(Z), 1)
    rx, ry = rankdata(x[ok]), rankdata(y[ok])
    RZ = np.c_[np.ones(ok.sum()), np.apply_along_axis(rankdata, 0, Z[ok])]
    def res(a): return a - RZ @ np.linalg.lstsq(RZ, a, rcond=None)[0]
    return np.corrcoef(res(rx), res(ry))[0, 1]

for which, ttl in [("k27m", "K27M"), ("ko13", "KO r1+r3")]:
    d = pd.read_csv(f"{W}/hsj019_{which}_dpt.csv.gz", index_col=0)
    S = pd.read_csv(f"{W}/hsj019_{which}_allmetab_resid.csv.gz", index_col=0).loc[d.index]
    var = pd.read_csv(f"{W}/hsj019_{which}_pathway_variability.csv", index_col=0)
    pt = d.dpt_pseudotime.values
    CC = np.c_[d.S_score.values, d.G2M_score.values]
    g1 = (d.phase == "G1").values
    pt_bin = pd.qcut(pt, NBIN, labels=False, duplicates="drop")
    pt_bin_g1 = pd.qcut(pt[g1], NBIN, labels=False, duplicates="drop")

    rows = []
    for c in S.columns:
        v = S[c].values
        rows.append(dict(pathway=c,
                         rho=var.rho.get(c), rho_partial_cc=partial(pt, v, CC),
                         eta2=var.eta2.get(c),
                         eta2_G1only=eta2(v[g1], pt_bin_g1, NBIN),
                         lib_consistency=var.lib_consistency.get(c),
                         z_null=var.z_null.get(c)))
    r = pd.DataFrame(rows).set_index("pathway")
    r["eta2_ratio_G1"] = r.eta2_G1only / r.eta2
    r.round(4).to_csv(f"{W}/hsj019_{which}_pathway_cellcycle_control.csv")

    print(f"\n=== {ttl}: does the signal survive cell-cycle control? ===")
    sel = r[(r.z_null > 3) & (r.lib_consistency > 0.5)].sort_values("eta2", ascending=False)
    print(sel[["rho", "rho_partial_cc", "eta2", "eta2_G1only", "eta2_ratio_G1"]]
          .round(3).to_string())
    print(f"\n  G1-only cells: {g1.sum()} of {len(g1)}")
    # what survives in G1 only
    r2 = r[(r.eta2_G1only > 0.08) & (r.lib_consistency > 0.5)].sort_values(
        "eta2_G1only", ascending=False)
    print(f"\n  top pathways by eta2 computed in G1 cells ONLY:")
    print(r2[["eta2", "eta2_G1only", "rho", "rho_partial_cc"]].head(12).round(3).to_string())
