#!/usr/bin/env python3
"""
Classify DIPGXIII K27M CRISPR screen hits by DepMap dependency class.

Ranks the metabolism-focused screen (beta_correlation.xlsx) by most-negative
K27M beta, takes the top N, and annotates each gene with:
  - DepMap "Common Essential"  (CRISPRInferredCommonEssentials.csv, 24Q4)
  - DepMap "Strongly Selective" (normLRT > 100, recomputed from CRISPRGeneEffect)
  - "Other dependency"          (>=1 line with dependency probability > 0.5)
  - "Not dependent"             (0 lines dependent)
  - "Not in DepMap"

normLRT = 2 * (loglik of Azzalini skew-t fit - loglik of normal fit) to the
gene's cross-cell-line Chronos gene-effect distribution. This is the statistic
the DepMap portal uses to flag selective dependencies (threshold 100).

Data: DepMap Public 24Q4 (figshare 27993248).
"""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import json
import numpy as np
import pandas as pd
from scipy import stats, optimize

DEPMAP = "/tmp/depmap"
SCREEN = REPO_ROOT + "/CRISPR screen/beta_correlation.xlsx"
OUT = REPO_ROOT + "/CRISPR screen/depmap_classification"
TOP_N = 200
NORMLRT_THRESH = 100.0
DEP_PROB_THRESH = 0.5

# --- Azzalini skew-t log-likelihood -----------------------------------------


def _skewt_nll(p, x):
    xi, log_om, alpha, log_nu = p
    om = np.exp(log_om)
    nu = np.exp(log_nu)
    if not np.isfinite(om) or om <= 0 or nu <= 0.5 or nu > 1e6:
        return 1e12
    z = (x - xi) / om
    lg = stats.t.logpdf(z, nu)
    w = alpha * z * np.sqrt((nu + 1.0) / (z * z + nu))
    lcdf = stats.t.logcdf(w, nu + 1.0)
    ll = np.log(2.0) - log_om + lg + lcdf
    if not np.all(np.isfinite(ll)):
        return 1e12
    return -ll.sum()


def norm_lrt(x):
    """2*(skew-t loglik - normal loglik) for a 1-D sample."""
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if x.size < 50:
        return np.nan
    mu, sd = x.mean(), x.std(ddof=0)
    if sd <= 0:
        return np.nan
    ll_norm = stats.norm.logpdf(x, mu, sd).sum()

    best = -np.inf
    g1 = stats.skew(x)
    starts = [
        (mu, np.log(sd), 0.0, np.log(10.0)),
        (mu, np.log(sd), np.sign(g1) * 3.0 if g1 else 1.0, np.log(6.0)),
        (np.median(x), np.log(sd), np.sign(g1) * 8.0 if g1 else -3.0, np.log(20.0)),
        (np.median(x), np.log(sd * 0.7), -5.0, np.log(4.0)),
        (np.median(x), np.log(sd * 0.7), 5.0, np.log(4.0)),
    ]
    for s0 in starts:
        try:
            r = optimize.minimize(_skewt_nll, s0, args=(x,), method="Nelder-Mead",
                                  options={"maxiter": 4000, "xatol": 1e-6, "fatol": 1e-6})
            if np.isfinite(r.fun) and -r.fun > best:
                best = -r.fun
        except Exception:
            pass
    if not np.isfinite(best):
        return np.nan
    return 2.0 * (best - ll_norm)


# --- load --------------------------------------------------------------------

def strip_id(c):
    return c.split(" (")[0]


def main():
    scr = pd.read_excel(SCREEN)
    scr = scr.sort_values("K27M.beta").reset_index(drop=True)
    top = scr.head(TOP_N).copy()

    # validation genes: canonical common-essentials and canonical selectives
    validation = ["RPL5", "POLR2A", "RAN", "PSMB2", "EIF3B",            # common essential
                  "CTNNB1", "KRAS", "BRAF", "SOX10", "PAX8", "MYB",     # strongly selective
                  "EGFR", "NRAS", "MDM2", "TP53", "BCAT1", "BCAT2"]
    want = sorted(set(top["Gene"]).union(validation))

    ce = pd.read_csv(f"{DEPMAP}/CRISPRInferredCommonEssentials.csv")
    ce_genes = set(ce["Essentials"].map(strip_id))

    print("reading CRISPRGeneEffect ...", flush=True)
    hdr = pd.read_csv(f"{DEPMAP}/CRISPRGeneEffect.csv", nrows=0)
    cols = list(hdr.columns)
    sym = {c: strip_id(c) for c in cols[1:]}
    use = [cols[0]] + [c for c in cols[1:] if sym[c] in want]
    eff = pd.read_csv(f"{DEPMAP}/CRISPRGeneEffect.csv", usecols=use, index_col=0)
    eff.columns = [strip_id(c) for c in eff.columns]

    print("reading CRISPRGeneDependency ...", flush=True)
    hdr2 = pd.read_csv(f"{DEPMAP}/CRISPRGeneDependency.csv", nrows=0)
    cols2 = list(hdr2.columns)
    use2 = [cols2[0]] + [c for c in cols2[1:] if strip_id(c) in want]
    dep = pd.read_csv(f"{DEPMAP}/CRISPRGeneDependency.csv", usecols=use2, index_col=0)
    dep.columns = [strip_id(c) for c in dep.columns]

    n_lines = eff.shape[0]
    print(f"DepMap lines: {n_lines}; genes matched: {eff.shape[1]}/{len(want)}", flush=True)

    rows = []
    for i, g in enumerate(want):
        if g not in eff.columns:
            rows.append(dict(Gene=g, in_depmap=False))
            continue
        e = eff[g].dropna().values
        d = dep[g].dropna().values if g in dep.columns else np.array([])
        n_dep = int((d > DEP_PROB_THRESH).sum())
        rows.append(dict(
            Gene=g, in_depmap=True,
            n_lines=len(e), mean_gene_effect=float(np.mean(e)),
            median_gene_effect=float(np.median(e)),
            n_dep_lines=n_dep, frac_dep_lines=n_dep / len(d) if len(d) else np.nan,
            normLRT=norm_lrt(e),
        ))
        if (i + 1) % 25 == 0:
            print(f"  {i+1}/{len(want)}", flush=True)

    ann = pd.DataFrame(rows)

    def classify(r):
        if not r["in_depmap"]:
            return "Not in DepMap"
        if r["Gene"] in ce_genes:
            return "Common essential"
        if np.isfinite(r["normLRT"]) and r["normLRT"] > NORMLRT_THRESH and r["n_dep_lines"] >= 1:
            return "Strongly selective"
        if r["n_dep_lines"] >= 1:
            return "Other dependency"
        return "Not dependent"

    ann["depmap_class"] = ann.apply(classify, axis=1)
    ann.to_csv(f"{OUT}/depmap_annotation_all.csv", index=False)

    merged = top.merge(ann, on="Gene", how="left")
    merged["depmap_class"] = merged["depmap_class"].fillna("Not in DepMap")
    merged["rank"] = np.arange(1, len(merged) + 1)
    merged.to_csv(f"{OUT}/top{TOP_N}_K27M_depmap_annotated.csv", index=False)

    print("\n--- validation genes ---")
    print(ann[ann.Gene.isin(validation)][
        ["Gene", "mean_gene_effect", "frac_dep_lines", "normLRT", "depmap_class"]
    ].to_string(index=False))
    print("\n--- top200 class counts ---")
    print(merged["depmap_class"].value_counts().to_string())


if __name__ == "__main__":
    main()
