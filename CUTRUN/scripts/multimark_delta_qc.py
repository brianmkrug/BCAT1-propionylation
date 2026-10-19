#!/usr/bin/env python
"""
Track QC on the Δ bigWigs produced by run_multimark_tracks.sh.

Two things, both about whether the tracks are fit to feed the promoter analysis:

1. WIRING CHECK -- confirm each delta_<mark>_<rep>.log2.bw really was built from the
   intended depth-matching tier, by recomputing log2((shBCAT1+1)/(shNT+1)) from the
   candidate CPM tracks and seeing which one matches. The run log records the output
   name but not the inputs, so this is checked rather than assumed.

2. REPLICATE CONCORDANCE OF Δ AT PROMOTERS -- the acceptance gate this project already
   uses (rho > 0.3 AND > 60% same-direction, pre-registered in
   H3K4me3/reports/phase6_findings_report.md). Promoter windows come from
   H3K4me3/reference/core.bed (TSS +/- 1 kb). A mark that fails this gate cannot have
   its Δ used to predict expression.
"""
from __future__ import annotations
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root

import os

import numpy as np
import pyBigWig
from scipy.stats import rankdata

ROOT = os.path.expanduser(REPO_ROOT + "/CUTRUN")
BW = os.path.join(ROOT, "multimark_tracks", "03_bw")
CMP = os.path.join(ROOT, "multimark_tracks", "04_cmp")
CORE = os.path.join(ROOT, "H3K4me3", "reference", "core.bed")
OUTDIR = os.path.join(ROOT, "analysis_QC_normalization", "multimark_2026-07", "data")

MARKS = ["H3K27ac", "H3K9pr", "ENL"]
NT = {"H3K9pr": {"r1": "shNT-H3K9pr-r1", "r2": "shNT-H3K9pr-r2"},
      "H3K27ac": {"r1": "shNT-H3K27ac-r1", "r2": "shNT-H3K27ac-r2"},
      "ENL": {"r1": "shNT-ENL-r1", "r2": "shNT-ENL-r2"}}
BC = {"H3K9pr": {"r1": "shBCAT1-H3K9pr-r1", "r2": "shBCAT1-H3K9pr-r2_60M"},
      "H3K27ac": {"r1": "shBCAT1-H3K27ac-r1", "r2": "shBCAT1-H3K27ac-r2"},
      "ENL": {"r1": "shBCAT1-ENL-r1", "r2": "shBCAT1-ENL-r2"}}


def sp(a, b):
    ra, rb = rankdata(a), rankdata(b)
    ra = ra - ra.mean(); rb = rb - rb.mean()
    return float(ra @ rb / (np.linalg.norm(ra) * np.linalg.norm(rb)))


def load_core():
    rows = []
    with open(CORE) as fh:
        for line in fh:
            f = line.split()
            rows.append((f[0], int(f[1]), int(f[2]), f[3]))
    return rows


def region_means(path, regions):
    bw = pyBigWig.open(path)
    ch = bw.chroms()
    out = np.full(len(regions), np.nan)
    for i, (c, a, b, _) in enumerate(regions):
        if c not in ch:
            continue
        a2, b2 = max(a, 0), min(b, ch[c])
        if b2 <= a2:
            continue
        try:
            v = bw.stats(c, a2, b2, type="mean")[0]
        except RuntimeError:
            v = None
        if v is not None:
            out[i] = v
    bw.close()
    return out


def wiring_check():
    """Which tier fed each delta track?"""
    print("=== 1. wiring check: which CPM tier fed each delta track ===")
    probe = [("chr1", 1_000_000, 1_400_000), ("chr12", 6_400_000, 6_700_000),
             ("chr20", 30_000_000, 30_400_000)]
    for m in MARKS:
        for rep in ("r1", "r2"):
            d = os.path.join(CMP, f"delta_{m}_{rep}.log2.bw")
            if not os.path.exists(d):
                print(f"  delta_{m}_{rep}: MISSING")
                continue
            dv = np.concatenate([np.array(pyBigWig.open(d).values(*p)) for p in probe])
            best, bestd = None, None
            for tag in ("tierA", "tierB"):
                pn = os.path.join(BW, f"{NT[m][rep]}.{tag}.cpm.bw")
                pb = os.path.join(BW, f"{BC[m][rep]}.{tag}.cpm.bw")
                if not (os.path.exists(pn) and os.path.exists(pb)):
                    continue
                nv = np.concatenate([np.array(pyBigWig.open(pn).values(*p)) for p in probe])
                bv = np.concatenate([np.array(pyBigWig.open(pb).values(*p)) for p in probe])
                exp = np.log2((bv + 1) / (nv + 1))
                ok = np.isfinite(dv) & np.isfinite(exp)
                md = float(np.abs(dv[ok] - exp[ok]).max())
                print(f"    delta_{m}_{rep} vs {tag}: max|diff| = {md:.2e}")
                if bestd is None or md < bestd:
                    best, bestd = tag, md
            verdict = "OK" if bestd is not None and bestd < 1e-4 else "NO EXACT MATCH"
            print(f"  delta_{m}_{rep:2s} -> built from {best}   [{verdict}]")


def concordance_aggregate(regions):
    """The fair version of the test.

    Averaging a per-50 bp log2-ratio track across a promoter is a poor Δ estimator:
    with genome-wide mean coverage ~0.07 CPM, bigwigCompare's `--pseudocount 1`
    dominates every low-coverage bin, so the per-bin ratio is mostly pseudocount and
    the promoter mean is mostly noise. The estimator this project actually uses is
    aggregate-then-ratio: total promoter signal per condition, then one ratio.

    Reported against a signal floor, because concordance in windows with no signal in
    either condition is meaningless and including them only dilutes rho. The floor is
    applied to the *mean of both conditions* so it cannot select for a direction.
    """
    print("\n=== 2b. replicate concordance of promoter Δ — aggregate-then-ratio ===")
    tiers = {}
    for m in MARKS:
        for rep in ("r1", "r2"):
            for tag in ("tierB", "tierA"):        # prefer tierB where it exists
                if os.path.exists(os.path.join(BW, f"{NT[m][rep]}.{tag}.cpm.bw")):
                    tiers[(m, rep)] = tag
                    break
    sets = {}
    for m in MARKS:
        for rep in ("r1", "r2"):
            tag = tiers[(m, rep)]
            a = region_means(os.path.join(BW, f"{NT[m][rep]}.{tag}.cpm.bw"), regions)
            b = region_means(os.path.join(BW, f"{BC[m][rep]}.{tag}.cpm.bw"), regions)
            sets[(m, rep)] = (a, b)
    p1 = os.path.join(ROOT, "H3K4me3", "phase1_tracks", "03_bw")
    for rep in ("r1", "r2"):
        a = region_means(os.path.join(p1, f"shNT-H3K4me3-{rep}.ds.cpm.bw"), regions)
        b = region_means(os.path.join(p1, f"shBCAT1-H3K4me3-{rep}.ds.cpm.bw"), regions)
        sets[("H3K4me3", rep)] = (a, b)

    lines = ["mark,signal_floor_cpm,n_promoters,spearman_r1_vs_r2,pct_same_direction,verdict"]
    for m in ("H3K4me3", "H3K27ac", "H3K9pr", "ENL"):
        a1, b1 = sets[(m, "r1")]
        a2, b2 = sets[(m, "r2")]
        print(f"  {m}")
        for floor in (0.0, 0.1, 0.3, 1.0):
            ok = (np.isfinite(a1) & np.isfinite(b1) & np.isfinite(a2) & np.isfinite(b2)
                  & (((a1 + b1) / 2) > floor) & (((a2 + b2) / 2) > floor))
            if ok.sum() < 200:
                print(f"    floor {floor:>4.1f} CPM: only {ok.sum()} promoters — skipped")
                continue
            eps = 1e-3
            d1 = np.log2((b1[ok] + eps) / (a1[ok] + eps))
            d2 = np.log2((b2[ok] + eps) / (a2[ok] + eps))
            r = sp(d1, d2)
            same = float(np.mean(np.sign(d1) == np.sign(d2)) * 100)
            verdict = "PASS" if (r > 0.3 and same > 60) else "FAIL"
            print(f"    floor {floor:>4.1f} CPM: n={ok.sum():6d}  rho = {r:+.3f}  "
                  f"same-direction {same:5.1f}%  -> {verdict}")
            lines.append(f"{m},{floor},{ok.sum()},{r:.4f},{same:.1f},{verdict}")
    p = os.path.join(OUTDIR, "delta_replicate_concordance_promoters_aggregate.csv")
    with open(p, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("[ok]", p)


def concordance(regions):
    print("\n=== 2a. replicate concordance of promoter Δ (mean of log2-ratio track) ===")
    print(f"    {len(regions):,} promoter windows")
    lines = ["mark,n_promoters_both_reps,spearman_r1_vs_r2,pct_same_direction,"
             "gate_rho_gt_0.3,gate_dir_gt_60pct,verdict"]
    # H3K4me3 deltas live with the phase 1 tracks
    extra = {"H3K4me3": [os.path.join(ROOT, "H3K4me3", "phase1_tracks", "04_cmp",
                                      f"delta_H3K4me3_{r}.log2.bw") for r in ("r1", "r2")]}
    todo = {m: [os.path.join(CMP, f"delta_{m}_{r}.log2.bw") for r in ("r1", "r2")]
            for m in MARKS}
    todo.update(extra)
    for m in ("H3K4me3", "H3K27ac", "H3K9pr", "ENL"):
        p1, p2 = todo[m]
        if not (os.path.exists(p1) and os.path.exists(p2)):
            print(f"  {m:8s} missing delta track(s) — skipped")
            continue
        a, b = region_means(p1, regions), region_means(p2, regions)
        ok = np.isfinite(a) & np.isfinite(b) & ((a != 0) | (b != 0))
        r = sp(a[ok], b[ok])
        same = float(np.mean(np.sign(a[ok]) == np.sign(b[ok])) * 100)
        g1, g2 = r > 0.3, same > 60
        verdict = "PASS" if (g1 and g2) else "FAIL"
        print(f"  {m:8s} n={ok.sum():6d}  rho(r1,r2) = {r:+.3f}   "
              f"same-direction {same:5.1f}%   -> {verdict}")
        lines.append(f"{m},{ok.sum()},{r:.4f},{same:.1f},{g1},{g2},{verdict}")
    p = os.path.join(OUTDIR, "delta_replicate_concordance_promoters.csv")
    with open(p, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("[ok]", p)


if __name__ == "__main__":
    wiring_check()
    regions = load_core()
    concordance(regions)
    concordance_aggregate(regions)
