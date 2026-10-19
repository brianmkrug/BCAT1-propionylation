#!/usr/bin/env python
"""
Per-sample normalization metrics for the multi-mark panel, reproducing the
columns of analysis_QC_normalization/data/normalization_metrics.csv so the new
H3K27ac libraries can be read against the 2026-07-08 numbers:

  spearman_raw_vs_IgG          how much of the track is copy-number background
  spearman_log2ratio_vs_IgG    residual after log2(IP/matched IgG), CPM-scaled
  pct_CNV_removed              1 - |residual| / |raw|

Plus, for this pass:
  replicate concordance per mark/condition (r1 vs r2), and the within-batch vs
  within-condition split for H3K27ac, which is the thing that decides whether
  the two H3K27ac replicates can be pooled.
"""
from __future__ import annotations
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root

import os

import numpy as np
from scipy.stats import rankdata

ROOT = os.path.expanduser(REPO_ROOT + "/CUTRUN")
OUT = os.path.join(ROOT, "analysis_QC_normalization", "multimark_2026-07")

MARKS = {"H3K4me3": "H3K4me3", "H3K27ac": "H3K27ac", "H3K9pr": "H3K9pr",
         "ENL": "ENL", "IgG": "IgG"}


def mark_of(l):
    for m in MARKS:
        if f"-{m}-" in l:
            return m
    raise ValueError(l)


def cond_of(l):
    return "shBCAT1" if l.startswith("shBCAT1") else "shNT"


def sp(a, b):
    ra, rb = rankdata(a), rankdata(b)
    ra = ra - ra.mean(); rb = rb - rb.mean()
    return float(ra @ rb / (np.linalg.norm(ra) * np.linalg.norm(rb)))


def main():
    d = np.load(os.path.join(OUT, "data", "multimark_counts_10kb.npz"), allow_pickle=True)
    labels = [str(x) for x in d["labels"]]
    keep = d["blacklist_keep"]
    C = d["M"][keep]
    cpm = C / C.sum(axis=0, keepdims=True) * 1e6
    col = {l: i for i, l in enumerate(labels)}

    igg = {c: cpm[:, col[f"{c}-IgG-r1"]] for c in ("shNT", "shBCAT1")}

    rows = []
    for l in labels:
        m, c = mark_of(l), cond_of(l)
        if m == "IgG":
            rows.append((l, m, c, "", "", ""))
            continue
        ip, bg = cpm[:, col[l]], igg[c]
        raw = sp(ip, bg)
        ratio = np.log2((ip + 1) / (bg + 1))
        res = sp(ratio, bg)
        pct = 100 * (1 - abs(res) / abs(raw)) if raw else float("nan")
        rows.append((l, m, c, f"{raw:.3f}", f"{res:.3f}", f"{pct:.0f}"))

    p = os.path.join(OUT, "data", "normalization_metrics.csv")
    with open(p, "w") as fh:
        fh.write("sample,mark,condition,spearman_raw_vs_IgG,"
                 "spearman_log2ratio_vs_IgG,pct_CNV_removed\n")
        for r in rows:
            fh.write(",".join(r) + "\n")
    print("[ok]", p)
    for r in rows:
        print(f"  {r[0]:24s} {r[1]:8s} raw={r[3]:>6s}  residual={r[4]:>6s}  removed={r[5]:>3s}%")

    # ---- replicate / batch concordance on log2 CPM ----------------------------
    X = np.log2(cpm + 1)
    print("\n[replicate concordance] Spearman on log2 CPM, 10 kb bins")
    lines = ["comparison,mark,type,spearman"]
    for m in ("H3K4me3", "H3K27ac", "H3K9pr", "ENL"):
        labs = [l for l in labels if mark_of(l) == m]
        for c in ("shNT", "shBCAT1"):
            pair = [l for l in labs if cond_of(l) == c]
            if len(pair) == 2:
                r = sp(X[:, col[pair[0]]], X[:, col[pair[1]]])
                print(f"  {m:8s} {c:8s} r1 vs r2            r = {r:.3f}")
                lines.append(f"{pair[0]} vs {pair[1]},{m},same-condition-across-rep,{r:.4f}")
        for rep in ("r1", "r2"):
            pair = [l for l in labs if l.split("-")[-1].startswith(rep)]
            if len(pair) == 2:
                r = sp(X[:, col[pair[0]]], X[:, col[pair[1]]])
                print(f"  {m:8s} {rep:8s} shNT vs shBCAT1     r = {r:.3f}")
                lines.append(f"{pair[0]} vs {pair[1]},{m},same-rep-across-condition,{r:.4f}")
    # ---- the same comparisons at MATCHED depth -------------------------------
    # A replicate pair with a 5-10x depth gap loses correlation to shot noise
    # alone, so a low raw value is not by itself evidence of a batch effect.
    # Thin the deeper library of each pair down to the shallower one and repeat:
    # what survives is real disagreement.
    rng = np.random.default_rng(42)
    tot = C.sum(axis=0)

    def thinned_pair(a, b):
        ia, ib = col[a], col[b]
        t = min(tot[ia], tot[ib])
        out = []
        for i in (ia, ib):
            p = t / tot[i]
            v = C[:, i] if p >= 0.999 else rng.binomial(C[:, i].astype(np.int64), p)
            out.append(np.log2(v / v.sum() * 1e6 + 1))
        return sp(out[0], out[1]), t

    print("\n[replicate concordance, depth-matched] deeper library thinned to the shallower")
    lines2 = ["comparison,mark,type,spearman_depth_matched,fragments_each"]
    for m in ("H3K4me3", "H3K27ac", "H3K9pr", "ENL"):
        labs = [l for l in labels if mark_of(l) == m]
        for c in ("shNT", "shBCAT1"):
            pair = [l for l in labs if cond_of(l) == c]
            if len(pair) == 2:
                r, t = thinned_pair(*pair)
                print(f"  {m:8s} {c:8s} r1 vs r2   r = {r:.3f}   ({t/1e6:.1f}M each)")
                lines2.append(f"{pair[0]} vs {pair[1]},{m},same-condition-across-rep,{r:.4f},{t:.0f}")
        for rep in ("r1", "r2"):
            pair = [l for l in labs if l.split("-")[-1].startswith(rep)]
            if len(pair) == 2:
                r, t = thinned_pair(*pair)
                print(f"  {m:8s} {rep:8s} NT vs BCAT1 r = {r:.3f}   ({t/1e6:.1f}M each)")
                lines2.append(f"{pair[0]} vs {pair[1]},{m},same-rep-across-condition,{r:.4f},{t:.0f}")

    p2 = os.path.join(OUT, "data", "replicate_concordance.csv")
    with open(p2, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("[ok]", p2)
    p3 = os.path.join(OUT, "data", "replicate_concordance_depthmatched.csv")
    with open(p3, "w") as fh:
        fh.write("\n".join(lines2) + "\n")
    print("[ok]", p3)


if __name__ == "__main__":
    main()
