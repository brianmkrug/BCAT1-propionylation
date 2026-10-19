#!/usr/bin/env python
"""
regions x marks DIRECT-RATIO delta matrix for any BED.

This is the canonical quantitative estimator for cross-condition change
(see multimark_2026-07/PLAN_direct_ratio_migration.md):

    replicates POOLED at full depth (counts summed)
    -> CPM against the GENOME-WIDE pooled depth
       (not the in-region total -- otherwise a mark that redistributes INTO the
        region set would normalise its own effect away)
    -> spike-in factor on the shBCAT1 side
    -> delta = log2( shBCAT1 / shNT )

No IgG term: copy number is identical in the two conditions and cancels exactly
in the ratio. IgG is still counted and reported, purely as a background reference
column for judging whether a region has signal at all.

Counts come from count_regions.py, which places each fragment once by midpoint,
so this matrix is immune to the fragment-length weighting that affects
coverage-based (bigWig) quantification.

Usage:
  build_region_delta_matrix.py --bed <file.bed> --name <label> [--jobs 3]
"""
from __future__ import annotations
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root

import argparse
import os
import subprocess
import sys
import time

import numpy as np
import pandas as pd

ROOT = os.path.expanduser(REPO_ROOT + "/CUTRUN")
QC = os.path.join(ROOT, "analysis_QC_normalization", "multimark_2026-07")
PROMCAT = os.path.join(QC, "promcat")
SCRIPTS = os.path.join(ROOT, "scripts")

SPIKEIN = {"ENL": 0.80, "H3K9pr": 0.62, "H3K27ac": 1.00, "H3K4me3": 1.00, "IgG": 1.00}
ASSUMED = {"H3K4me3"}      # no measured factor supplied; carried at 1.00
MARKS = ["H3K4me3", "H3K27ac", "H3K9pr", "ENL", "IgG"]


def manifest():
    return pd.read_csv(os.path.join(QC, "sample_manifest.tsv"), sep="\t")


def genome_depth(label):
    """Library total on chr1-22,X at MAPQ30 -- the CPM denominator."""
    p = os.path.join(PROMCAT, f"{label}.category.tsv")
    if not os.path.exists(p):
        sys.exit(f"[fatal] missing {p}; run count_promoter_and_category.py first")
    return int(pd.read_csv(p, sep="\t", index_col=0).fragments["total"])


def run_counts(man, bed, outdir, jobs):
    todo = []
    for _, r in man.iterrows():
        if not os.path.exists(os.path.join(outdir, f"{r.label}.regions.tsv.gz")):
            todo.append((r.label, r.naming, os.path.join(ROOT, "hg38", r.source_bam)))
    if not todo:
        print("[counts] all libraries already counted")
        return
    print(f"[counts] {len(todo)} libraries to count, {jobs} at a time")
    running, allp = [], []
    for label, naming, bam in todo:
        while len(running) >= jobs:
            time.sleep(0.5)
            running = [p for p in running if p.poll() is None]
        p = subprocess.Popen(
            [sys.executable, os.path.join(SCRIPTS, "count_regions.py"),
             label, naming, bam, bed, outdir])
        running.append(p)
        allp.append(p)
    for p in allp:
        p.wait()
    bad = [p.args[2] for p in allp if p.returncode != 0]
    if bad:
        sys.exit(f"[fatal] counting failed for: {', '.join(bad)}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bed", required=True)
    ap.add_argument("--name", required=True)
    ap.add_argument("--jobs", type=int, default=3)
    a = ap.parse_args()

    man = manifest()
    outdir = os.path.join(QC, "regions", a.name, "counts")
    os.makedirs(outdir, exist_ok=True)
    run_counts(man, os.path.abspath(a.bed), outdir, a.jobs)

    bed = pd.read_csv(a.bed, sep="\t", header=None, comment="#")
    bed = bed.rename(columns={0: "chrom", 1: "start", 2: "end"})
    n = len(bed)
    out = bed.iloc[:, :min(4, bed.shape[1])].copy()
    if bed.shape[1] > 3:
        out = out.rename(columns={3: "name"})

    pooled, depth = {}, {}
    for m in MARKS:
        for cond in ("shNT", "shBCAT1"):
            labs = man[(man["mark"] == m) & (man.condition == cond)].label.tolist()
            v = np.zeros(n, dtype=np.int64)
            d = 0
            for lab in labs:
                t = pd.read_csv(os.path.join(outdir, f"{lab}.regions.tsv.gz"), sep="\t")
                v[t.region_row.to_numpy() - 1] += t.fragments.to_numpy()
                d += genome_depth(lab)
            pooled[(m, cond)], depth[(m, cond)] = v, d
            print(f"[pool] {cond}-{m:8s} {'+'.join(labs):45s} "
                  f"in-region {v.sum():>12,}  genome depth {d:>12,}")

    eps = 1e-2
    for m in MARKS:
        f = SPIKEIN[m]
        nt = pooled[(m, "shNT")] / depth[(m, "shNT")] * 1e6
        bc = pooled[(m, "shBCAT1")] / depth[(m, "shBCAT1")] * 1e6 * f
        out[f"{m}_shNT_cpm"] = nt
        out[f"{m}_shBCAT1_cpm"] = bc
        if m != "IgG":
            out[f"{m}_delta"] = np.log2((bc + eps) / (nt + eps))
            out[f"{m}_delta_nospike"] = np.log2((bc / f + eps) / (nt + eps))

    dest = os.path.join(QC, "data", f"{a.name}_delta_matrix.csv")
    out.to_csv(dest, index=False)
    print(f"\n[ok] {dest}  ({n:,} regions x {len(MARKS)} marks)")
    print("[note] spike-in: " + ", ".join(
        f"{m} x{SPIKEIN[m]:.2f}" + ("*" if m in ASSUMED else "") for m in MARKS)
        + "   (* = assumed, not measured)")
    print("[note] delta = log2(shBCAT1/shNT), direct ratio, no IgG term; "
          "IgG columns are a background reference only")
    for m in MARKS:
        if m == "IgG":
            continue
        d = out[f"{m}_delta"]
        sig = (out[f"{m}_shNT_cpm"] + out[f"{m}_shBCAT1_cpm"]) > 0
        print(f"  {m:8s} median delta {d[sig].median():+.3f}  "
              f"({int(sig.sum()):,} regions with signal)")


if __name__ == "__main__":
    main()
