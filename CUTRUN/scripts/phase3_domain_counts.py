#!/usr/bin/env python
"""
PHASE 3 (prep) — per-LIBRARY counts on each mark's union domains, plus the
sample table edgeR needs.

Unlike build_region_delta_matrix.py, this keeps the four libraries of a mark
SEPARATE. Pooling throws away exactly the replicate-level variation that the
differential model needs in order to say anything about uncertainty.

Writes, per mark:
  domains/<mark>.counts.tsv    domains x 4 libraries (raw fragment counts)
  domains/<mark>.samples.tsv   label, condition, replicate, genome_depth,
                               spikein, edger_lib_size

`edger_lib_size = genome_depth / spikein` is how the spike-in enters the model.
edgeR normalises counts by lib.size, so dividing the library size by the factor
makes the modelled abundance equal count/depth*spikein -- the spike-in-scaled
quantity -- while the DISPERSION is still estimated from raw counts, which is
the whole point of using an offset instead of multiplying afterwards.
"""
from __future__ import annotations
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root

import os
import subprocess
import sys
import time

import numpy as np
import pandas as pd

ROOT = os.path.expanduser(REPO_ROOT + "/CUTRUN")
QC = os.path.join(ROOT, "analysis_QC_normalization", "multimark_2026-07")
DOM = os.path.join(QC, "domains")
PROMCAT = os.path.join(QC, "promcat")
SCRIPTS = os.path.join(ROOT, "scripts")

SPIKEIN = {"H3K9pr": 0.62, "H3K27ac": 1.00, "H3K4me3": 1.00}
# ENL dropped: at full combined depth against pooled IgG, SEACR called 63 regions
# genome-wide with an empirical FDR of 0.433 -- it has no signal separable from
# background, consistent with 1.64x promoter enrichment (IgG 1.18x), the closest
# genome-wide correlation to IgG of any mark, and a failed replicate-concordance
# gate. Testing 54 domains would manufacture a result. ENL stays on the fixed
# region sets (promoters, gene bodies) as a bounded negative.
MARKS = ["H3K4me3", "H3K27ac", "H3K9pr"]
JOBS = 3


def genome_depth(label):
    return int(pd.read_csv(os.path.join(PROMCAT, f"{label}.category.tsv"),
                           sep="\t", index_col=0).fragments["total"])


def main():
    man = pd.read_csv(os.path.join(QC, "sample_manifest.tsv"), sep="\t")
    for mark in MARKS:
        bed = os.path.join(DOM, f"{mark}.union.bed")
        if not os.path.exists(bed):
            print(f"[skip] {mark}: no union domains yet")
            continue
        outdir = os.path.join(DOM, "counts", mark)
        os.makedirs(outdir, exist_ok=True)
        libs = man[man["mark"] == mark]

        running, allp = [], []
        for _, r in libs.iterrows():
            if os.path.exists(os.path.join(outdir, f"{r.label}.regions.tsv.gz")):
                continue
            while len(running) >= JOBS:
                time.sleep(0.5)
                running = [p for p in running if p.poll() is None]
            p = subprocess.Popen(
                [sys.executable, os.path.join(SCRIPTS, "count_regions.py"),
                 r.label, r.naming, os.path.join(ROOT, "hg38", r.source_bam),
                 bed, outdir])
            running.append(p); allp.append(p)
        for p in allp:
            p.wait()
        if any(p.returncode != 0 for p in allp):
            sys.exit(f"[fatal] counting failed for {mark}")

        dom = pd.read_csv(bed, sep="\t", header=None,
                          names=["chrom", "start", "end", "name"])
        n = len(dom)
        mat = dom[["chrom", "start", "end", "name"]].copy()
        rows = []
        for _, r in libs.iterrows():
            v = np.zeros(n, dtype=np.int64)
            t = pd.read_csv(os.path.join(outdir, f"{r.label}.regions.tsv.gz"), sep="\t")
            v[t.region_row.to_numpy() - 1] = t.fragments.to_numpy()
            mat[r.label] = v
            d = genome_depth(r.label)
            f = SPIKEIN[mark] if r.condition == "shBCAT1" else 1.0
            rows.append({"label": r.label, "condition": r.condition,
                         "replicate": r.rep, "genome_depth": d, "spikein": f,
                         "edger_lib_size": d / f,
                         "in_domain": int(v.sum()),
                         "frip_pct": round(100 * v.sum() / d, 2)})
        mat.to_csv(os.path.join(DOM, f"{mark}.counts.tsv"), sep="\t", index=False)
        st = pd.DataFrame(rows)
        st.to_csv(os.path.join(DOM, f"{mark}.samples.tsv"), sep="\t", index=False)
        print(f"[{mark}] {n:,} domains")
        print(st[["label", "condition", "replicate", "genome_depth",
                  "spikein", "frip_pct"]].to_string(index=False))
        print()


if __name__ == "__main__":
    main()
