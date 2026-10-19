#!/usr/bin/env python
"""
Per-region fragment counts for ANY BED, one pass per library.

Generalises count_promoter_and_category.py (which was hard-wired to the 2 kb
promoter windows) so the same estimator serves gene bodies, peaks, enhancers or
any other region set.

Fragment definition is identical to every other count in this project:
  -f 66 -F 1804, MAPQ >= 30, chr1-22 + chrX, GenBank->UCSC rename for H3K27ac,
  placed by MIDPOINT = min(POS, PNEXT) - 1 + |TLEN| // 2.

A fragment counts toward EVERY region containing its midpoint. That is correct
for region quantification (overlapping promoters legitimately share a fragment)
and is why these totals can exceed the library total.

--------------------------------------------------------------------------
The one non-obvious bit: how far to walk back when regions overlap
--------------------------------------------------------------------------
Regions are sorted by start and located with `searchsorted`, which gives the LAST
region whose start <= mid. Other regions containing mid sit earlier in the array,
so the search has to walk backwards -- but by how many steps?

"Max overlap depth" is NOT a valid bound once region widths differ. Counter-example:
A=[0,100), B=[10,20), C=[30,40), D=[50,60), mid=90. Max depth is 2, searchsorted
lands on D, and only A actually contains 90 -- four steps back. Walking `depth`
steps would silently miss it.

The valid bound is width-based: any region containing mid must start after
mid - max_region_width. So K is computed as the largest number of region starts
falling inside any window of that width, which is exact for any BED. For the
equal-width promoter case this reduces to the old depth bound, which is why that
version was correct.

Usage: count_regions.py <label> <naming: ucsc|gb> <bam> <bed> <outdir>
"""
from __future__ import annotations

import os
import subprocess
import sys

import numpy as np
import pandas as pd

BIG = 300_000_000
CHUNK = 4_000_000

GB2UCSC = {f"CM0006{n:02d}.2": f"chr{c}" for n, c in zip(range(63, 85), range(1, 23))}
GB2UCSC["CM000685.2"] = "chrX"
MAIN = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
CIDX = {c: i for i, c in enumerate(MAIN)}


def load_regions(bed):
    df = pd.read_csv(bed, sep="\t", header=None, usecols=[0, 1, 2],
                     names=["chrom", "start", "end"], comment="#")
    df["rid"] = np.arange(len(df), dtype=np.int64)      # row order in the BED
    n_all = len(df)
    df = df[df.chrom.isin(CIDX) & (df.end > df.start)]
    if not len(df):
        sys.exit(f"[fatal] no usable regions on chr1-22,X in {bed}")
    ci = df.chrom.map(CIDX).to_numpy(np.int64)
    ks = ci * BIG + df.start.to_numpy(np.int64)
    ke = ci * BIG + df.end.to_numpy(np.int64)
    o = np.argsort(ks, kind="stable")
    ks, ke, rid = ks[o], ke[o], df.rid.to_numpy(np.int64)[o]

    width = int((df.end - df.start).max())
    # largest number of starts inside any window of `width` -> exact walk-back bound
    lo = np.searchsorted(ks, ks - width, side="right")
    K = int((np.arange(len(ks)) - lo + 1).max())
    print(f"[regions] {len(df):,}/{n_all:,} usable  max width {width:,} bp  "
          f"walk-back K={K}", flush=True)
    return ks, ke, rid, n_all, K


def main():
    label, naming, bam, bed, outdir = sys.argv[1:6]
    os.makedirs(outdir, exist_ok=True)
    out = os.path.join(outdir, f"{label}.regions.tsv.gz")
    if os.path.exists(out):
        print(f"[skip] {label}")
        return

    ks, ke, rid, n_all, K = load_regions(bed)
    counts = np.zeros(n_all, dtype=np.int64)

    regions = list(GB2UCSC) if naming == "gb" else MAIN
    proc = subprocess.Popen(
        ["samtools", "view", "-@", "3", "-q", "30", "-f", "66", "-F", "1804", bam]
        + regions, stdout=subprocess.PIPE, bufsize=1 << 22)
    print(f"[count] {label} ({naming})", flush=True)

    total = 0
    reader = pd.read_csv(proc.stdout, sep="\t", header=None, usecols=[2, 3, 7, 8],
                         names=["rname", "pos", "pnext", "tlen"],
                         dtype={2: str, 3: np.int64, 7: np.int64, 8: np.int64},
                         chunksize=CHUNK, quoting=3, engine="c")
    for chunk in reader:
        chrom = chunk.rname.map(GB2UCSC) if naming == "gb" else chunk.rname
        ci = chrom.map(CIDX)
        ok = ci.notna().to_numpy()
        if not ok.any():
            continue
        ci = ci.to_numpy()[ok].astype(np.int64)
        pos = chunk.pos.to_numpy()[ok]
        pnext = chunk.pnext.to_numpy()[ok]
        tlen = np.abs(chunk.tlen.to_numpy()[ok])
        mid = np.minimum(pos, pnext) - 1 + (tlen // 2)
        key = ci * BIG + mid
        total += key.size

        j = np.searchsorted(ks, key, side="right") - 1
        for k in range(K):
            jj = j - k
            valid = jj >= 0
            if not valid.any():
                break
            jc = np.maximum(jj, 0)
            m = valid & (key < ke[jc]) & (ks[jc] <= key)
            if m.any():
                counts += np.bincount(rid[jj[m]], minlength=n_all)

    proc.stdout.close()
    if proc.wait() != 0:
        sys.exit(f"[fatal] samtools failed for {label}")

    nz = np.nonzero(counts)[0]
    pd.DataFrame({"region_row": nz + 1, "fragments": counts[nz]}).to_csv(
        out + ".partial", sep="\t", index=False, compression="gzip")
    os.replace(out + ".partial", out)
    print(f"[done] {label}  fragments={total:,}  in-region sum={counts.sum():,}  "
          f"regions hit={len(nz):,}", flush=True)


if __name__ == "__main__":
    main()
