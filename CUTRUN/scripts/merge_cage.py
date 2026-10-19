#!/usr/bin/env python
"""
Merge stranded CAGE CPM bigWigs into one unstranded initiation profile, and bin
it onto the two grids used for the correlation analysis.

THE MINUS STRAND IS STORED AS NEGATIVE CPM (CAGE_bigwig_tracks/README.md, and
track_summary.tsv shows minus spanning -5335..0 against plus 0..8861). So the
merge is plus + |minus|. A plain `bigwigCompare --operation add` would SUBTRACT
the minus strand from the plus and silently cancel signal at every bidirectional
promoter -- exactly the regions of interest.

CAGE is 1 bp point signal covering ~2.65 Mb, i.e. 0.088% of the genome, so it is
read via bw.intervals() (a few million non-zero runs) rather than by
materialising per-base arrays for 3 Gb.

Outputs
  <out>.merged.bw            plus + |minus|, unstranded, CPM
  cage_10kb.tsv.gz           summed CPM per 10 kb bin, on the project grid
  cage_promoters.tsv.gz      summed CPM per core.bed window (TSS +/- 1 kb)

Usage: merge_cage.py [--sample S.K27M]
"""
from __future__ import annotations

import argparse
import gzip
import os

import numpy as np
import pyBigWig

HOME = os.path.expanduser("~")
CAGE = os.path.join(HOME, "Desktop/BCAT1/BCAT1/CAGE_bigwig_tracks")
ROOT = os.path.join(HOME, "Desktop/BCAT1/BCAT1/CUTRUN")
QC = os.path.join(ROOT, "analysis_QC_normalization", "multimark_2026-07")
GRID = os.path.join(ROOT, "analysis_QC_normalization", "data", "counts_10kb_matrix.npz")
CORE = os.path.join(ROOT, "H3K4me3", "reference", "core.bed")
BIN = 10_000
BIG = 300_000_000
MAIN = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
CIDX = {c: i for i, c in enumerate(MAIN)}


def load_stranded(sample):
    """Return (chrom_idx, start, end, value) for plus + |minus| combined."""
    parts = []
    for strand in ("plus", "minus"):
        p = os.path.join(CAGE, f"{sample}.{strand}.cpm.bw")
        bw = pyBigWig.open(p)
        tot = 0.0
        for c in bw.chroms():
            if c not in CIDX:
                continue
            iv = bw.intervals(c)
            if not iv:
                continue
            a = np.array(iv, dtype=np.float64)
            v = np.abs(a[:, 2])                      # minus is stored negative
            parts.append((np.full(len(a), CIDX[c], dtype=np.int64),
                          a[:, 0].astype(np.int64), a[:, 1].astype(np.int64), v))
            tot += float((v * (a[:, 1] - a[:, 0])).sum())
        print(f"[load] {sample}.{strand}: total |CPM|*bp = {tot:,.0f}")
        bw.close()
    ci = np.concatenate([p[0] for p in parts])
    st = np.concatenate([p[1] for p in parts])
    en = np.concatenate([p[2] for p in parts])
    va = np.concatenate([p[3] for p in parts])
    print(f"[load] {len(ci):,} non-zero intervals, "
          f"{int((en-st).sum()):,} bp covered")
    return ci, st, en, va


def write_merged_bw(ci, st, en, va, out):
    """plus + |minus| on one track. Overlapping strands are summed."""
    src = pyBigWig.open(os.path.join(CAGE, "S.K27M.plus.cpm.bw"))
    chroms = {c: l for c, l in src.chroms().items() if c in CIDX}
    src.close()
    bw = pyBigWig.open(out, "w")
    bw.addHeader([(c, chroms[c]) for c in MAIN if c in chroms])
    for c in MAIN:
        if c not in chroms:
            continue
        m = ci == CIDX[c]
        if not m.any():
            continue
        s, e, v = st[m], en[m], va[m]
        o = np.argsort(s, kind="stable")
        s, e, v = s[o], e[o], v[o]
        # both strands can mark the same base; collapse duplicates by summing
        key = s * 4 + (e - s)
        uk, inv = np.unique(key, return_inverse=True)
        vs = np.bincount(inv, weights=v)
        us = (uk // 4).astype(np.int64)
        ue = us + (uk % 4).astype(np.int64)
        ok = ue > us
        bw.addEntries([c] * int(ok.sum()), us[ok].tolist(),
                      ends=ue[ok].tolist(), values=vs[ok].astype(np.float64).tolist())
    bw.close()
    print("[bw]", out)


def bin_10kb(ci, st, va, out):
    d = np.load(GRID, allow_pickle=True)
    chroms = np.asarray(d["chroms"], dtype="U6")
    starts = np.asarray(d["starts"], dtype=np.int64)
    idx = {(c, s): i for i, (c, s) in enumerate(zip(chroms, starts))}
    acc = np.zeros(len(chroms))
    hit = 0
    for i, c in enumerate(MAIN):
        m = ci == i
        if not m.any():
            continue
        b = (st[m] // BIN) * BIN
        v = va[m]
        for bb, vv in zip(b, v):
            k = idx.get((c, int(bb)))
            if k is not None:
                acc[k] += vv
                hit += 1
    with gzip.open(out, "wt") as fh:
        for (c, s), a in zip(zip(chroms, starts), acc):
            if a > 0:
                fh.write(f"{c}\t{s}\t{a:.6f}\n")
    nz = int((acc > 0).sum())
    print(f"[10kb] {nz:,}/{len(acc):,} bins non-zero ({100*nz/len(acc):.1f}%) -> {out}")
    return acc


def bin_regions(ci, st, va, bed, out):
    rows = []
    with open(bed) as fh:
        for i, line in enumerate(fh):
            f = line.split()
            if f[0] in CIDX:
                rows.append((CIDX[f[0]], int(f[1]), int(f[2]), i))
    n_all = i + 1
    rows.sort(key=lambda r: (r[0], r[1]))
    ks = np.array([r[0] * BIG + r[1] for r in rows], dtype=np.int64)
    ke = np.array([r[0] * BIG + r[2] for r in rows], dtype=np.int64)
    rid = np.array([r[3] for r in rows], dtype=np.int64)
    width = int((ke - ks).max())
    lo = np.searchsorted(ks, ks - width, side="right")
    K = int((np.arange(len(ks)) - lo + 1).max())

    key = ci * BIG + st
    acc = np.zeros(n_all)
    j = np.searchsorted(ks, key, side="right") - 1
    for k in range(K):
        jj = j - k
        m = (jj >= 0) & (key < ke[np.maximum(jj, 0)]) & (ks[np.maximum(jj, 0)] <= key)
        if m.any():
            np.add.at(acc, rid[jj[m]], va[m])
    with gzip.open(out, "wt") as fh:
        for i, a in enumerate(acc):
            if a > 0:
                fh.write(f"{i+1}\t{a:.6f}\n")
    nz = int((acc > 0).sum())
    print(f"[regions] {nz:,}/{n_all:,} windows with CAGE ({100*nz/n_all:.1f}%), "
          f"walk-back K={K} -> {out}")
    return acc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", default="S.K27M")
    a = ap.parse_args()
    outdir = os.path.join(QC, "cage")
    os.makedirs(outdir, exist_ok=True)

    ci, st, en, va = load_stranded(a.sample)
    write_merged_bw(ci, st, en, va, os.path.join(outdir, f"{a.sample}.merged.cpm.bw"))
    bin_10kb(ci, st, va, os.path.join(outdir, f"{a.sample}.cage_10kb.tsv.gz"))
    bin_regions(ci, st, va, CORE, os.path.join(outdir, f"{a.sample}.cage_promoters.tsv.gz"))
    # same disjoint promoter/exonic/intronic/intergenic partition used for the
    # mark category panel, so CAGE can be shown on the same axes
    part = os.path.join(QC, "annotation", "partition.bed")
    if os.path.exists(part):
        acc = bin_regions(ci, st, va, part,
                          os.path.join(outdir, f"{a.sample}.cage_partition.tsv.gz"))
        cats = {}
        with open(part) as fh:
            for i, line in enumerate(fh):
                cats.setdefault(int(line.split()[3]), []).append(i)
        tot = acc.sum()
        names = {1: "promoter", 2: "exonic", 3: "intronic", 4: "intergenic"}
        with open(os.path.join(outdir, f"{a.sample}.cage_categories.tsv"), "w") as fh:
            fh.write("category\tcpm\tpct\n")
            for k in sorted(cats):
                v = acc[cats[k]].sum()
                fh.write(f"{names[k]}\t{v:.4f}\t{100*v/tot:.4f}\n")
                print(f"[cage-cat] {names[k]:11s} {100*v/tot:6.2f}% of CAGE signal")


if __name__ == "__main__":
    main()
