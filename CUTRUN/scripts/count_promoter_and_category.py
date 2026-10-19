#!/usr/bin/env python
"""
Per-library fragment counts in (a) the 36,689 core.bed promoter windows and
(b) the 4 disjoint genomic categories, in ONE pass over each BAM.

Fragment definition matches every other count in this project: -f 66 -F 1804,
MAPQ >= 30, chr1-22 + chrX, GenBank->UCSC rename for the H3K27ac libraries.

Fragments are placed by MIDPOINT. At 10 kb the read1 5' end was fine; over a 2 kb
promoter window a ~130 bp offset is not. The midpoint is exact integer arithmetic
with no CIGAR/SEQ parsing, because for a properly-paired pair the fragment spans
[min(POS,PNEXT), min(POS,PNEXT)+|TLEN|):

    mid = min(POS, PNEXT) - 1 + |TLEN| // 2

Categories are disjoint, so each fragment lands in exactly one and the four counts
sum to the library total -- asserted, not assumed. Promoter windows overlap (up to
7 deep), so a fragment may count toward several promoters; that is correct for
per-promoter quantification, and it is why promoter totals exceed the promoter
category count.

Chromosome and position are packed into a single int64 key (chrom_idx * 3e8 + pos)
so lookup is one global np.searchsorted rather than a per-chromosome loop.

Usage: count_promoter_and_category.py <label> <naming: ucsc|gb> <bam> <outdir>
"""
from __future__ import annotations
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root

import os
import subprocess
import sys

import numpy as np
import pandas as pd

ROOT = os.path.expanduser(REPO_ROOT + "/CUTRUN")
ANN = os.path.join(ROOT, "analysis_QC_normalization", "multimark_2026-07", "annotation")
CORE = os.path.join(ROOT, "H3K4me3", "reference", "core.bed")
BIG = 300_000_000                     # > longest chromosome (248,956,422)
MAXDEPTH = 8                          # max overlapping promoters is 7; one to spare
CHUNK = 4_000_000

GB2UCSC = {f"CM0006{n:02d}.2": f"chr{c}" for n, c in zip(
    range(63, 85), list(range(1, 23)))}
GB2UCSC["CM000685.2"] = "chrX"
MAIN = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
CIDX = {c: i for i, c in enumerate(MAIN)}
CATS = ["promoter", "exonic", "intronic", "intergenic"]


def load_partition():
    df = pd.read_csv(os.path.join(ANN, "partition.bed"), sep="\t", header=None,
                     names=["chrom", "start", "end", "cat"])
    df = df[df.chrom.isin(CIDX)]
    key = df.chrom.map(CIDX).to_numpy(np.int64) * BIG + df.start.to_numpy(np.int64)
    o = np.argsort(key, kind="stable")
    return (key[o],
            df.chrom.map(CIDX).to_numpy(np.int64)[o] * BIG + df.end.to_numpy(np.int64)[o],
            df.cat.to_numpy(np.int8)[o])


def load_promoters():
    df = pd.read_csv(CORE, sep="\t", header=None, usecols=[0, 1, 2],
                     names=["chrom", "start", "end"])
    df["pid"] = np.arange(len(df), dtype=np.int64)     # row order == core.bed order
    df = df[df.chrom.isin(CIDX)]
    ci = df.chrom.map(CIDX).to_numpy(np.int64)
    ks = ci * BIG + df.start.to_numpy(np.int64)
    ke = ci * BIG + df.end.to_numpy(np.int64)
    o = np.argsort(ks, kind="stable")
    return ks[o], ke[o], df.pid.to_numpy(np.int64)[o], len(df)


def main():
    label, naming, bam, outdir = sys.argv[1:5]
    os.makedirs(outdir, exist_ok=True)
    out_p = os.path.join(outdir, f"{label}.promoter.tsv.gz")
    out_c = os.path.join(outdir, f"{label}.category.tsv")
    if os.path.exists(out_p) and os.path.exists(out_c):
        print(f"[skip] {label}")
        return

    pk_s, pk_e, pcat = load_partition()
    qk_s, qk_e, qid, n_prom = load_promoters()

    regions = ([g for g in GB2UCSC] if naming == "gb" else MAIN)
    cmd = ["samtools", "view", "-@", "3", "-q", "30", "-f", "66", "-F", "1804",
           bam] + regions
    print(f"[count] {label} ({naming})", flush=True)
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, bufsize=1 << 22)

    cat_counts = np.zeros(len(CATS) + 1, dtype=np.int64)     # 1..4 used
    prom_counts = np.zeros(n_prom, dtype=np.int64)
    total = 0
    unmapped_chrom = 0

    reader = pd.read_csv(proc.stdout, sep="\t", header=None, usecols=[2, 3, 7, 8],
                         names=["rname", "pos", "pnext", "tlen"],
                         dtype={2: str, 3: np.int64, 7: np.int64, 8: np.int64},
                         chunksize=CHUNK, quoting=3, engine="c")
    for chunk in reader:
        chrom = chunk.rname.map(GB2UCSC) if naming == "gb" else chunk.rname
        ci = chrom.map(CIDX)
        ok = ci.notna().to_numpy()
        unmapped_chrom += int((~ok).sum())
        if not ok.any():
            continue
        ci = ci.to_numpy()[ok].astype(np.int64)
        pos = chunk.pos.to_numpy()[ok]
        pnext = chunk.pnext.to_numpy()[ok]
        tlen = np.abs(chunk.tlen.to_numpy()[ok])
        mid = np.minimum(pos, pnext) - 1 + (tlen // 2)
        key = ci * BIG + mid
        total += key.size

        # ---- category: disjoint, exactly one hit per fragment ----
        j = np.searchsorted(pk_s, key, side="right") - 1
        good = (j >= 0) & (key < pk_e[np.maximum(j, 0)])
        cat_counts += np.bincount(pcat[j[good]], minlength=len(CATS) + 1)

        # ---- promoters: overlapping, walk back up to MAXDEPTH ----
        j = np.searchsorted(qk_s, key, side="right") - 1
        for k in range(MAXDEPTH):
            jj = j - k
            m = (jj >= 0) & (key < qk_e[np.maximum(jj, 0)]) & (qk_s[np.maximum(jj, 0)] <= key)
            if not m.any():
                continue
            prom_counts += np.bincount(qid[jj[m]], minlength=n_prom)

    proc.stdout.close()
    if proc.wait() != 0:
        sys.exit(f"[fatal] samtools failed for {label}")

    got = int(cat_counts[1:].sum())
    if got != total:
        sys.exit(f"[fatal] {label}: categories sum to {got:,} but {total:,} fragments "
                 f"were placed — the partition is not covering the genome")

    with open(out_c + ".partial", "w") as fh:
        fh.write("category\tfragments\n")
        for i, c in enumerate(CATS, start=1):
            fh.write(f"{c}\t{cat_counts[i]}\n")
        fh.write(f"total\t{total}\n")
        fh.write(f"offchrom_skipped\t{unmapped_chrom}\n")
    os.replace(out_c + ".partial", out_c)

    nz = np.nonzero(prom_counts)[0]
    pd.DataFrame({"promoter_row": nz + 1, "fragments": prom_counts[nz]}).to_csv(
        out_p + ".partial", sep="\t", index=False, compression="gzip")
    os.replace(out_p + ".partial", out_p)

    print(f"[done] {label}  total={total:,}  promoter={cat_counts[1]:,} "
          f"({100*cat_counts[1]/total:.2f}%)  in-window sum={prom_counts.sum():,}",
          flush=True)


if __name__ == "__main__":
    main()
