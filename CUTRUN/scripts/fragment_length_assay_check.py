#!/usr/bin/env python
"""
CORRECTION (2026-10-05): CTCF and H3K4me1 are CUT&RUN from shNT DIPG-XIII cells
(Diagenode C15410210; CST 5326). The xChIP/ChIPRx/DMSO/K27M-EV strings in the original
fastq names are legacy labels, not sample metadata; text below that calls them ChIP or a
different treatment arm is superseded.

Read-level test of assay type: fragment-length distribution.

The assignment of CTCF and H3K4me1 to "ChIP" originally came from METADATA -- the
fastq names preserved in bwa's @PG line say `xChIP` and `ChIPRx`. That is a label,
not evidence from the data. This script asks whether the reads themselves say the
same thing, independently of any filename.

The discriminator: CUT&RUN cuts with tethered MNase, which is size-constrained --
it leaves a sub-nucleosomal population and a tight distribution truncated near
mononucleosome size. Crosslinked ChIP fragments by sonication, which is
structure-agnostic and leaves a broad right-shifted smear.

Statistics reported (from the TLEN histogram, chr20):
  P(TLEN < 120)   sub-nucleosomal -- higher in CUT&RUN
  P(TLEN > 250)   sonication tail -- THE cleanest separator here
  mode, SD        distribution position and spread

A nucleosomal-periodicity autocorrelation was tried first and is NOT reported: it
returned 0 for every library including the known CUT&RUN ones, because these
particular libraries show no strong di-nucleosome peak. It discriminated nothing,
so keeping it would have been decoration.
"""
from __future__ import annotations
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root

import os
import subprocess

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = os.path.expanduser(REPO_ROOT + "/CUTRUN")
QC = os.path.join(ROOT, "analysis_QC_normalization", "multimark_2026-07")
CHROM = "chr20"
MAXLEN = 600

SAMPLES = [
    ("H3K4me3 (CUT&RUN)",  "DIPG13-shNT-H3K4me3-r1_nodup_sorted_unique.bam", "#2166ac"),
    ("H3K27ac (CUT&RUN)",  "shNT-H3K27ac-r1_nodup_sorted_unique.bam",        "#e8862a"),
    ("H3K9pr (CUT&RUN)",   "DIPG13-shNT-H3K9pr-r1_nodup_sorted_unique.bam",  "#6c3483"),
    ("IgG (CUT&RUN)",      "DIPG13-shNT-IgG-r1_nodup_sorted_unique.bam",     "#808080"),
    ("H3K4me1 (?)",        "K4me1_nodup_sorted_unique.bam",                  "#d94a8c"),
    ("CTCF (?)",           "CTCF_nodup_sorted_unique.bam",                   "#7d6608"),
]
# H3K27ac r1 is GenBank-named
GB = {"shNT-H3K27ac-r1_nodup_sorted_unique.bam": "CM000682.2"}


def tlen_hist(bam):
    region = GB.get(bam, CHROM)
    p = subprocess.run(
        ["samtools", "view", "-q", "30", "-f", "66", "-F", "1804",
         os.path.join(ROOT, "hg38", bam), region],
        capture_output=True, text=True)
    h = np.zeros(MAXLEN + 1)
    for line in p.stdout.splitlines():
        t = abs(int(line.split("\t", 9)[8]))
        if 0 < t <= MAXLEN:
            h[t] += 1
    return h / max(h.sum(), 1)


def shape(h):
    """Sonication smear vs MNase ladder, without assuming a visible ladder.

    An autocorrelation-based periodicity score was tried first and returned 0 for
    every library including the known CUT&RUN ones -- these particular libraries
    do not show a strong di-nucleosome peak, so that statistic discriminates
    nothing and is not reported. What does separate them is where the mass sits:
    MNase digestion leaves a large sub-nucleosomal population and a mode at or
    below mononucleosome size; sonication leaves a broad right-shifted smear.
    """
    x = np.arange(len(h))
    sub = float(h[:120].sum())
    longf = float(h[250:].sum())
    mode = int(np.argmax(h))
    mu = float((x * h).sum())
    sd = float(np.sqrt(((x - mu) ** 2 * h).sum()))
    return sub, longf, mode, sd


def main():
    fig, axes = plt.subplots(1, 2, figsize=(13.0, 4.9),
                             gridspec_kw={"width_ratios": [1.55, 1]})
    rows = []
    ax = axes[0]
    for name, bam, col in SAMPLES:
        h = tlen_hist(bam)
        sub, longf, mode, sd = shape(h)
        rows.append((name, sub, longf, mode, sd))
        ls = "--" if "(?)" in name else "-"
        ax.plot(np.arange(MAXLEN + 1), h, lw=1.7, color=col, ls=ls, label=name)
        print(f"  {name:20s} <120bp {100*sub:5.1f}%   >250bp {100*longf:5.1f}%   "
              f"mode {mode:3d} bp   SD {sd:5.1f}")
    ax.set_xlim(0, MAXLEN)
    ax.set_xlabel("fragment length (TLEN, bp)")
    ax.set_ylabel("fraction of fragments")
    ax.set_title(f"Fragment-length distribution ({CHROM})\n"
                 "solid = known CUT&RUN, dashed = the two new files", fontsize=10)
    ax.legend(fontsize=7.6, frameon=False)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    ax = axes[1]
    for name, sub, longf, mode, sd in rows:
        col = dict((n, c) for n, _, c in SAMPLES)[name]
        ax.scatter(100 * sub, 100 * longf, s=150, color=col,
                   marker="^" if "(?)" in name else "o",
                   edgecolor="k", linewidth=0.7, zorder=3)
        ax.annotate(name.split(" (")[0], (100 * sub, 100 * longf), fontsize=8,
                    textcoords="offset points", xytext=(8, 4))
    ax.set_xlabel("sub-nucleosomal fraction, TLEN < 120 bp (%)")
    ax.set_ylabel("long-fragment fraction, TLEN > 250 bp (%)")
    ax.set_title("Assay signature\nMNase (CUT&RUN): short-shifted\n"
                 "sonication (ChIP): long-shifted", fontsize=10)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    out = os.path.join(QC, "plots", "fragment_length_assay_check.png")
    fig.savefig(out, dpi=200); fig.savefig(out.replace(".png", ".pdf"))
    print("[plot]", out)

    with open(os.path.join(QC, "data", "fragment_length_assay_check.csv"), "w") as fh:
        fh.write("sample,pct_under_120bp,pct_over_250bp,mode_bp,sd_bp\n")
        for name, sub, longf, mode, sd in rows:
            fh.write(f"\"{name}\",{100*sub:.2f},{100*longf:.2f},{mode},{sd:.1f}\n")


if __name__ == "__main__":
    main()
