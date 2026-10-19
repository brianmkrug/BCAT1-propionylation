#!/usr/bin/env python
"""
CORRECTION (2026-10-05): CTCF and H3K4me1 are CUT&RUN from shNT DIPG-XIII cells
(Diagenode C15410210; CST 5326). The xChIP/ChIPRx/DMSO/K27M-EV strings in the original
fastq names are legacy labels, not sample metadata; text below that calls them ChIP or a
different treatment arm is superseded.

Full genomic-category panel: promoter / exonic / intronic / intergenic for every
target, now including CTCF, H3K4me1 and CAGE.

Categories are the disjoint partition from build_genome_partition.sh, so the four
shares are one pie per target and a promoter base that is also exonic counts once,
as promoter. Enrichment = (% of signal) / (% of genome), so 1.0 is chance.

CUT&RUN marks are shown as the pooled shNT profile. CTCF and H3K4me1 are single
shNT CUT&RUN libraries (no replicate, no shBCAT1 counterpart) -- they position a
target in this space but carry no delta.
CAGE (S.K27M) is 1 bp 5'-CTSS signal, so its "share" is a share of tag mass.
"""
from __future__ import annotations
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = os.path.expanduser(REPO_ROOT + "/CUTRUN")
QC = os.path.join(ROOT, "analysis_QC_normalization", "multimark_2026-07")
PC, ANN = os.path.join(QC, "promcat"), os.path.join(QC, "annotation")
CATS = ["promoter", "exonic", "intronic", "intergenic"]
CAT_COLOR = {"promoter": "#c0392b", "exonic": "#e8862a",
             "intronic": "#2166ac", "intergenic": "#95a5a6"}
COLOR = {"H3K4me3": "#2166ac", "ENL": "#1e8449", "H3K9pr": "#6c3483",
         "H3K27ac": "#e8862a", "H3K4me1": "#d94a8c", "CTCF": "#7d6608",
         "IgG": "#808080", "CAGE": "#c0392b"}
# display order: promoter-focused -> background
ORDER = ["CAGE", "H3K4me3", "H3K9pr", "CTCF", "H3K27ac", "H3K4me1", "ENL", "IgG"]
POOLS = {"H3K4me3": ["shNT-H3K4me3-r1", "shNT-H3K4me3-r2"],
         "H3K27ac": ["shNT-H3K27ac-r1", "shNT-H3K27ac-r2"],
         "H3K9pr": ["shNT-H3K9pr-r1", "shNT-H3K9pr-r2"],
         "ENL": ["shNT-ENL-r1", "shNT-ENL-r2"],
         "IgG": ["shNT-IgG-r1"], "CTCF": ["CTCF"], "H3K4me1": ["H3K4me1"]}
ASSAY = {"H3K4me3": "CUT&RUN", "ENL": "CUT&RUN", "H3K9pr": "CUT&RUN",
         "H3K27ac": "CUT&RUN", "IgG": "CUT&RUN",
         "H3K4me1": "CUT&RUN", "CTCF": "CUT&RUN", "CAGE": "CAGE"}


def genome_fractions():
    p = pd.read_csv(os.path.join(ANN, "partition.bed"), sep="\t", header=None,
                    names=["c", "s", "e", "cat"])
    bp = p.assign(w=p.e - p.s).groupby("cat").w.sum()
    return {c: bp[i + 1] / bp.sum() for i, c in enumerate(CATS)}


def main():
    g = genome_fractions()
    print("[genome] " + "  ".join(f"{c} {100*g[c]:.2f}%" for c in CATS))
    pct = {}
    for t in ORDER:
        if t == "CAGE":
            d = pd.read_csv(os.path.join(QC, "cage", "S.K27M.cage_categories.tsv"),
                            sep="\t", index_col=0)
            pct[t] = {c: float(d.loc[c, "pct"]) for c in CATS}
            continue
        agg, tot = {c: 0 for c in CATS}, 0
        for lab in POOLS[t]:
            d = pd.read_csv(os.path.join(PC, f"{lab}.category.tsv"), sep="\t",
                            index_col=0).fragments
            for c in CATS:
                agg[c] += int(d[c])
            tot += int(d["total"])
        pct[t] = {c: 100 * agg[c] / tot for c in CATS}

    rows = ["target,assay," + ",".join(f"pct_{c}" for c in CATS) + "," +
            ",".join(f"enrich_{c}" for c in CATS)]
    print("\n=== signal composition | enrichment vs genome ===")
    for t in ORDER:
        e = {c: pct[t][c] / (100 * g[c]) for c in CATS}
        print(f"  {t:8s} ({ASSAY[t]:8s})  " +
              "  ".join(f"{c[:4]} {pct[t][c]:5.2f}% ({e[c]:5.2f}x)" for c in CATS))
        rows.append(f"{t},{ASSAY[t]}," + ",".join(f"{pct[t][c]:.4f}" for c in CATS)
                    + "," + ",".join(f"{e[c]:.4f}" for c in CATS))
    with open(os.path.join(QC, "data", "category_composition_full.csv"), "w") as fh:
        fh.write("\n".join(rows) + "\n")

    fig, axes = plt.subplots(1, 3, figsize=(17.6, 5.2),
                             gridspec_kw={"width_ratios": [1.25, 1.35, 1.0]})
    x = np.arange(len(ORDER))
    lab = [f"{t}\n({ASSAY[t]})" for t in ORDER]

    ax = axes[0]
    bottom = np.zeros(len(ORDER))
    for c in CATS:
        v = np.array([pct[t][c] for t in ORDER])
        ax.bar(x, v, bottom=bottom, color=CAT_COLOR[c], label=c,
               edgecolor="white", linewidth=0.6)
        bottom += v
    for i, c in enumerate(CATS):
        ax.axhline(sum(100 * g[CATS[j]] for j in range(i + 1)), color="k",
                   ls=":", lw=1.0, zorder=5)
    ax.set_xticks(x); ax.set_xticklabels(lab, fontsize=8)
    ax.set_ylabel("% of signal"); ax.set_ylim(0, 100)
    ax.set_title("Signal composition\ndotted = genome composition (chance)", fontsize=10)
    ax.legend(fontsize=7.4, frameon=False, ncol=4, loc="upper center",
              bbox_to_anchor=(0.5, -0.13))

    ax = axes[1]
    w = 0.2
    for k, c in enumerate(CATS):
        ax.bar(x + (k - 1.5) * w, [pct[t][c] / (100 * g[c]) for t in ORDER], w,
               color=CAT_COLOR[c], label=c)
    ax.axhline(1, color="k", lw=1.0, ls="--")
    ax.set_yscale("log")
    ax.set_xticks(x); ax.set_xticklabels(lab, fontsize=8)
    ax.set_ylabel("enrichment vs genome (×, log)")
    ax.set_title("Category enrichment — 1.0 = chance\n"
                 "log scale: CAGE promoter enrichment dwarfs the marks", fontsize=10)
    ax.legend(fontsize=7.4, frameon=False, ncol=2)

    ax = axes[2]
    pe = [pct[t]["promoter"] / (100 * g["promoter"]) for t in ORDER]
    ax.barh(x, pe, color=[COLOR[t] for t in ORDER], edgecolor="k", linewidth=0.5)
    ax.axvline(1, color="k", ls="--", lw=1.0)
    for i, v in enumerate(pe):
        ax.text(v * 1.03, i, f"{v:.1f}×", va="center", fontsize=8)
    ax.set_yticks(x); ax.set_yticklabels(lab, fontsize=8)
    ax.invert_yaxis(); ax.set_xscale("log")
    ax.set_xlabel("promoter enrichment (×, log)")
    ax.set_title("Promoter focus, ranked", fontsize=10)
    for a in axes:
        for s in ("top", "right"):
            a.spines[s].set_visible(False)
    fig.tight_layout()
    out = os.path.join(QC, "plots", "genomic_category_composition_full.png")
    fig.savefig(out, dpi=200); fig.savefig(out.replace(".png", ".pdf"))
    print("\n[plot]", out)


if __name__ == "__main__":
    main()
