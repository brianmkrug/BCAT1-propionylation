#!/usr/bin/env python
"""
Revision of the genome-wide 10 kb correlation heatmap, restricted to the six
chromatin-state targets: CAGE and IgG dropped, magma colour scale, no dendrogram.

Subsetting the stored matrix is exact, not an approximation: a Spearman
correlation between two profiles depends only on those two profiles, so dropping
CAGE and IgG cannot change any remaining entry. Nothing is recomputed.

Row/column order is still the average-linkage clustering order -- the dendrogram
is omitted from the drawing, not from the ordering, and the subtitle says so.
Pass --fixed-order to use the hypothesis order instead
(H3K4me3, ENL, H3K9pr | H3K27ac, H3K4me1, CTCF) if you would rather the figure
not let the algorithm choose the arrangement.
"""
from __future__ import annotations
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root

import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform

QC = os.path.expanduser(
    REPO_ROOT + "/CUTRUN/analysis_QC_normalization/multimark_2026-07")
SRC = os.path.join(QC, "data", "spearman_withCTCF_K4me1_CAGE_10kb.csv")
DROP = ["CAGE", "IgG"]
FIXED = ["H3K4me3", "ENL", "H3K9pr", "H3K27ac", "H3K4me1", "CTCF"]
ASSAY = {"H3K4me3": "CUT&RUN", "ENL": "CUT&RUN", "H3K9pr": "CUT&RUN",
         "H3K27ac": "CUT&RUN", "H3K4me1": "CUT&RUN", "CTCF": "CUT&RUN"}
LABCOL = {"H3K4me3": "#2166ac", "ENL": "#1e8449", "H3K9pr": "#6c3483",
          "H3K27ac": "#e8862a", "H3K4me1": "#d94a8c", "CTCF": "#7d6608"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fixed-order", action="store_true")
    ap.add_argument("--cmap", default="magma")
    a = ap.parse_args()

    S = pd.read_csv(SRC, index_col=0)
    S = S.drop(index=DROP, columns=DROP)
    labels = list(S.index)
    M = S.to_numpy()

    if a.fixed_order:
        order = [labels.index(x) for x in FIXED if x in labels]
        note = "order: hypothesis groups (H3K4me3/ENL/H3K9pr | H3K27ac/H3K4me1/CTCF)"
    else:
        Z = linkage(squareform(np.clip(1 - M, 0, None), checks=False), method="average")
        order = dendrogram(Z, no_plot=True)["leaves"]
        note = "order: average-linkage clustering (dendrogram omitted)"
    Mx = M[np.ix_(order, order)]
    lx = [labels[i] for i in order]

    fig, ax = plt.subplots(figsize=(7.4, 6.6))
    vmin = np.floor(Mx[~np.eye(len(lx), dtype=bool)].min() * 20) / 20
    im = ax.imshow(Mx, cmap=a.cmap, vmin=vmin, vmax=1.0)

    cmap = matplotlib.colormaps[a.cmap]
    for i in range(len(lx)):
        for j in range(len(lx)):
            r, g, b, _ = cmap((Mx[i, j] - vmin) / max(1.0 - vmin, 1e-9))
            lum = 0.299 * r + 0.587 * g + 0.114 * b   # pick text colour by cell luminance
            ax.text(j, i, f"{Mx[i, j]:.2f}", ha="center", va="center", fontsize=10,
                    color="black" if lum > 0.55 else "white")

    # target name only: no assay suffix, no per-mark colour
    ax.set_xticks(range(len(lx))); ax.set_xticklabels(lx, rotation=90, fontsize=10,
                                                      color="black")
    ax.set_yticks(range(len(lx))); ax.set_yticklabels(lx, fontsize=10, color="black")
    ax.set_xticks(np.arange(-.5, len(lx), 1), minor=True)
    ax.set_yticks(np.arange(-.5, len(lx), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.4)
    ax.tick_params(which="minor", length=0)
    for s in ax.spines.values():
        s.set_visible(False)

    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    cb.set_label("Spearman r (10 kb bins)", fontsize=9)
    cb.ax.tick_params(labelsize=8)
    # no title above, no caption below -- the ordering note is printed to stdout
    # instead so it stays on the record without appearing on the figure
    fig.tight_layout()

    suffix = "_fixedorder" if a.fixed_order else ""
    out = os.path.join(QC, "plots", f"corr_heatmap_chromatin_states_10kb{suffix}.png")
    fig.savefig(out, dpi=200)
    fig.savefig(out.replace(".png", ".pdf"))
    fig.savefig(out.replace(".png", ".svg"))
    plt.close(fig)
    print("[plot]", out, "(+ .pdf, .svg)")
    print("[note]", note)
    print("\nmatrix as drawn:")
    print(pd.DataFrame(Mx, index=lx, columns=lx).round(3).to_string())


if __name__ == "__main__":
    main()
