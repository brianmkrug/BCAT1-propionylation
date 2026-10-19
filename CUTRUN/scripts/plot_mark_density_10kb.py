#!/usr/bin/env python3
"""Figures 6B and S6A -- genome-wide co-occurrence of H3K9pr with ENL and H3K4me3.

Input:  CUTRUN/analysis_QC_normalization/data/counts_10kb_matrix.npz
Output: analysis_QC_normalization/plots/H3K9pr_vs_ENL_density_10kb.{png,svg,pdf}     Fig 6B
        analysis_QC_normalization/plots/H3K9pr_vs_H3K4me3_density_10kb.{png,svg,pdf} Fig S6A

Counts in 303,114 genome-wide 10 kb bins are CPM-normalized, a pseudocount of 1.0 CPM
is added, and each mark is expressed as log2(mark / IgG) within the shNT condition.
Every bin is used -- no coverage filter. Density is shown as a log-scaled 2D histogram
with the least-squares line overlaid; the annotated r is the Pearson correlation.

VERIFIED: this recipe returns r = 0.750 (6B) and r = 0.831 (S6A), matching the values
printed on the published panels. The script asserts both.

Run from anywhere:  python CUTRUN/scripts/plot_mark_density_10kb.py
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.stats import pearsonr
from _cutrun_common import QC, set_style, save, bins_log2_vs_igg

set_style()
marks = bins_log2_vs_igg(os.path.join(QC, "data/counts_10kb_matrix.npz"), "shNT")

PANELS = [("ENL", "H3K9pr_vs_ENL_density_10kb", 0.75, "magma"),
          ("H3K4me3", "H3K9pr_vs_H3K4me3_density_10kb", 0.83, "magma")]

for other, stem, expect, cmap in PANELS:
    x, y = marks["H3K9pr"], marks[other]
    ok = np.isfinite(x) & np.isfinite(y)
    r = pearsonr(x[ok], y[ok])[0]
    assert abs(r - expect) < 0.01, \
        f"{other}: r={r:.3f} but the published panel prints {expect}"
    fig, ax = plt.subplots(figsize=(3.6, 3.2))
    h = ax.hist2d(x[ok], y[ok], bins=180, norm=LogNorm(vmin=1), cmap=cmap)
    m, b = np.polyfit(x[ok], y[ok], 1)
    xs = np.array([x[ok].min(), x[ok].max()])
    ax.plot(xs, m * xs + b, color="#8fd3f4", lw=1.0, zorder=4)
    ax.text(0.04, 0.95, f"r={r:.2f}", transform=ax.transAxes, va="top", fontsize=9)
    ax.set_xlabel("mean log2(H3K9pr / IgG)")
    ax.set_ylabel(f"mean log2({other} / IgG)")
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(width=0.9, length=3)
    cb = fig.colorbar(h[3], ax=ax, pad=0.02, fraction=0.046)
    cb.set_label("10 kb bins", fontsize=8)
    cb.ax.tick_params(labelsize=7.5, width=0.8, length=2.5)
    save(fig, os.path.join(QC, "plots", stem))
    plt.close(fig)
    print(f"   {other}: n={ok.sum()} bins, Pearson r={r:.3f}")
