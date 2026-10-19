#!/usr/bin/env python3
"""Figure 7E -- promoter H3K9pr in BCAAr vs control mouse allograft tumours.

Input:  CUTRUN/mm10_analysis/data/FAmetab_H3K9pr.csv
        CUTRUN/mm10_analysis/data/counts_10kb_matrix.npz  (background density)
Output: mm10_analysis/plots/mm10_FAmetab_H3K9pr_2HIT.{png,svg,pdf}

Grey hexagons are genome-wide 10 kb bins, log2(H3K9pr/IgG) in control (x) against BCAAr
(y). Coloured points are the 198 lipid-metabolism gene promoters in FAmetab_H3K9pr.csv,
coloured by their RNA log2 fold change in the 4-hit allograft. The dashed line is y = x;
points below it lost H3K9pr on BCAA restriction.

Run from anywhere:  python CUTRUN/scripts/plot_mm10_FAmetab_H3K9pr.py
"""
import os
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from _cutrun_common import MM10, set_style, save, bins_log2_vs_igg

set_style()
LABEL = ["Elovl6", "Ptprg", "Ube2l6"]

marks = bins_log2_vs_igg(os.path.join(MM10, "data/counts_10kb_matrix.npz"), "control")
marks_b = bins_log2_vs_igg(os.path.join(MM10, "data/counts_10kb_matrix.npz"), "BCAAr")
xk = [k for k in marks if "H3K9pr" in k][0]
bx, by = marks[xk], marks_b[[k for k in marks_b if "H3K9pr" in k][0]]
ok = np.isfinite(bx) & np.isfinite(by)

d = pd.read_csv(os.path.join(MM10, "data/FAmetab_H3K9pr.csv"))
norm = TwoSlopeNorm(vmin=-1.5, vcenter=0.0, vmax=1.5)

fig, ax = plt.subplots(figsize=(4.0, 3.6))
ax.hexbin(bx[ok], by[ok], gridsize=60, cmap="Greys", bins="log", mincnt=1,
          linewidths=0, zorder=1)
sc = ax.scatter(d.enrich_ctrl, d.enrich_BCAAr, c=d.RNA_log2FC_4HIT, cmap="RdBu_r",
                norm=norm, s=34, edgecolor="black", linewidth=0.5, zorder=3)
lim = (-2.2, 4.0)
ax.plot(lim, lim, ls="--", color="#4a6fa5", lw=0.9, zorder=2)
for g in LABEL:
    r = d[d.gene == g]
    if len(r):
        r = r.iloc[0]
        ax.annotate(g, (r.enrich_ctrl, r.enrich_BCAAr), textcoords="offset points",
                    xytext=(7, 2), fontsize=8, style="italic")
ax.set_xlim(*lim); ax.set_ylim(*lim)
ax.set_xlabel("control H3K9pr log2(H3K9pr/IgG)")
ax.set_ylabel("BCAAr H3K9pr log2(H3K9pr/IgG)")
ax.spines[["top", "right"]].set_visible(False)
ax.tick_params(width=0.9, length=3)
cb = fig.colorbar(sc, ax=ax, pad=0.02, fraction=0.046)
cb.set_label("Lipid metabolism pathway RNA log2FC", fontsize=7.5)
cb.ax.tick_params(labelsize=7.5, width=0.8, length=2.5)
save(fig, os.path.join(MM10, "plots", "mm10_FAmetab_H3K9pr_2HIT"))
print(f"   {len(d)} lipid-metabolism promoters; "
      f"{int((d.enrich_BCAAr < d.enrich_ctrl).sum())} below the diagonal "
      f"(lost H3K9pr on BCAAr); background bins n={int(ok.sum())}")
