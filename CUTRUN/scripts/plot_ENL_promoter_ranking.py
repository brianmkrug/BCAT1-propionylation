#!/usr/bin/env python3
"""Figure 6G -- ranked ENL occupancy at protein-coding promoters ("super-ENL").

Input:  CUTRUN/analysis_QC_normalization/data/ENL_promoter_ranking_proteincoding.csv
Output: analysis_QC_normalization/plots/ENL_promoter_ranking_proteincoding.{png,svg,pdf}

Promoters ranked by IgG-subtracted ENL CPM. Promoters above the ROSE-style inflection
are flagged `super_ENL` in the input table and drawn in colour; glial-progenitor genes
(blue) and mevalonate/isoprenoid genes (orange) among them are labelled.

NOTE: the input table flags 473 super-ENL promoters. The manuscript text states 472 --
reconcile before submission.

Run from anywhere:  python CUTRUN/scripts/plot_ENL_promoter_ranking.py
"""
import os
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from _cutrun_common import QC, set_style, save

set_style()
GLIAL = ["MEIS2", "TUBA1C", "SOX2", "FABP7", "VIM"]
MEVAL = ["HMGCS1", "GGPS1", "HMGCR", "IDI1", "PDSS2"]
C_GLIAL, C_MEVAL, C_SUPER, C_REST = "#2b7bba", "#e2582a", "#9b6fc4", "#c9c9c9"

d = pd.read_csv(os.path.join(QC, "data/ENL_promoter_ranking_proteincoding.csv"))
d = d.sort_values("ENL_signal_IgGsub_CPM").reset_index(drop=True)
d["x"] = np.arange(len(d))
n_super = int(d.super_ENL.sum())
cut = d.loc[d.super_ENL.astype(bool), "x"].min()

fig, ax = plt.subplots(figsize=(4.2, 3.4))
rest = d[~d.super_ENL.astype(bool)]
ax.scatter(rest.x, rest.ENL_signal_IgGsub_CPM, s=3, color=C_REST,
           edgecolor="none", zorder=1)
sup = d[d.super_ENL.astype(bool)]
ax.scatter(sup.x, sup.ENL_signal_IgGsub_CPM, s=5, color=C_SUPER,
           edgecolor="none", zorder=2)
ax.axvline(cut, ls="--", color="0.45", lw=0.9, zorder=1)

for genes, col in ((GLIAL, C_GLIAL), (MEVAL, C_MEVAL)):
    sel = d[d.gene.isin(genes)].sort_values("ENL_signal_IgGsub_CPM", ascending=False)
    for k, (_, r) in enumerate(sel.iterrows()):
        ax.scatter(r.x, r.ENL_signal_IgGsub_CPM, s=34, color=col,
                   edgecolor="black", linewidth=0.5, zorder=4)
        ax.annotate(r.gene, (r.x, r.ENL_signal_IgGsub_CPM),
                    xytext=(-70, 6 - 11 * k if col == C_GLIAL else -12 - 11 * k),
                    textcoords="offset points", fontsize=7.5, color=col,
                    ha="right", va="center",
                    arrowprops=dict(arrowstyle="-", lw=0.6, color=col,
                                    shrinkA=0, shrinkB=2))

ax.set_xlabel("Promoter rank")
ax.set_ylabel("Normalized ENL signal")
ax.set_xlim(0, len(d) * 1.02)
ax.spines[["top", "right"]].set_visible(False)
ax.tick_params(width=0.9, length=3)
h = [plt.Line2D([], [], ls="", marker="o", ms=5, color=C_GLIAL),
     plt.Line2D([], [], ls="", marker="o", ms=5, color=C_MEVAL)]
ax.legend(h, ["Glial progenitor", "Mevalonate/isoprenoid"], frameon=False,
          fontsize=8, loc="upper left", handletextpad=0.3)
save(fig, os.path.join(QC, "plots", "ENL_promoter_ranking_proteincoding"))
print(f"   promoters={len(d)}  super-ENL={n_super} (manuscript says 472)  "
      f"inflection at rank {cut}")
