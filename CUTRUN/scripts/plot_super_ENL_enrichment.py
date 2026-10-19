#!/usr/bin/env python3
"""Figure S6F -- pathway enrichment among super-ENL promoters.

Input:  CUTRUN/analysis_QC_normalization/data/enrich_metabolic_inducibleSets.csv
Output: analysis_QC_normalization/plots/
        ENL_super_metabolic_enrichment_inducibleSets.{png,svg,pdf}

Two-sided Fisher exact test of each metabolic / Hallmark gene set against the
super-ENL promoter set, restricted to sets with more than ten members in the universe.
Bars are odds ratios; green marks FDR < 0.05 (Benjamini-Hochberg), grey n.s. The
dashed line is an odds ratio of 1.

Run from anywhere:  python CUTRUN/scripts/plot_super_ENL_enrichment.py
"""
import os
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from _cutrun_common import QC, set_style, save

set_style()
SIG, NS = "#2e8b3d", "#b8b8b8"
TOPN = 14

d = pd.read_csv(os.path.join(QC, "data/enrich_metabolic_inducibleSets.csv"))
d = d[d.n_set_in_universe > 10].sort_values("odds_ratio", ascending=False).head(TOPN)
d = d.iloc[::-1]

fig, ax = plt.subplots(figsize=(4.6, 3.6))
y = np.arange(len(d))
ax.barh(y, d.odds_ratio, height=0.68, edgecolor="none",
        color=[SIG if q < 0.05 else NS for q in d.FDR], zorder=2)
ax.axvline(1, ls="--", color="0.4", lw=0.9)
ax.set_yticks(y)
ax.set_yticklabels([s.split(":", 1)[-1] if ":" in s else s for s in d.set],
                   fontsize=7.5)
ax.set_xlabel("odds ratio (super-ENL vs background,\nFisher exact, gene sets over 10)",
              fontsize=8)
ax.spines[["top", "right"]].set_visible(False)
ax.tick_params(width=0.9, length=3)
h = [plt.Rectangle((0, 0), 1, 1, fc=SIG, ec="none"),
     plt.Rectangle((0, 0), 1, 1, fc=NS, ec="none")]
ax.legend(h, ["FDR<0.05", "n.s."], frameon=False, fontsize=8, loc="lower right",
          handlelength=1.1)
save(fig, os.path.join(QC, "plots", "ENL_super_metabolic_enrichment_inducibleSets"))
print(d[["set", "n_set_in_universe", "overlap_obs", "odds_ratio", "FDR"]]
      .iloc[::-1].to_string(index=False))
