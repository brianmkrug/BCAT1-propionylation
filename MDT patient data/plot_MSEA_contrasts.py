#!/usr/bin/env python3
"""Figures 4G, S4D and S4E -- metabolite-set enrichment across patient tissue contrasts.

Input:  MDT patient data/contrast_comparison/MSEA_NES_comparison.csv
        MDT patient data/WTGBM_vs_Normal/MSEA_WTGBM_vs_Normal.csv
        MDT patient data/WTGBM_vs_IDHMUT/MSEA_WTGBM_vs_IDHMUT.csv
Output: contrast_comparison/MSEA_NES_scatter.{png,svg,pdf}    Figure 4G
        contrast_comparison/MSEA_bars_composite.{png,svg,pdf} Figures S4D + S4E

4G plots the NES of each Metabolon sub-pathway in WT HGG vs normal brain (x) against
WT HGG vs IDH-mutant HGG (y). Sets enriched in WT HGG in BOTH contrasts are highlighted
purple, sets depleted in both are highlighted green; the rest are grey.

S4D/S4E are the two contrasts drawn as ranked NES bars with the FDR q value printed on
each bar.

Run from anywhere:  python "MDT patient data/plot_MSEA_contrasts.py"
"""
import os
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from _mdt_common import HERE, set_style, save

set_style()
UP, DOWN, GREY = "#5b4b9e", "#2e8b3d", "#c8c8c8"

# ---------------- Figure 4G ----------------
cmp_ = pd.read_csv(os.path.join(HERE, "contrast_comparison/MSEA_NES_comparison.csv"))
x, y = cmp_["NES_WTGBM_vs_Normal"], cmp_["NES_WTGBM_vs_IDHMUT"]
both_up = (x > 0) & (y > 0)
both_dn = (x < 0) & (y < 0)

fig, ax = plt.subplots(figsize=(3.6, 3.4))
ax.axhline(0, color="0.5", lw=0.8); ax.axvline(0, color="0.5", lw=0.8)
lim = max(abs(np.r_[x, y])) * 1.12
ax.plot([-lim, lim], [-lim, lim], ls=":", color="0.75", lw=0.9, zorder=1)
ax.add_patch(plt.Rectangle((0, 0), lim, lim, fc=UP, alpha=0.07, zorder=0))
ax.add_patch(plt.Rectangle((-lim, -lim), lim, lim, fc=DOWN, alpha=0.07, zorder=0))
ax.scatter(x[~(both_up | both_dn)], y[~(both_up | both_dn)], s=22, color=GREY,
           edgecolor="none", zorder=2)
ax.scatter(x[both_up], y[both_up], s=26, color=UP, edgecolor="none", zorder=3)
ax.scatter(x[both_dn], y[both_dn], s=26, color=DOWN, edgecolor="none", zorder=3)
ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
ax.set_xlabel("NES: WT HGG vs Normal brain")
ax.set_ylabel("NES: WT HGG vs IDH-Mut HGG")
ax.text(0.97, 0.97, "IDH-WT HGG\nenriched", transform=ax.transAxes, ha="right",
        va="top", fontsize=7.5, color=UP, linespacing=1.3)
ax.text(0.03, 0.03, "IDH-WT HGG\ndepleted", transform=ax.transAxes, ha="left",
        va="bottom", fontsize=7.5, color=DOWN, linespacing=1.3)
ax.spines[["top", "right"]].set_visible(False)
ax.tick_params(width=0.9, length=3)
save(fig, os.path.join(HERE, "contrast_comparison/MSEA_NES_scatter"))
plt.close(fig)
print("  both-enriched sets: %d | both-depleted: %d | total: %d"
      % (both_up.sum(), both_dn.sum(), len(cmp_)))

# ---------------- Figures S4D + S4E ----------------
PANELS = [("WTGBM_vs_Normal/MSEA_WTGBM_vs_Normal.csv",
           "Subpathway MSEA: WT HGG vs Normal brain", "Normal brain", DOWN),
          ("WTGBM_vs_IDHMUT/MSEA_WTGBM_vs_IDHMUT.csv",
           "Subpathway MSEA: WT HGG vs IDH-Mut HGG", "IDH-Mut HGG", "#e07b39")]
TOPN = 10
fig, axes = plt.subplots(2, 1, figsize=(5.4, 6.4),
                         gridspec_kw=dict(hspace=0.42))
for ax, (path, title, comparator, ccol) in zip(axes, PANELS):
    d = pd.read_csv(os.path.join(HERE, path)).sort_values("NES", ascending=False)
    sel = pd.concat([d.head(TOPN // 2 + (TOPN % 2)), d.tail(TOPN // 2)])
    sel = sel.sort_values("NES")
    yy = np.arange(len(sel))
    ax.barh(yy, sel["NES"], height=0.68, edgecolor="none",
            color=[UP if v > 0 else ccol for v in sel["NES"]], zorder=2)
    for i, (_, r) in enumerate(sel.iterrows()):
        off = 0.06 if r["NES"] > 0 else -0.06
        ax.text(r["NES"] + off, i, f"FDR={r['FDR q-val']:.3f}", va="center",
                ha="left" if r["NES"] > 0 else "right", fontsize=6.5)
    ax.axvline(0, color="black", lw=0.9)
    ax.set_yticks(yy); ax.set_yticklabels(sel["Term"], fontsize=7.5)
    ax.set_title(title, fontsize=9, pad=4)
    ax.set_xlim(-3.4, 3.4)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(width=0.9, length=3)
    h = [plt.Rectangle((0, 0), 1, 1, fc=UP, ec="none"),
         plt.Rectangle((0, 0), 1, 1, fc=ccol, ec="none")]
    ax.legend(h, ["WT HGG", comparator], frameon=False, fontsize=7.5,
              loc="lower right", handlelength=1.1)
axes[1].set_xlabel("NES  (>0 enriched in WT HGG,  <0 enriched in comparator)",
                   fontsize=8.5)
save(fig, os.path.join(HERE, "contrast_comparison/MSEA_bars_composite"))
