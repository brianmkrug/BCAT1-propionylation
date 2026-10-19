#!/usr/bin/env python3
"""Figures 6H, 6I and S6G -- promoter mark loss on BCAT1 knockdown.

Input:  CUTRUN/analysis_QC_normalization/data/delta_promoter_table.csv
        CUTRUN/analysis_QC_normalization/data/DEG_down_padj05.csv
Output: analysis_QC_normalization/plots/deltaENL_spikein_DEGdown_H3K9pr.{png,svg,pdf}  Fig 6H
        analysis_QC_normalization/plots/human_deltaK9pr_violin_DEG.{png,svg,pdf}       Fig 6I
        analysis_QC_normalization/plots/deltaENL_spikein_DEGdown_H3K4me3.{png,svg,pdf} Fig S6G
        analysis_QC_normalization/plots/delta_promoter_verification.csv

Spike-in scale factors are applied as printed on the published axes:
    H3K9pr   shBCAT1 x 0.62
    ENL      shBCAT1 x 0.80
i.e. log2(0.62) and log2(0.80) are added to the stored delta columns. H3K4me3 is
plotted unscaled, matching the S6G axis label.

    *** PROVENANCE WARNING -- READ BEFORE USING ***
The stored delta_promoter_table.csv does NOT reproduce the statistics printed on the
published panels. With this table Figure 6H's both-down quadrant contains 40 genes
(the panel prints 130) and Figure 6I gives p = 1.8e-04 by Welch t (the panel prints
2e-07). The medians also differ. The published figures were therefore made from a
different version of this table -- most likely a larger promoter universe -- which is
not in the repository. This script prints a verification block on every run; do not
use its output to replace the published panels until that is resolved.

Run from anywhere:  python CUTRUN/scripts/plot_delta_promoter_panels.py
"""
import os
import numpy as np, pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from _cutrun_common import QC, set_style, save

set_style()
SCALE = {"dH3K9pr_log2FC": 0.62, "dENL_log2FC": 0.80, "dH3K4me3_log2FC": 1.0}
RED, GREY = "#e03131", "#efe9df"

d = pd.read_csv(os.path.join(QC, "data/delta_promoter_table.csv"))
deg = set(pd.read_csv(os.path.join(QC, "data/DEG_down_padj05.csv")).gene)
d["down"] = d.gene.isin(deg)
for c, f in SCALE.items():
    d[c] = d[c] + np.log2(f)

PUB = {"6H_both_down": 130, "6I_p": 2e-07, "S6G_n_red": 146, "super_ENL": 472}
got = {"6H_both_down": int(((d.dENL_log2FC < 0) & (d.dH3K9pr_log2FC < 0) & d.down).sum()),
       "S6G_n_red": int(d.down.sum()),
       "super_ENL": int(d.super_ENL.sum())}
a = d.loc[d.down, "dH3K9pr_log2FC"].dropna()
b = d.loc[~d.down, "dH3K9pr_log2FC"].dropna()
got["6I_p"] = stats.ttest_ind(a, b, equal_var=False).pvalue

print("=" * 66)
print("VERIFICATION against the values printed on the published panels")
for k in PUB:
    ok = "OK " if (abs(np.log10(got[k] / PUB[k])) < 0.05 if k == "6I_p"
                   else got[k] == PUB[k]) else "!! MISMATCH"
    print(f"  {ok}  {k:16s} published={PUB[k]:<10} this table={got[k]}")
print("  -> see the provenance warning in this script's docstring")
print("=" * 66)
pd.DataFrame([{"quantity": k, "published": PUB[k], "from_stored_table": got[k]}
              for k in PUB]).to_csv(
    os.path.join(QC, "plots", "delta_promoter_verification.csv"), index=False)


def hexpanel(ycol, ylab, xlab, stem, nlab):
    fig, ax = plt.subplots(figsize=(3.7, 3.4))
    s = d.dropna(subset=["dH3K9pr_log2FC" if "H3K9pr" in xlab else "dH3K4me3_log2FC",
                         ycol])
    xc = "dH3K9pr_log2FC" if "H3K9pr" in xlab else "dH3K4me3_log2FC"
    hb = ax.hexbin(s[xc], s[ycol], gridsize=45, bins="log", cmap="viridis",
                   mincnt=1, linewidths=0, zorder=2)
    r = s[s.down]
    ax.scatter(r[xc], r[ycol], s=16, color=RED, edgecolor="black", linewidth=0.35,
               zorder=4, label=f"downregulated (padj<0.05), n={len(r)}")
    ax.axhline(0, color="0.3", lw=0.9); ax.axvline(0, color="0.3", lw=0.9)
    lim = 1.6
    ax.add_patch(plt.Rectangle((-lim, -lim), lim, lim, fc=GREY, zorder=0))
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_xlabel(xlab); ax.set_ylabel(ylab)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(width=0.9, length=3)
    if nlab:
        ax.text(0.04, 0.04, nlab, transform=ax.transAxes, color=RED, fontsize=8)
    ax.legend(frameon=False, fontsize=7.5, loc="lower center",
              bbox_to_anchor=(0.5, -0.30), markerscale=1.0)
    cb = fig.colorbar(hb, ax=ax, pad=0.02, fraction=0.046)
    cb.set_label("promoters per hexagon (log)", fontsize=7.5)
    cb.ax.tick_params(labelsize=7, width=0.8, length=2.5)
    save(fig, os.path.join(QC, "plots", stem)); plt.close(fig)


hexpanel("dENL_log2FC", "Δ ENL log2(shBCAT1·0.80 / shNT)",
         "Δ H3K9pr log2(shBCAT1·0.62 / shNT)",
         "deltaENL_spikein_DEGdown_H3K9pr", f"both↓ {got['6H_both_down']}")
hexpanel("dENL_log2FC", "Δ ENL log2(shBCAT1·0.80 / shNT)",
         "Δ H3K4me3 log2(shBCAT1 / shNT) [CPM]",
         "deltaENL_spikein_DEGdown_H3K4me3", None)

# ---------------- Figure 6I ----------------
fig, ax = plt.subplots(figsize=(2.6, 3.2))
parts = ax.violinplot([a.values, b.values], positions=[0, 1], showextrema=False,
                      widths=0.8)
for pc, c in zip(parts["bodies"], [RED, "#b9b9b9"]):
    pc.set_facecolor(c); pc.set_alpha(1.0); pc.set_edgecolor("black")
    pc.set_linewidth(0.7)
for i, v in enumerate([a.values, b.values]):
    q1, med, q3 = np.percentile(v, [25, 50, 75])
    ax.vlines(i, q1, q3, color="black", lw=4, zorder=3)
    ax.scatter(i, med, s=14, color="white", zorder=4, edgecolor="black",
               linewidth=0.5)
p = got["6I_p"]
ax.text(0.5, 1.01, f"p={p:.0e}", transform=ax.transAxes, ha="center", fontsize=8)
ax.set_xticks([0, 1]); ax.set_xticklabels(["down", "unchanged"])
ax.set_ylabel("Δ H3K9pr log2(shBCAT1·0.62 / shNT)")
ax.spines[["top", "right"]].set_visible(False)
ax.tick_params(width=0.9, length=3); ax.tick_params(axis="x", length=0)
save(fig, os.path.join(QC, "plots", "human_deltaK9pr_violin_DEG"))
print(f"\n   6I: down n={len(a)} median {a.median():+.3f} | "
      f"unchanged n={len(b)} median {b.median():+.3f} | Welch p={p:.3g}")
