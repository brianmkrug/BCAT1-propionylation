#!/usr/bin/env python3
"""Figure 4I -- BCAA-derived metabolites and acylcarnitines across patient tissue groups.

Input:  MDT patient data/PFA_ tissue edited copy.xlsx, sheet 'MassNormImpData (tissue)'
Output: barplots/MDT_ggIle_ketoIle_carnitines_panel.{png,svg,pdf}
        barplots/MDT_ggIle_ketoIle_carnitines_stats.csv

Three blocks, left to right:
  1. gamma-glutamylisoleucine and ketoisoleucine (3-methyl-2-oxovalerate)
  2. acylcarnitines assigned to 'Leucine, Isoleucine and Valine Metabolism'
  3. acylcarnitines assigned to 'Fatty Acid Metabolism (also BCAA Metabolism)'

Bars are group means of the mass-normalized, imputed abundance with SEM error bars and
individual patients overlaid. n = 5 normal brain, 3 IDH-Mut HGG, 3 WT HGG. p values are
Welch unequal-variance two-tailed t-tests of WT HGG against normal brain on LOG2-
transformed abundances, uncorrected. (The log2 step matters: testing on the raw scale
gives systematically larger p values and does not reproduce the published figure.)

Run from anywhere:  python "MDT patient data/plot_ggIle_ketoIle_carnitines.py"
"""
import os
import numpy as np, pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
from _mdt_common import HERE, ORDER, COL, load_tissue, set_style, save

set_style(8.5)
RNG = np.random.default_rng(9)

BLOCKS = [
    (None, ["gamma-glutamylisoleucine*", "3-methyl-2-oxovalerate"],
     ["γ-glutamyl-\nisoleucine", "ketoisoleucine"]),
    ("BCAA catabolism",
     ["2-methylbutyrylcarnitine (C5)", "tiglylcarnitine (C5:1-DC)",
      "isobutyrylcarnitine (C4)", "isovalerylcarnitine (C5)",
      "beta-hydroxyisovaleroylcarnitine"],
     ["2-methylbutyryl (C5)", "tiglyl (C5:1-DC)", "isobutyryl (C4)",
      "isovaleryl (C5)", "β-OH-isovaleroyl"]),
    ("FA / also BCAA",
     ["acetylcarnitine (C2)", "propionylcarnitine (C3)", "butyrylcarnitine (C4)"],
     ["acetyl (C2)", "propionyl (C3)", "butyryl (C4)"]),
]

long = load_tissue()
have = set(long.metabolite)
rows, panels = [], []
for header, mets, labels in BLOCKS:
    for m, lab in zip(mets, labels):
        if m not in have:
            print("  !! not found in workbook:", m); continue
        panels.append((header, m, lab))

fig, axes = plt.subplots(1, len(panels), figsize=(1.05 * len(panels) + 0.9, 2.5),
                         gridspec_kw=dict(wspace=0.72))
for ax, (header, m, lab) in zip(np.atleast_1d(axes), panels):
    sub = long[long.metabolite == m]
    vals = {g: sub[sub.group == g]["value"].to_numpy(float) for g in ORDER}
    for x, g in enumerate(ORDER):
        v = vals[g]
        ax.bar(x, v.mean(), width=0.68, color=COL[g], edgecolor="none",
               yerr=v.std(ddof=1) / np.sqrt(len(v)) if len(v) > 1 else None,
               error_kw=dict(elinewidth=0.9, capsize=2.5, capthick=0.9), zorder=2)
        ax.scatter(x + RNG.uniform(-0.15, 0.15, len(v)), v, s=11,
                   facecolor="white", edgecolor="black", linewidth=0.5, zorder=4)
    # Welch on log2-transformed abundances -- this reproduces the published p values
    p = stats.ttest_ind(np.log2(vals["WT HGG"]), np.log2(vals["Normal brain"]),
                        equal_var=False).pvalue
    hi = max(v.max() for v in vals.values())
    ax.text(0.5, 1.02, f"p={p:.1g}" if p < 0.01 else f"p={p:.3f}",
            transform=ax.transAxes, ha="center", va="bottom", fontsize=7)
    ax.set_ylim(0, hi * 1.28)
    ax.set_xticks(range(len(ORDER))); ax.set_xticklabels([])
    ax.set_xlim(-0.65, len(ORDER) - 0.35)
    ax.set_title(lab, fontsize=7.5, pad=13)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(width=0.9, length=3, labelsize=7)
    ax.tick_params(axis="x", length=0)
    if ax is np.atleast_1d(axes)[0]:
        ax.set_ylabel("norm. abund.", fontsize=8)
    rows.append(dict(figure="Figure 4I", block=header or "BCAA metabolites",
                     metabolite=m,
                     **{f"mean_{g}": vals[g].mean() for g in ORDER},
                     **{f"n_{g}": len(vals[g]) for g in ORDER},
                     p_WTHGG_vs_Normal=p))

h = [plt.Rectangle((0, 0), 1, 1, fc=COL[g], ec="none") for g in ORDER]
fig.legend(h, ORDER, loc="lower center", bbox_to_anchor=(0.5, -0.14), ncol=3,
           frameon=False, fontsize=8, handlelength=1.1, columnspacing=1.6)
save(fig, os.path.join(HERE, "barplots/MDT_ggIle_ketoIle_carnitines_panel"))

st = pd.DataFrame(rows)
st["test"] = ("Welch unequal-variance two-tailed t-test on log2-transformed values, "
              "WT HGG vs Normal brain; uncorrected. Bars show mass-normalized "
              "imputed abundance on the linear scale.")
st.to_csv(os.path.join(HERE, "barplots/MDT_ggIle_ketoIle_carnitines_stats.csv"),
          index=False)
print(st[["metabolite", "mean_Normal brain", "mean_IDH-Mut HGG", "mean_WT HGG",
          "p_WTHGG_vs_Normal"]].round(4).to_string(index=False))
