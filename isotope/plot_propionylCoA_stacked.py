#!/usr/bin/env python3
"""Stacked 13C-labelled / unlabelled propionyl-CoA bars.

Generates three panels from three tidy sheets of the tracing workbook:

    Figure 4D   sh-NT vs sh-BCAT1          -> shBCAT1_propionylCoA_stacked.{png,svg,pdf}
    Figure 4E   stem-cell factors +/-      -> SCM_propionylCoA_stacked.{png,svg,pdf}
    Figure S4B  K27M vs KO, HPLM media     -> HPLM_propionylCoA_stacked.{png,svg,pdf}

Method (recovered and verified against the published values, see PROVENANCE_AUDIT.md):
each replicate's labelled and unlabelled propionyl-CoA is divided by the mean TOTAL
(labelled + unlabelled) of the control group in that comparison. The bar is the group
mean of the normalized total, split into a labelled (black) and unlabelled (grey)
segment; error bars are SEM of the normalized total; points are individual replicates.
Two Welch unequal-variance two-tailed t-tests are reported per panel -- one on the
summed pool, one on the labelled fraction -- with no multiple-comparison correction.

Run from anywhere:  python isotope/plot_propionylCoA_stacked.py
"""
import numpy as np, pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
from _isotope_common import (load, set_style, save, pfmt, C_LABELED, C_UNLABELED, HERE)
import os

set_style()
RNG = np.random.default_rng(11)

# sheet, group column, (control, test), x labels, y label, output stem, figure id
PANELS = [
    ("Plotting shBCAT1", "Construct", ("sh-NT", "sh-BCAT1"),
     ("sh-NT", "sh-BCAT1"), "Propionyl-CoA fractional abundance",
     "shBCAT1_propionylCoA_stacked", "Figure 4D", None),
    ("Plot SCM", "Media",
     ("Stem cell factors (+)", "Stem cell factors (-)"),
     ("SCM (+)", "SCM (−)"), "Propionyl-CoA fractional abundance",
     "SCM_propionylCoA_stacked", "Figure 4E", None),
    ("Plot HPLM", "Genotype", ("K27M", "KO"),
     ("K27M", "KO"), "Normalized propionyl-CoA abundance",
     "HPLM_propionylCoA_stacked", "Figure S4B", ("Basal media", "HPLM")),
]

LAB, UNLAB = "Isoleucine-C13 labeled", "Unlabeled"
rows = []

for sheet, gcol, (ctrl, test), xlabels, ylab, stem, figid, subset in PANELS:
    df = load(sheet)
    if subset:                      # Plot HPLM holds both media; S4B uses HPLM only
        df = df[df[subset[0]] == subset[1]]

    g = {k: df[df[gcol] == k][[LAB, UNLAB]].to_numpy(float)
         for k in (ctrl, test)}
    missing = [k for k, v in g.items() if len(v) == 0]
    if missing:
        raise SystemExit(f"{sheet}: no rows for {missing}")

    norm = g[ctrl].sum(1).mean()            # control-group mean total
    tot = {k: v.sum(1) / norm for k, v in g.items()}
    lab = {k: v[:, 0] / norm for k, v in g.items()}

    p_sum = stats.ttest_ind(tot[ctrl], tot[test], equal_var=False).pvalue
    p_lab = stats.ttest_ind(lab[ctrl], lab[test], equal_var=False).pvalue

    fig, ax = plt.subplots(figsize=(2.35, 3.0))
    for x, k in enumerate((ctrl, test)):
        mL, mT = lab[k].mean(), tot[k].mean()
        ax.bar(x, mL, width=0.62, color=C_LABELED, edgecolor="black",
               linewidth=0.8, zorder=2)
        ax.bar(x, mT - mL, bottom=mL, width=0.62, color=C_UNLABELED,
               edgecolor="black", linewidth=0.8, zorder=2)
        ax.errorbar(x, mT, yerr=tot[k].std(ddof=1) / np.sqrt(len(tot[k])),
                    fmt="none", ecolor="black", capsize=3, lw=0.9,
                    capthick=0.9, zorder=4)
        jit = RNG.uniform(-0.13, 0.13, len(tot[k]))
        ax.scatter(x + jit, tot[k], s=15, facecolor="white", edgecolor="black",
                   linewidth=0.6, zorder=5)
        ax.scatter(x + jit, lab[k], s=15, facecolor="black", edgecolor="black",
                   linewidth=0.6, zorder=5)

    top = max(tot[ctrl].max(), tot[test].max())
    y = top * 1.12
    ax.plot([0, 0, 1, 1], [y, y * 1.03, y * 1.03, y], color="black", lw=0.9)
    ax.text(0.5, y * 1.06, f"Sum {pfmt(p_sum)}\nLabel {pfmt(p_lab)}",
            ha="center", va="bottom", fontsize=7.5, linespacing=1.35)

    ax.set_xticks([0, 1]); ax.set_xticklabels(xlabels)
    ax.set_xlim(-0.62, 1.62); ax.set_ylim(0, top * 1.40)
    ax.set_ylabel(ylab)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(width=0.9, length=3)
    ax.tick_params(axis="x", length=0)

    h = [plt.Rectangle((0, 0), 1, 1, fc=C_LABELED, ec="black", lw=0.8),
         plt.Rectangle((0, 0), 1, 1, fc=C_UNLABELED, ec="black", lw=0.8)]
    fig.legend(h, ["Labeled ($^{13}$C)", "Unlabeled"], loc="upper center",
               bbox_to_anchor=(0.55, 1.06), ncol=1, frameon=False,
               fontsize=8, handlelength=1.1)

    save(fig, stem); plt.close(fig)

    for k in (ctrl, test):
        rows.append(dict(figure=figid, sheet=sheet, group=k, n=len(tot[k]),
                         mean_total=tot[k].mean(),
                         sem_total=tot[k].std(ddof=1) / np.sqrt(len(tot[k])),
                         mean_labeled=lab[k].mean(),
                         sem_labeled=lab[k].std(ddof=1) / np.sqrt(len(lab[k])),
                         p_sum=p_sum, p_label=p_lab))

st = pd.DataFrame(rows)
st["normalization"] = "each replicate / control-group mean total"
st["test"] = ("Welch unequal-variance two-tailed t-test, control vs test; "
              "reported separately for the summed pool and the labeled fraction; "
              "no multiple-comparison correction")
st.to_csv(os.path.join(HERE, "propionylCoA_stacked_stats.csv"), index=False)
print()
print(st.to_string(index=False))
