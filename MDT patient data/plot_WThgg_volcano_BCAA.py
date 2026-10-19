#!/usr/bin/env python3
"""Figure 4H -- tissue metabolite volcano, IDH-WT HGG vs pooled normal brain + IDH-mutant.

Input:  MDT patient data/WThgg_vs_pooled_classifier/WThgg_vs_pooled_ranked.csv
Output: WThgg_vs_pooled_classifier/volcano_WThgg_vs_pooled_BCAA.{png,svg,pdf}

x = log2 fold change (WT HGG / pooled comparator), y = -log10 Welch p.
BCAA-pathway and BCAA-adjacent metabolites are drawn in green: filled where they pass
both thresholds, open where they do not. Everything else is grey. Dashed guides mark
p = 0.05 and |log2FC| = 1.

Run from anywhere:  python "MDT patient data/plot_WThgg_volcano_BCAA.py"
"""
import os
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from _mdt_common import HERE, set_style, save

set_style()
GREEN, GREY = "#2e8b3d", "#c4c4c4"
P_CUT, FC_CUT = 0.05, 1.0

# sub-pathways counted as BCAA-pathway / BCAA-adjacent
BCAA_SUB = ["Leucine, Isoleucine and Valine Metabolism",
            "Fatty Acid Metabolism (also BCAA Metabolism)",
            "Gamma-glutamyl Amino Acid"]
LABEL = ["gamma-glutamylisoleucine*", "propionylcarnitine (C3)",
         "2-methylbutyrylcarnitine (C5)"]

d = pd.read_csv(os.path.join(HERE,
    "WThgg_vs_pooled_classifier/WThgg_vs_pooled_ranked.csv")).dropna(
    subset=["log2FC_WT_vs_pooled", "p_welch"])
d["nlp"] = -np.log10(d["p_welch"])
d["is_bcaa"] = d.SUB_PATHWAY.isin(BCAA_SUB)
d["sig"] = (d.p_welch < P_CUT) & (d.log2FC_WT_vs_pooled.abs() > FC_CUT)

fig, ax = plt.subplots(figsize=(4.0, 3.6))
o = d[~d.is_bcaa]
ax.scatter(o.log2FC_WT_vs_pooled, o.nlp, s=9, color=GREY, edgecolor="none", zorder=2)
ns = d[d.is_bcaa & ~d.sig]
ax.scatter(ns.log2FC_WT_vs_pooled, ns.nlp, s=26, facecolor="none",
           edgecolor=GREEN, linewidth=0.9, zorder=3, label="BCAA, n.s.")
s = d[d.is_bcaa & d.sig]
ax.scatter(s.log2FC_WT_vs_pooled, s.nlp, s=30, color=GREEN, edgecolor="none",
           zorder=4, label="BCAA, sig")

for name in LABEL:
    r = d[d.BIOCHEMICAL == name]
    if len(r):
        r = r.iloc[0]
        ax.annotate(name.rstrip("*"), (r.log2FC_WT_vs_pooled, r.nlp),
                    textcoords="offset points", xytext=(-8, 7), ha="right",
                    fontsize=7.5, color=GREEN)

ax.axhline(-np.log10(P_CUT), ls="--", color="0.6", lw=0.8)
for v in (-FC_CUT, FC_CUT):
    ax.axvline(v, ls="--", color="0.6", lw=0.8)
lim = np.ceil(d.log2FC_WT_vs_pooled.abs().max())
ax.set_xlim(-lim, lim)
ax.set_xlabel("log2 fold-change (WT HGG / pooled Normal+IDH-Mut)")
ax.set_ylabel("-log10 p (Welch)")
ax.spines[["top", "right"]].set_visible(False)
ax.tick_params(width=0.9, length=3)
ax.legend(frameon=False, fontsize=8, loc="upper left", handletextpad=0.3)
save(fig, os.path.join(HERE, "WThgg_vs_pooled_classifier/volcano_WThgg_vs_pooled_BCAA"))
print("  metabolites: %d | BCAA-pathway: %d (%d significant at p<%.2f & |log2FC|>%.0f)"
      % (len(d), d.is_bcaa.sum(), (d.is_bcaa & d.sig).sum(), P_CUT, FC_CUT))
print(s[["BIOCHEMICAL", "SUB_PATHWAY", "log2FC_WT_vs_pooled", "p_welch"]]
      .sort_values("p_welch").head(10).to_string(index=False))
