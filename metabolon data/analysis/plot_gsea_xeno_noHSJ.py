#!/usr/bin/env python3
"""Revised xenograft (outlier-removed) sub-pathway GSEA bar plot.
Drops Food Component/Plant; shows top 7 up + top 7 down by NES; FDR printed in the
blank space opposite each bar. Muted directional colors — no significance highlight
yet (user will choose an FDR threshold to bold)."""
import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42


HERE = os.path.dirname(os.path.abspath(__file__))
EN = os.path.join(HERE, "enrichment")
FIG = os.path.join(HERE, "figures")
CM = 1 / 2.54
plt.rcParams.update({"svg.fonttype": "none"})

d = pd.read_csv(os.path.join(EN, "GSEA_xenograft_subpathway_noHSJ.csv"))
d = d[~d["pathway"].str.contains("Food Component/Plant", case=False, na=False)]

up = d[d["NES"] > 0].nlargest(7, "NES")
dn = d[d["NES"] < 0].nsmallest(7, "NES")
sel = pd.concat([up, dn]).sort_values("NES").reset_index(drop=True)  # bottom=most neg

SIG = 0.05
UP_BOLD, DN_BOLD = "#d1495b", "#2e86ab"     # saturated: FDR < 0.05
UP_MUTE, DN_MUTE = "#ecc4ca", "#c3dae7"     # muted: not significant


def fmt(f):
    if f < 0.001:
        return "FDR < 0.001"
    return f"FDR = {f:.3f}"


FS = 6  # base text size (pt)
fig, ax = plt.subplots(figsize=(10 * CM, 6 * CM))
y = np.arange(len(sel))
sigmask = sel["FDR"] < SIG
colors = [(UP_BOLD if s else UP_MUTE) if v > 0 else (DN_BOLD if s else DN_MUTE)
          for v, s in zip(sel["NES"], sigmask)]
ax.barh(y, sel["NES"], color=colors, edgecolor="none", height=.72)

ax.set_yticks(y)
ax.set_yticklabels(sel["pathway"], fontsize=FS)
for lab, s in zip(ax.get_yticklabels(), sigmask):
    lab.set_fontweight("bold" if s else "normal")
ax.axvline(0, color="k", lw=.4)

xmax = np.ceil(sel["NES"].abs().max() * 10) / 10 + 0.55
ax.set_xlim(-xmax, xmax)
pad = 0.06
for yi, (nes, fdr, s) in enumerate(zip(sel["NES"], sel["FDR"], sigmask)):
    weight = "bold" if s else "normal"
    col = "#000000" if s else "#666666"
    ha = "right" if nes > 0 else "left"
    ax.text(-pad if nes > 0 else pad, yi, fmt(fdr), va="center", ha=ha,
            fontsize=FS - 0.5, color=col, fontweight=weight)

ax.set_xlabel("NES   (> 0 K27M,   < 0 KO)\nbold = FDR < 0.05", fontsize=FS)
ax.set_title("Xenograft sub-pathway GSEA (K27M vs KO)", fontweight="bold", fontsize=FS + 0.5)
ax.tick_params(axis="both", labelsize=FS, width=0.3, length=1.3)
ax.set_xticks([-2, -1, 0, 1, 2])
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)
for sp in ("left", "bottom"):
    ax.spines[sp].set_linewidth(0.3)
fig.tight_layout()
for ext in ("png", "svg"):
    out = os.path.join(FIG, f"GSEA_xenograft_subpathway_noHSJ.{ext}")
    fig.savefig(out, dpi=600, bbox_inches="tight")
    print("wrote", out)
