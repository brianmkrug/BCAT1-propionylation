#!/usr/bin/env python3
"""Cell-culture sub-pathway GSEA bar plot, styled to match the xenograft-noHSJ figure:
top 7 up + top 7 down by NES, FDR printed opposite each bar, bars bold/saturated when
FDR<0.05 (muted otherwise), simplified title. Exports png + svg."""
import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# keep text as real, editable strings in SVG/PDF (not outlined paths)
plt.rcParams["svg.fonttype"] = "none"   # SVG: write <text> elements
plt.rcParams["pdf.fonttype"] = 42       # PDF: embed TrueType (editable in Illustrator)

HERE = os.path.dirname(os.path.abspath(__file__))
EN = os.path.join(HERE, "enrichment")
FIG = os.path.join(HERE, "figures")

d = pd.read_csv(os.path.join(EN, "GSEA_cellculture_subpathway.csv"))
d = d[~d["pathway"].str.contains("Food Component/Plant", case=False, na=False)]
up = d[d["NES"] > 0].nlargest(7, "NES")
dn = d[d["NES"] < 0].nsmallest(7, "NES")
sel = pd.concat([up, dn]).sort_values("NES").reset_index(drop=True)

SIG = 0.05
UP_BOLD, DN_BOLD = "#d1495b", "#2e86ab"
UP_MUTE, DN_MUTE = "#ecc4ca", "#c3dae7"


def fmt(f):
    return "FDR < 0.001" if f < 0.001 else f"FDR = {f:.3f}"


fig, ax = plt.subplots(figsize=(10, 7))
y = np.arange(len(sel))
sigmask = sel["FDR"] < SIG
colors = [(UP_BOLD if s else UP_MUTE) if v > 0 else (DN_BOLD if s else DN_MUTE)
          for v, s in zip(sel["NES"], sigmask)]
ax.barh(y, sel["NES"], color=colors, edgecolor="none", height=.7)

ax.set_yticks(y)
ax.set_yticklabels(sel["pathway"], fontsize=9)
for lab, s in zip(ax.get_yticklabels(), sigmask):
    lab.set_fontweight("bold" if s else "normal")
ax.axvline(0, color="k", lw=.9)

xmax = np.ceil(sel["NES"].abs().max() * 10) / 10 + 0.55
ax.set_xlim(-xmax, xmax)
pad = 0.06
for yi, (nes, fdr, s) in enumerate(zip(sel["NES"], sel["FDR"], sigmask)):
    weight = "bold" if s else "normal"
    col = "#000000" if s else "#666666"
    if nes > 0:
        ax.text(-pad, yi, fmt(fdr), va="center", ha="right", fontsize=8,
                color=col, fontweight=weight)
    else:
        ax.text(pad, yi, fmt(fdr), va="center", ha="left", fontsize=8,
                color=col, fontweight=weight)

ax.set_xlabel("NES   (> 0 enriched in K27M,   < 0 enriched in KO)\nbold = FDR < 0.05",
              fontsize=10)
ax.set_title("Cell culture sub-pathway GSEA (K27M vs KO)", fontweight="bold", fontsize=12)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.tight_layout()
for ext in ("png", "svg", "pdf"):
    out = os.path.join(FIG, f"GSEA_cellculture_subpathway.{ext}")
    fig.savefig(out, dpi=150)
    print("wrote", out)
