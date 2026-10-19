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

d = pd.read_csv(os.path.join(EN, "GSEA_xenograft_subpathway_noHSJ.csv"))
d = d[~d["pathway"].str.contains("Food Component/Plant", case=False, na=False)]

up = d[d["NES"] > 0].nlargest(7, "NES")
dn = d[d["NES"] < 0].nsmallest(7, "NES")
sel = pd.concat([up, dn]).sort_values("NES").reset_index(drop=True)  # bottom=most neg
SPECIAL = {
    "Leucine, Isoleucine and Valine Metabolism": "BCAA Metab",
    "Long Chain Polyunsaturated Fatty Acid (n3 and n6)": "PUFA",
    "Glycolysis, Gluconeogenesis, and Pyruvate Metabolism": "Glycol, Gluconeogen & Pyruvate Metab",
}
def simp(n):
    return SPECIAL.get(n, n.replace("Metabolism", "Metab"))
sel["label"] = sel["pathway"].map(simp)

SIG = 0.05
UP_BOLD, DN_BOLD = "#d1495b", "#2e86ab"     # saturated: FDR < 0.05
UP_MUTE, DN_MUTE = "#ecc4ca", "#c3dae7"     # muted: not significant


def fmt(f):
    if f < 0.001:
        return "FDR < 0.001"
    return f"FDR = {f:.3f}"


fig, ax = plt.subplots(figsize=(10/2.54, 5/2.54))
y = np.arange(len(sel))
sigmask = sel["FDR"] < SIG
colors = [(UP_BOLD if s else UP_MUTE) if v > 0 else (DN_BOLD if s else DN_MUTE)
          for v, s in zip(sel["NES"], sigmask)]
ax.barh(y, sel["NES"], color=colors, edgecolor="#555555", linewidth=.35, height=.74)

ax.set_yticks(y)
ax.set_yticklabels(sel["label"], fontsize=4)
for lab, s in zip(ax.get_yticklabels(), sigmask):
    lab.set_fontweight("bold" if s else "normal")
ax.axvline(0, color="k", lw=.6)

xmax = sel["NES"].abs().max() + 0.08
ax.set_xlim(-xmax, xmax)
pad = 0.04
for yi, (nes, fdr, s) in enumerate(zip(sel["NES"], sel["FDR"], sigmask)):
    weight = "bold" if s else "normal"
    col = "#000000" if s else "#666666"
    if nes > 0:                                     # bar right -> FDR text on left
        ax.text(-pad, yi, fmt(fdr), va="center", ha="right", fontsize=3.3,
                color=col, fontweight=weight)
    else:                                           # bar left -> FDR text on right
        ax.text(pad, yi, fmt(fdr), va="center", ha="left", fontsize=3.3,
                color=col, fontweight=weight)

ax.set_xlabel("NES  (>0 K27M,  <0 KO;  bold = FDR<0.05)", fontsize=4.5)
ax.tick_params(labelsize=4, length=1.5, width=0.4)
ax.set_title("Xenograft subpathway MSEA", fontweight="bold", fontsize=5.5, pad=3)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.tight_layout()
out = os.path.join(FIG, "GSEA_xenograft_subpathway_noHSJ.png")
fig.savefig(out, dpi=600, bbox_inches="tight"); fig.savefig(out.replace(".png",".pdf"), bbox_inches="tight")
print("wrote", out)
