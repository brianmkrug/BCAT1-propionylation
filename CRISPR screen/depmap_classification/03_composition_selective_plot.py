#!/usr/bin/env python3
"""
Condensed v2 of fig_K27M_top200_depmap_waterfall: two panels only.

  A  Class composition by rank bin (stacked %), with the full class colour legend
  B  Top selective (non-pan-essential) dependencies

Dropped from v1: the 200-gene stratified waterfall and the
"screen effect vs pan-cancer breadth" scatter.

Inputs : top200_annotated.csv  (from 01_classify_depmap.py)
Outputs: fig_K27M_depmap_composition_selective.{png,pdf,svg}
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

# keep SVG/PDF text editable in Illustrator
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["pdf.fonttype"] = 42

OUT = Path(__file__).resolve().parent

CLASS_ORDER = ["Common essential", "Strongly selective", "Other dependency",
               "Not dependent", "Not in DepMap"]
CLASS_COLOR = {
    "Common essential":   "#5a5a5a",
    "Strongly selective": "#d94801",
    "Other dependency":   "#3182bd",
    "Not dependent":      "#9ecae1",
    "Not in DepMap":      "#d9d9d9",
}
CLASS_DESC = {
    "Common essential":   "pan-essential across DepMap lines",
    "Strongly selective": "normLRT > 100, selective dependency",
    "Other dependency":   "dependent in ≥1 line, not selective",
    "Not dependent":      "no dependent DepMap line",
    "Not in DepMap":      "not screened in DepMap 24Q4",
}
BCAT1_C = "#c0392b"
N_TOP_SELECTIVE = 18

df = (pd.read_csv(OUT / "top200_annotated.csv")
        .sort_values("K27M.beta").reset_index(drop=True))
df["x"] = np.arange(len(df))
counts = df["depmap_class"].value_counts()

fig = plt.figure(figsize=(11.0, 6.2))
gs = fig.add_gridspec(
    2, 2, height_ratios=[1.0, 0.26], width_ratios=[1.0, 1.0],
    hspace=0.30, wspace=0.26,
    left=0.075, right=0.98, top=0.90, bottom=0.09,
)

# ==================== A: class composition by rank bin ======================
axa = fig.add_subplot(gs[0, 0])
bins = [(0, 25), (25, 50), (50, 100), (100, 150), (150, 200)]
labels = ["1–25", "26–50", "51–100", "101–150", "151–200"]

bottom = np.zeros(len(bins))
for c in CLASS_ORDER:
    vals = np.array([((df["x"] >= a) & (df["x"] < b) & (df["depmap_class"] == c)).sum()
                     / (b - a) * 100 for a, b in bins])
    if vals.sum() == 0:
        continue
    axa.bar(labels, vals, bottom=bottom, color=CLASS_COLOR[c],
            width=0.76, linewidth=0)
    bottom += vals

axa.set_ylabel("% of genes", fontsize=10)
axa.set_xlabel("Rank bin (by K27M $\\beta$)", fontsize=10)
axa.set_ylim(0, 100)
axa.set_yticks([0, 20, 40, 60, 80, 100])
axa.set_title("Class composition by rank", fontsize=11.5, fontweight="bold", pad=8)
axa.spines[["top", "right"]].set_visible(False)
axa.tick_params(labelsize=9)

# ---- class colour legend, directly under the panel it explains -------------
axleg = fig.add_subplot(gs[1, 0])
axleg.axis("off")
handles = [Patch(facecolor=CLASS_COLOR[c],
                 label=f"{c}  (n={counts.get(c, 0)}) — {CLASS_DESC[c]}")
           for c in CLASS_ORDER if counts.get(c, 0) > 0]
axleg.legend(handles=handles, fontsize=7.6, loc="upper left",
             bbox_to_anchor=(0.0, 1.25), frameon=False,
             handlelength=1.2, handleheight=1.0, labelspacing=0.42,
             borderaxespad=0.0, ncol=1)

# ============ B: top selective (non-pan-essential) dependencies =============
axb = fig.add_subplot(gs[0:, 1])
sel = df[df["depmap_class"] == "Strongly selective"].sort_values("K27M.beta")
topsel = sel.head(N_TOP_SELECTIVE).iloc[::-1]

axb.barh(np.arange(len(topsel)), topsel["K27M.beta"],
         color=CLASS_COLOR["Strongly selective"], height=0.74, linewidth=0)
axb.set_yticks(np.arange(len(topsel)))
axb.set_yticklabels(topsel["Gene"], fontsize=9)
for t, g in zip(axb.get_yticklabels(), topsel["Gene"]):
    if g == "BCAT1":
        t.set_color(BCAT1_C)
        t.set_fontweight("bold")

for i, (_, r) in enumerate(topsel.iterrows()):
    axb.text(r["K27M.beta"] - 0.04, i, f"{r['frac_dep_lines'] * 100:.0f}%",
             ha="right", va="center", fontsize=7.6, color="#555555")

axb.set_ylim(-0.7, len(topsel) - 0.3)
axb.set_xlim(topsel["K27M.beta"].min() - 0.45, 0)
axb.set_xlabel("K27M $\\beta$   (grey = % DepMap lines dependent)", fontsize=10)
axb.set_title("Top selective (non-pan-essential) dependencies",
              fontsize=11.5, fontweight="bold", pad=8)
axb.spines[["top", "right", "left"]].set_visible(False)
axb.tick_params(axis="x", labelsize=9)
axb.tick_params(axis="y", length=0)

fig.text(0.006, 0.006,
         "DepMap Public 24Q4 (Chronos, 1,178 cell lines).  Common essential = CRISPRInferredCommonEssentials list;  "
         "strongly selective = normLRT > 100 (Azzalini skew-t vs normal fit to the cross-line gene-effect distribution) "
         "with ≥1 dependent line;  dependent line = gene dependency probability > 0.5.",
         fontsize=6.6, color="#555555")

for ext in ("png", "pdf", "svg"):
    p = OUT / f"fig_K27M_depmap_composition_selective.{ext}"
    fig.savefig(p, dpi=300 if ext == "png" else None,
                bbox_inches="tight", facecolor="white")
    print("wrote", p)
