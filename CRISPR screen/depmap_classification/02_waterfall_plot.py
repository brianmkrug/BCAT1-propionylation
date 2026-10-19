#!/usr/bin/env python3
"""
Waterfall of the top 200 K27M fitness dependencies (most-negative MAGeCK beta),
coloured by DepMap dependency class.

Inputs : top200_annotated.csv  (from 01_classify_depmap.py)
Outputs: fig_K27M_top200_depmap_waterfall.{png,pdf,svg}
"""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

OUT = REPO_ROOT + "/CRISPR screen/depmap_classification"

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

df = pd.read_csv(f"{OUT}/top200_annotated.csv").sort_values("K27M.beta").reset_index(drop=True)
df["x"] = np.arange(len(df))
colors = df["depmap_class"].map(CLASS_COLOR)

fig = plt.figure(figsize=(14.0, 8.0))
gs = fig.add_gridspec(2, 3, height_ratios=[2.2, 1.5], width_ratios=[1.25, 1.0, 1.0],
                      hspace=0.30, wspace=0.34)

# ============================ A: waterfall ==================================
ax = fig.add_subplot(gs[0, :])
ax.bar(df["x"], df["K27M.beta"], color=colors, width=0.88, linewidth=0)
ax.axhline(0, color="black", lw=0.8)

ymin = df["K27M.beta"].min()
ax.set_ylim(ymin - 0.10, 0.78)

# greedy, non-overlapping labels for the strongest selective hits (+ BCAT1)
sel = df[df["depmap_class"] == "Strongly selective"].sort_values("K27M.beta")
chosen, taken = [], []
for _, r in sel.iterrows():
    if all(abs(r["x"] - t) >= 6 for t in taken):
        chosen.append(r); taken.append(r["x"])
    if len(chosen) >= 16:
        break
b1 = df[df["Gene"] == "BCAT1"]
if len(b1) and "BCAT1" not in [c["Gene"] for c in chosen]:
    chosen.append(b1.iloc[0])

for r in chosen:
    is_b = r["Gene"] == "BCAT1"
    ax.annotate(r["Gene"], xy=(r["x"], 0), xytext=(r["x"], 0.07),
                rotation=90, ha="center", va="bottom", fontsize=7.6,
                color="#c0392b" if is_b else "#8c2d04",
                fontweight="bold" if is_b else "normal")
    ax.plot([r["x"], r["x"]], [0, 0.06], color="#8c2d04", lw=0.55, alpha=0.75)

ax.set_xlim(-1.5, len(df) + 0.5)
ax.set_xlabel("Genes ranked by K27M fitness effect  (top 200 of 2,978 screened)", fontsize=10)
ax.set_ylabel("K27M MAGeCK MLE $\\beta$", fontsize=10.5)
ax.set_yticks([0, -0.5, -1.0, -1.5, -2.0])
ax.set_title("DIPGXIII K27M fitness dependencies stratified by DepMap dependency class",
             fontsize=13, fontweight="bold", pad=8)
ax.spines[["top", "right"]].set_visible(False)
ax.spines["bottom"].set_position(("data", ymin - 0.10))
ax.tick_params(labelsize=9)

counts = df["depmap_class"].value_counts()
handles = [Patch(facecolor=CLASS_COLOR[c], label=f"{c}  (n={counts.get(c, 0)}) — {CLASS_DESC[c]}")
           for c in CLASS_ORDER if counts.get(c, 0) > 0]
ax.legend(handles=handles, fontsize=8.6, loc="lower right", frameon=False,
          handlelength=1.3, borderaxespad=0.6)

# ================= B: top selective dependencies (readable) =================
axb = fig.add_subplot(gs[1, 0])
topsel = sel.head(18).iloc[::-1]
axb.barh(np.arange(len(topsel)), topsel["K27M.beta"],
         color=CLASS_COLOR["Strongly selective"], height=0.74, linewidth=0)
axb.set_yticks(np.arange(len(topsel)))
axb.set_yticklabels(topsel["Gene"], fontsize=7.4)
for t, g in zip(axb.get_yticklabels(), topsel["Gene"]):
    if g == "BCAT1":
        t.set_color("#c0392b"); t.set_fontweight("bold")
for i, (_, r) in enumerate(topsel.iterrows()):
    axb.text(r["K27M.beta"] - 0.03, i, f"{r['frac_dep_lines']*100:.0f}%",
             ha="right", va="center", fontsize=6.4, color="#555555")
axb.set_xlim(topsel["K27M.beta"].min() - 0.42, 0)
axb.set_xlabel("K27M $\\beta$   (grey = % DepMap lines dependent)", fontsize=8.6)
axb.set_title("Top selective (non-pan-essential) dependencies", fontsize=9.6, fontweight="bold")
axb.spines[["top", "right", "left"]].set_visible(False)
axb.tick_params(axis="x", labelsize=8)

# =================== C: class composition by rank bin =======================
axc = fig.add_subplot(gs[1, 1])
bins = [(0, 25), (25, 50), (50, 100), (100, 150), (150, 200)]
labels = ["1–25", "26–50", "51–100", "101–150", "151–200"]
bottom = np.zeros(len(bins))
for c in CLASS_ORDER:
    vals = np.array([((df["x"] >= a) & (df["x"] < b) & (df["depmap_class"] == c)).sum()
                     / (b - a) * 100 for a, b in bins])
    if vals.sum() == 0:
        continue
    axc.bar(labels, vals, bottom=bottom, color=CLASS_COLOR[c], width=0.74, linewidth=0)
    bottom += vals
axc.set_ylabel("% of genes", fontsize=9)
axc.set_xlabel("Rank bin (by K27M $\\beta$)", fontsize=9)
axc.set_ylim(0, 100)
axc.set_title("Class composition by rank", fontsize=9.6, fontweight="bold")
axc.spines[["top", "right"]].set_visible(False)
axc.tick_params(labelsize=8)
plt.setp(axc.get_xticklabels(), rotation=25, ha="right")

# ============ D: screen effect vs pan-cancer breadth ========================
axd = fig.add_subplot(gs[1, 2])
d2 = df[df["in_depmap"] == True]
axd.scatter(d2["frac_dep_lines"] * 100, d2["K27M.beta"], s=17,
            c=d2["depmap_class"].map(CLASS_COLOR), edgecolor="none", alpha=0.9)
axd.axvline(50, color="#999999", lw=0.7, ls="--")
if len(b1):
    axd.scatter(b1["frac_dep_lines"] * 100, b1["K27M.beta"], s=58, facecolor="none",
                edgecolor="#c0392b", lw=1.6, zorder=5)
    axd.annotate("BCAT1", (float(b1["frac_dep_lines"].iloc[0]) * 100,
                           float(b1["K27M.beta"].iloc[0])),
                 xytext=(7, 3), textcoords="offset points", fontsize=8,
                 color="#c0392b", fontweight="bold")
axd.text(0.02, 0.04, "DIPGXIII-selective", transform=axd.transAxes,
         fontsize=7.4, color="#8c2d04", style="italic")
axd.text(0.98, 0.04, "pan-essential", transform=axd.transAxes, ha="right",
         fontsize=7.4, color="#5a5a5a", style="italic")
axd.set_xlabel("% DepMap lines dependent", fontsize=9)
axd.set_ylabel("K27M $\\beta$", fontsize=9)
axd.set_title("Screen effect vs pan-cancer breadth", fontsize=9.6, fontweight="bold")
axd.spines[["top", "right"]].set_visible(False)
axd.tick_params(labelsize=8)

fig.text(0.006, 0.004,
         "DepMap Public 24Q4 (Chronos, 1,178 cell lines).  Common essential = CRISPRInferredCommonEssentials list;  "
         "strongly selective = normLRT > 100 (Azzalini skew-t vs normal fit to the cross-line gene-effect distribution) "
         "with ≥1 dependent line;  dependent line = gene dependency probability > 0.5.",
         fontsize=7, color="#555555")

fig.subplots_adjust(left=0.065, right=0.985, top=0.935, bottom=0.085)
for ext in ("png", "pdf", "svg"):
    fig.savefig(f"{OUT}/fig_K27M_top200_depmap_waterfall.{ext}",
                dpi=300 if ext == "png" else None, bbox_inches="tight")
print("saved")
