#!/usr/bin/env python3
"""Figure 7A (charts) -- amino-acid and macronutrient composition of the two diets.

Input:  diet data/BCAAs diet composition.xlsx
          sheet 'AAs'          Amino acid (%), Control, BCAAr
          sheet 'composition'  Macronutrient (%), Control, BCAAr
Output: diet_composition_piecharts.{png,svg,pdf}

Two nested donuts, one per diet: the outer ring is the amino-acid breakdown of the
protein fraction, the inner ring is the protein / fat / carbohydrate split. The BCAAs
(leucine, isoleucine, valine) are drawn in the three greens of the palette and go to
zero in the BCAAr diet, with the removed mass redistributed over the other amino acids
(chiefly glycine) so that total protein is matched at 15.2%.

Run from anywhere:  python "diet data/plot_diet_composition.py"
"""
import os, json
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def repo_root(start=None):
    d = os.path.abspath(start or os.path.dirname(__file__))
    while not os.path.exists(os.path.join(d, "BCAT1_figure_palette.json")):
        p = os.path.dirname(d)
        if p == d:
            raise RuntimeError("repository root not found")
        d = p
    return d


ROOT = repo_root(); HERE = os.path.join(ROOT, "diet data")
XL = os.path.join(HERE, "BCAAs diet composition.xlsx")
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans"],
                     "font.size": 8, "svg.fonttype": "none", "pdf.fonttype": 42,
                     "ps.fonttype": 42})

aa = pd.read_excel(XL, sheet_name="AAs").dropna(how="all")
mac = pd.read_excel(XL, sheet_name="composition").dropna(how="all")
aa.columns = ["amino_acid", "Control", "BCAAr"]
mac.columns = ["macronutrient", "Control", "BCAAr"]
for c in ("Control", "BCAAr"):
    aa[c] = pd.to_numeric(aa[c]); mac[c] = pd.to_numeric(mac[c])

# BCAAs first (green ramp), then the rest in sheet order
BCAA = ["Leucine", "Isoleucine", "Valine"]
BCAA_COL = ["#1b7a3a", "#3fa45c", "#8ed08b"]
rest = [a for a in aa.amino_acid if a not in BCAA]
order = BCAA + rest
cmap = plt.get_cmap("tab20")
COL = dict(zip(BCAA, BCAA_COL))
COL.update({a: cmap(i % 20) for i, a in enumerate(rest)})
MAC_COL = {"Protein": "#e6c260", "Fat": "#e8e8e8", "Carbohydrates": "#3d4a5c"}

fig, axes = plt.subplots(1, 2, figsize=(6.0, 3.2))
for ax, diet in zip(axes, ("Control", "BCAAr")):
    a = aa.set_index("amino_acid").loc[order, diet]
    m = mac.set_index("macronutrient")[diet]
    ax.pie(a.values, radius=1.0, colors=[COL[x] for x in order],
           wedgeprops=dict(width=0.30, edgecolor="white", linewidth=0.6),
           startangle=90, counterclock=False)
    ax.pie(m.values, radius=0.68, colors=[MAC_COL[x] for x in m.index],
           wedgeprops=dict(width=0.30, edgecolor="white", linewidth=0.6),
           startangle=90, counterclock=False)
    ax.set_title(diet, fontsize=11, fontweight="bold",
                 color="#2e8b3d" if diet == "Control" else "#f0a020", pad=6)
    ax.set(aspect="equal")

h_aa = [plt.Rectangle((0, 0), 1, 1, fc=COL[a], ec="none") for a in order]
h_m = [plt.Rectangle((0, 0), 1, 1, fc=MAC_COL[x], ec="none") for x in mac.macronutrient]
l1 = fig.legend(h_m, list(mac.macronutrient), loc="upper right",
                bbox_to_anchor=(1.34, 0.95), frameon=False, fontsize=8,
                title="Macronutrients", title_fontsize=8.5, handlelength=1.1)
fig.add_artist(l1)
fig.legend(h_aa, order, loc="lower center", bbox_to_anchor=(0.5, -0.34), ncol=4,
           frameon=False, fontsize=7.5, handlelength=1.1, columnspacing=1.2)

for e, kw in (("png", dict(dpi=400)), ("svg", {}), ("pdf", {})):
    fig.savefig(os.path.join(HERE, f"diet_composition_piecharts.{e}"),
                bbox_inches="tight", **kw)
print("wrote diet_composition_piecharts.{png,svg,pdf}")
print("\nprotein total (%%): Control %.3f  BCAAr %.3f"
      % (aa.Control.sum() * 100, aa.BCAAr.sum() * 100))
print(mac.to_string(index=False))
