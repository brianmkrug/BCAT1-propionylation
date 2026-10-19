#!/usr/bin/env python3
"""Column dot plot of the top-depleted Braindex signature from each progenitor/glial and
myeloid family in 2HIT BCAAr vs control.

Single column of dots (no x-axis); dot colour = NES with a horizontal NES scale below.
Glial / progenitor signatures on top, immune / myeloid below (group labels sit next to the
signatures). All signatures are at the permutation FDR floor (q<0.001) — note in the legend.
Values from the comprehensive Braindex run. PNG + SVG + PDF.
"""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))  # publication: repo-relative root
import json
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["pdf.fonttype"] = 42
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize

ROOT = REPO_ROOT + ""
FIG = f"{ROOT}/BCAAr_analysis/figures"
with open(f"{ROOT}/BCAT1_figure_palette.json") as f:
    CMAP = json.load(f)["heatmap"]["colormap"]        # RdBu_r

GLIAL = [("P-e13_GLIP", "glial progenitor"),
         ("F-e18_OPC",  "oligo precursor"),
         ("F-e16_RGC-P", "radial glia")]
IMMUNE = [("P-p6_MGL", "microglia"),
          ("F-p0_MAC", "macrophage")]

bx = pd.read_csv(f"{ROOT}/BCAAr_analysis/gsea/Braindex_mouse_2HIT.csv").set_index("Term")
def nes(t): return float(bx.loc[t, "NES"])

rows, y = [], 6.0
for sig, ct in GLIAL:
    rows.append((sig, ct, y, "glial")); y -= 1.0
y -= 0.8
for sig, ct in IMMUNE:
    rows.append((sig, ct, y, "immune")); y -= 1.0
df = pd.DataFrame(rows, columns=["sig", "ct", "y", "group"])
df["NES"] = df.sig.map(nes)
vmax = float(df.NES.abs().max())
norm = Normalize(-vmax, vmax); sm = cm.ScalarMappable(norm=norm, cmap=CMAP)

fig = plt.figure(figsize=(5.6, 4.3))
ax = fig.add_axes([0.03, 0.20, 0.95, 0.76]); ax.axis("off")
XD = 0.0
ax.scatter([XD] * len(df), df.y, s=440, c=df.NES, cmap=CMAP, vmin=-vmax, vmax=vmax,
           edgecolor="0.25", linewidth=0.7, zorder=3)

LBL_X = -0.10
for _, r in df.iterrows():
    ax.text(LBL_X, r.y, f"{r.sig}  ({r.ct})", ha="right", va="center", fontsize=10)

# group labels sit just right of the dots (compact)
GX = 0.17
ax.text(GX, df[df.group == "glial"].y.mean(), "Glial /\nprogenitor", ha="left", va="center",
        fontsize=9, fontweight="bold", color="#555", linespacing=0.95)
ax.text(GX, df[df.group == "immune"].y.mean(), "Immune /\nmyeloid", ha="left", va="center",
        fontsize=9, fontweight="bold", color="#555", linespacing=0.95)
ax.axhline((df[df.group == "glial"].y.min() + df[df.group == "immune"].y.max()) / 2,
           color="0.82", lw=0.8, ls=(0, (4, 3)), zorder=0)

ax.set_xlim(-1.75, 0.55)
ax.set_ylim(df.y.min() - 0.7, df.y.max() + 0.7)

# horizontal NES scale below the dot column
cax = fig.add_axes([0.52, 0.11, 0.34, 0.040])
cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
cb.set_label("NES (BCAAr vs control)", fontsize=9)
cb.ax.tick_params(labelsize=8)

for ext in ("png", "svg", "pdf"):
    fig.savefig(f"{FIG}/Braindex_depleted_dotplot_2HIT.{ext}",
                dpi=200 if ext == "png" else None, bbox_inches="tight")
plt.close(fig)
print(df[["sig", "ct", "group", "NES"]].to_string(index=False))
print("saved", f"{FIG}/Braindex_depleted_dotplot_2HIT.png / .svg / .pdf")
