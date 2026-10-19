#!/usr/bin/env python3
"""
Stacked panel figure: ranked MAGeCK-MLE beta scores (top) + KEGG GSEA NES bars (bottom)
for the DIPGXIII metabolic CRISPR screen, K27M vs K27M-KO.

Rebuild of gsea_out/stacked_dotplot_NES_K27M_KO.png with:
  - BCAT2 added alongside BCAT1 in both rank panels
  - BCAT1 / BCAT2 colored from BCAT1_figure_palette.json -> "genes"
    (ColorBrewer Dark2 pair, #1b9e77 / #7570b3), matching human_data/figB6_BCAT_simple,
    tabula_muris_OPC_oligo/figM_mouse_BCAT_simple and BCAT1_DepMap_analysis/fig_bcat1_vs_bcat2
  - bar edge strokes removed
  - vector SVG export with editable text for Illustrator

Inputs : ../beta_correlation.xlsx, gsea_K27M.csv, gsea_KO.csv
Outputs: stacked_dotplot_NES_K27M_KO_v2.svg / .png / .pdf
"""

import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# --- keep text as text in the SVG/PDF (Illustrator-editable, not outlined) -----
mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["font.family"] = "sans-serif"
mpl.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]

HERE = Path(__file__).resolve().parent
SCREEN = HERE.parent
ROOT = SCREEN.parent

# ----------------------------------------------------------------------------
# Colors
# ----------------------------------------------------------------------------
with open(ROOT / "BCAT1_figure_palette.json") as fh:
    PAL = json.load(fh)

# project-canonical per-gene colors (ColorBrewer Dark2): BCAT1 #1b9e77, BCAT2 #7570b3
GENE_COLOR = {g: PAL["genes"][g]["hex"] for g in ("BCAT1", "BCAT2")}
TERPENOID_C = "#c0392b"   # mevalonate / terpenoid backbone leading edge
FATTYACID_C = "#2166ac"   # fatty acid biosynthesis leading edge
GREY = "#b3b3b3"

# ----------------------------------------------------------------------------
# Data
# ----------------------------------------------------------------------------
betas = pd.read_excel(SCREEN / "beta_correlation.xlsx")
gsea = {
    "K27M": pd.read_csv(HERE / "gsea_K27M.csv").set_index("Term"),
    "KO": pd.read_csv(HERE / "gsea_KO.csv").set_index("Term"),
}

PANELS = [("K27M", "K27M.beta"), ("K27M-KO", "KO.beta")]
GSEA_KEY = {"K27M": "K27M", "K27M-KO": "KO"}

# genes annotated on the K27M rank panel (leading edge of the two significant sets)
TERPENOID = ["HMGCR", "FDPS", "ACAT2", "DHDDS", "MVD"]
FATTYACID = ["ACSL1", "HSD17B8", "MECR", "ACACA", "OXSM"]
BCATS = ["BCAT1", "BCAT2"]

TERMS = [
    "Terpenoid backbone biosynthesis",
    "Fatty acid biosynthesis",
    "Glycolysis / Gluconeogenesis",
    "Inositol phosphate metabolism",
    "Thiamine metabolism",
    "Pentose and glucuronate interconversions",
]

FDR_SIG = 0.05

# ----------------------------------------------------------------------------
# Figure
# ----------------------------------------------------------------------------
fig = plt.figure(figsize=(11.3, 9.4))
gs = fig.add_gridspec(
    2, 2, height_ratios=[1.0, 0.95], hspace=0.28, wspace=0.18,
    left=0.235, right=0.965, top=0.955, bottom=0.075,
)


def rank_panel(ax, label, col):
    d = betas[["Gene", col]].dropna().sort_values(col).reset_index(drop=True)
    d["rank"] = np.arange(1, len(d) + 1)
    pos = d.set_index("Gene")

    ax.axhline(0, color="black", lw=0.9, zorder=1)
    ax.scatter(d["rank"], d[col], s=11, color=GREY, linewidths=0, zorder=2)

    ax.set_xlim(-60, len(d) + 60)
    ax.set_ylim(-2, 2)
    ax.set_xlabel("Gene rank", fontsize=12)
    ax.set_title(label, fontsize=15, fontweight="bold", pad=10)
    ax.tick_params(labelsize=10)

    def mark(gene, color, marker, size):
        if gene not in pos.index:
            return None
        x, y = pos.loc[gene, "rank"], pos.loc[gene, col]
        ax.scatter([x], [y], s=size, color=color, marker=marker,
                   edgecolor="black", linewidths=0.8, zorder=5)
        return x, y

    def direct_label(gene, dx, dy, ha="left"):
        """Marker + short leader line to a label placed near the point."""
        xy = mark(gene, GENE_COLOR[gene], "D", 95)
        if xy is None:
            return
        x, y = xy
        ax.annotate(
            gene, xy=(x, y), xytext=(x + dx, y + dy), color=GENE_COLOR[gene],
            fontsize=13, fontweight="bold", va="center", ha=ha,
            arrowprops=dict(arrowstyle="-", color=GENE_COLOR[gene], lw=1.1,
                            shrinkA=0, shrinkB=4),
            zorder=6,
        )

    if label == "K27M":
        # ---- pathway leading-edge genes, leader lines fanned to the right ----
        entries = [(g, TERPENOID_C, "o", 58) for g in TERPENOID]
        entries += [(g, FATTYACID_C, "o", 58) for g in FATTYACID]
        entries += [("BCAT1", GENE_COLOR["BCAT1"], "D", 95)]

        placed = []
        for g, c, m, s in entries:
            xy = mark(g, c, m, s)
            if xy:
                placed.append((g, c, xy))

        # order labels by y so leader lines don't cross
        placed.sort(key=lambda t: -t[2][1])
        y_slots = np.linspace(-0.36, -1.60, len(placed))
        x_short, x_long = 620, 1500
        for i, ((g, c, (x, y)), ylab) in enumerate(zip(placed, y_slots)):
            xlab = x_short if i % 2 == 0 else x_long
            ax.annotate(
                g, xy=(x, y), xytext=(xlab, ylab), color=c,
                fontsize=12.5, fontweight="bold", va="center", ha="left",
                arrowprops=dict(arrowstyle="-", color=c, lw=1.0,
                                shrinkA=0, shrinkB=3,
                                connectionstyle="angle,angleA=0,angleB=90,rad=0"),
                zorder=6,
            )
        # BCAT2 is near-neutral and well clear of the cluster -> label directly
        direct_label("BCAT2", 260, 0.45)
    else:
        # ---- KO panel: BCAT1 / BCAT2 only ----
        direct_label("BCAT1", -260, 0.45, ha="right")
        direct_label("BCAT2", 260, 0.45)


def nes_panel(ax, label, show_ylabels):
    g = gsea[GSEA_KEY[label]]
    nes = np.array([g.loc[t, "NES"] for t in TERMS])
    fdr = np.array([g.loc[t, "FDR q-val"] for t in TERMS])
    y = np.arange(len(TERMS))

    colors = []
    for v, q in zip(nes, fdr):
        base = "#c0392b" if v > 0 else "#2166ac"
        colors.append(base if q < FDR_SIG else _lighten(base, 0.70))

    # linewidth=0 -> no bar edge stroke
    ax.barh(y, nes, height=0.62, color=colors, linewidth=0, zorder=3)
    ax.axvline(0, color="black", lw=1.1, zorder=4)

    for yi, v, q in zip(y, nes, fdr):
        sig = q < FDR_SIG
        txt = f"FDR={q:.3f}"
        if v < 0:
            ax.text(0.12, yi, txt, va="center", ha="left", fontsize=10.5,
                    color="black" if sig else "#8c8c8c",
                    fontweight="bold" if sig else "normal", zorder=5)
        else:
            ax.text(-0.12, yi, txt, va="center", ha="right", fontsize=10.5,
                    color="black" if sig else "#8c8c8c",
                    fontweight="bold" if sig else "normal", zorder=5)

    ax.set_yticks(y)
    ax.set_yticklabels(TERMS if show_ylabels else [], fontsize=12)
    ax.set_ylim(-0.65, len(TERMS) - 0.35)
    ax.set_xlim(-3, 3)
    ax.set_xlabel("NES", fontsize=12)
    ax.tick_params(labelsize=10, left=False)


def _lighten(hexcolor, amount):
    """Blend a hex color toward white by `amount` (0 = unchanged, 1 = white)."""
    rgb = np.array(mpl.colors.to_rgb(hexcolor))
    return mpl.colors.to_hex(rgb + (1.0 - rgb) * amount)


for j, (label, col) in enumerate(PANELS):
    rank_panel(fig.add_subplot(gs[0, j]), label, col)
    nes_panel(fig.add_subplot(gs[1, j]), label, show_ylabels=(j == 0))

fig.axes[0].set_ylabel("beta score", fontsize=12)

for ext in ("svg", "png", "pdf"):
    out = HERE / f"stacked_dotplot_NES_K27M_KO_v2.{ext}"
    fig.savefig(out, dpi=300, bbox_inches="tight",
                facecolor="white", transparent=False)
    print("wrote", out)
