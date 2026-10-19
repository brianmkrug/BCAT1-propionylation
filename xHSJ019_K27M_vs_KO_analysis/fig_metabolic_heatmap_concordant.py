#!/usr/bin/env python3
"""Concordant K27M metabolic program heatmap (HSJ019 & BT245), rebuilt in the
standardized BCAT1 figure style.

Layout (matches BCAAr_analysis/05_2HIT_LE_heatmap.py conventions):
  * category colour strip on the far left, then gene names, then the Z-score heatmap;
  * top annotation panel colours each sample column by condition (K27M / KO) with the
    cell-line (HSJ019 / BT245) labelled above it;
  * thin WHITE gaps separate gene blocks (horizontal) and sample groups (vertical);
  * no plot title, no bottom sample names;
  * the log2FC heatmap is placed to the RIGHT of the Z-score heatmap, with its own
    colour scale.
All colours are loaded from BCAT1_figure_palette.json (single source of truth).
"""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))  # publication: repo-relative root
import json
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["svg.fonttype"] = "none"   # keep text editable in Illustrator
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.cm as cm
from matplotlib.colors import Normalize
from pydeseq2.dds import DeseqDataSet

ROOT = REPO_ROOT + ""
HSJ = f"{ROOT}/xHSJ019_K27M_vs_KO_analysis"
BT = f"{ROOT}/xBT245_K27M_vs_KO_analysis"
FIG = f"{HSJ}/figures/fig_metabolic_heatmap_HSJ019_BT245_concordant.png"

with open(f"{ROOT}/BCAT1_figure_palette.json") as f:
    PAL = json.load(f)
CAT_COLORS = {k: v["hex"] for k, v in PAL["metabolic_categories"].items()}
COND_COLORS = {k: PAL["conditions"][k]["hex"] for k in ("K27M", "KO")}
HM = PAL["heatmap"]; LHM = PAL["log2fc_heatmap"]
GAP = HM["block_separator"]
DIV = PAL["annotations"]["subgroup_divider"]          # within-block dashed divider
HL = PAL["annotations"]["highlight_gene"]["hex"]      # emphasized gene-label color
HIGHLIGHT = "BCAT1"   # gene to emphasise

# ---- genes / groups / log2FC (concordant table) ----
tbl = pd.read_csv(f"{HSJ}/figures/fig_metabolic_heatmap_HSJ019_BT245_concordant_table.csv",
                  index_col=0)
genes = list(tbl.index)                                   # display order (as built originally)
groups = list(dict.fromkeys(tbl["group"]))                # block order
lfc = tbl[["lfc_HSJ019", "lfc_BT245"]].copy()

# ---- per-line VST, then z-score within line ----
LINES = {
    "HSJ019": (f"{HSJ}/gene_counts_raw.csv",
               ["K27M-167", "K27M-172", "K27M-616"], ["KO-613", "KO-614", "KO-618"]),
    "BT245":  (f"{BT}/gene_counts_raw.csv",
               ["K27M_663A", "K27M_663B"], ["KO_619", "KO_646", "KO_650"]),
}

def vst_z(path, k27m, ko):
    c = pd.read_csv(path, index_col=0)
    cols = k27m + ko
    c = c[cols]
    c = c[c.sum(1) >= 10]
    meta = pd.DataFrame({"condition": ["K27M"] * len(k27m) + ["KO"] * len(ko)}, index=cols)
    dds = DeseqDataSet(counts=c.T, metadata=meta, design="~condition", quiet=True)
    dds.deseq2(); dds.vst(use_design=False)
    v = pd.DataFrame(dds.layers["vst_counts"], index=cols, columns=c.index).T  # genes x samples
    v = v.reindex(genes)
    z = v.sub(v.mean(1), axis=0).div(v.std(1, ddof=0).replace(0, np.nan), axis=0)
    return z, cols

zmats, order_cols, cond_seq, line_spans = [], [], [], []
col0 = 0
for line, (path, k27m, ko) in LINES.items():
    z, cols = vst_z(path, k27m, ko)
    zmats.append(z.values)
    order_cols += [f"{line}:{s}" for s in cols]
    cond_seq += ["K27M"] * len(k27m) + ["KO"] * len(ko)
    line_spans.append((line, col0, col0 + len(cols)))
    col0 += len(cols)
Z = np.hstack(zmats)                       # genes x all-samples
ncols = Z.shape[1]
nrows = len(genes)

# ---- block boundaries (row) ----
row_bounds = []            # (group, start, end)
s = 0
for g in groups:
    n = int((tbl["group"] == g).sum())
    row_bounds.append((g, s, s + n)); s += n
# up/down split inside the last block (down-regulated genes at the bottom)
down_start = next((i for i, gname in enumerate(genes) if lfc.loc[gname, "lfc_HSJ019"] < 0), None)

# ---- vertical sample-group gaps ----
cond_gaps = [i for i in range(1, ncols) if cond_seq[i] != cond_seq[i - 1]]  # K27M|KO
line_gaps = [span[1] for span in line_spans[1:]]                             # HSJ|BT

# =====================================================================
fig = plt.figure(figsize=(11.5, 13))
# axis rectangles [left, bottom, width, height]
b, h = 0.055, 0.80
ax_cat = fig.add_axes([0.150, b, 0.018, h])                 # category colour strip
ax_z   = fig.add_axes([0.250, b, 0.400, h])                 # z-score heatmap
ax_l   = fig.add_axes([0.690, b, 0.085, h])                 # log2FC heatmap
cax_z  = fig.add_axes([0.845, 0.36, 0.016, 0.30])           # z-score colourbar
cax_l  = fig.add_axes([0.940, 0.36, 0.016, 0.30])           # log2FC colourbar

# ---------- Z-score heatmap ----------
im = ax_z.imshow(Z, cmap=HM["colormap"], vmin=HM["vmin"], vmax=HM["vmax"], aspect="auto")
ax_z.set_xticks([])
ax_z.set_yticks(range(nrows))
ax_z.set_yticklabels(genes, fontsize=7)
ax_z.tick_params(axis="y", length=0, pad=2)
for lbl in ax_z.get_yticklabels():
    if lbl.get_text() == HIGHLIGHT:
        lbl.set_color(HL); lbl.set_fontweight("bold")
for sp in ax_z.spines.values():
    sp.set_visible(False)

# vertical white gaps
for x in cond_gaps:
    ax_z.axvline(x - 0.5, color=GAP["color"], lw=GAP["linewidth"], zorder=5)
for x in line_gaps:
    ax_z.axvline(x - 0.5, color=GAP["color"], lw=GAP["linewidth"] + 1.4, zorder=6)
# horizontal white gaps between blocks
for g, s0, e0 in row_bounds[1:]:
    ax_z.axhline(s0 - 0.5, color=GAP["color"], lw=GAP["linewidth"], zorder=5)
if down_start is not None:
    ax_z.axhline(down_start - 0.5, color=DIV["color"], lw=DIV["linewidth"],
                 ls=(0, (4, 3)), zorder=5)

# ---------- top annotation: condition panel + cell-line labels ----------
bar_h = 0.9
bar_top = -0.5 - bar_h
# condition-coloured segments (respect white gaps at condition/line breaks)
breaks = sorted(set([0, ncols] + cond_gaps + line_gaps))
for a, bnd in zip(breaks[:-1], breaks[1:]):
    ax_z.add_patch(Rectangle((a - 0.5, bar_top), bnd - a, bar_h,
                             color=COND_COLORS[cond_seq[a]], clip_on=False, lw=0, zorder=4))
    ax_z.text((a + bnd - 1) / 2, bar_top + bar_h / 2, cond_seq[a], ha="center", va="center",
              fontsize=8, fontweight="bold", color="white", clip_on=False, zorder=5)
# cell-line labels above the panel
line_lbl_y = bar_top - 0.9
for line, s0, e0 in line_spans:
    ax_z.text((s0 + e0 - 1) / 2, line_lbl_y, line, ha="center", va="center",
              fontsize=12, fontweight="bold", clip_on=False)
ax_z.set_xlim(-0.5, ncols - 0.5)
ax_z.set_ylim(nrows - 0.5, line_lbl_y - 0.6)

# ---------- category colour strip (left) + vertical labels ----------
ax_cat.set_xlim(0, 1); ax_cat.set_ylim(nrows - 0.5, line_lbl_y - 0.6)
ax_cat.axis("off")
for g, s0, e0 in row_bounds:
    ax_cat.add_patch(Rectangle((0.55, s0 - 0.5), 0.45, e0 - s0,
                               color=CAT_COLORS[g], clip_on=False, lw=0))
    ax_cat.text(0.35, (s0 + e0) / 2 - 0.5, g, rotation=90, va="center", ha="center",
                fontsize=9, fontweight="bold", color=CAT_COLORS[g], clip_on=False)

# ---------- log2FC heatmap (right) ----------
Lm = lfc.values                                    # genes x 2  (HSJ, BT)
lmax = np.nanmax(np.abs(Lm))
imL = ax_l.imshow(Lm, cmap=LHM["colormap"], vmin=-lmax, vmax=lmax, aspect="auto")
ax_l.set_yticks([]); ax_l.set_xticks([0, 1])
ax_l.set_xticklabels(["HSJ", "BT"], fontsize=8)
ax_l.tick_params(axis="x", length=0, pad=2)
ax_l.xaxis.set_ticks_position("top")
for sp in ax_l.spines.values():
    sp.set_visible(False)
for g, s0, e0 in row_bounds[1:]:
    ax_l.axhline(s0 - 0.5, color=GAP["color"], lw=GAP["linewidth"], zorder=5)
if down_start is not None:
    ax_l.axhline(down_start - 0.5, color=DIV["color"], lw=DIV["linewidth"],
                 ls=(0, (4, 3)), zorder=5)
ax_l.set_ylim(nrows - 0.5, line_lbl_y - 0.6)
ax_l.text(0.5, line_lbl_y - 0.2, "log2FC", ha="center", va="center",
          fontsize=10, fontweight="bold", clip_on=False)

# ---------- colourbars ----------
cb = fig.colorbar(im, cax=cax_z); cb.set_label(HM["value_label"], fontsize=9)
cbL = fig.colorbar(imL, cax=cax_l); cbL.set_label(LHM["value_label"], fontsize=9)

# ---------- source-data reference table (exact values plotted) ----------
src = pd.DataFrame(Z, index=genes, columns=[c.replace(":", "_") for c in order_cols])
src.insert(0, "category", tbl["group"].values)
src["lfc_HSJ019"] = lfc["lfc_HSJ019"].values
src["lfc_BT245"] = lfc["lfc_BT245"].values
src["padj_HSJ019"] = tbl["padj_HSJ019"].values
src["padj_BT245"] = tbl["padj_BT245"].values
src.index.name = "gene"
SRC = f"{HSJ}/figures/fig_metabolic_heatmap_HSJ019_BT245_concordant_SOURCEDATA.csv"
src.to_csv(SRC)

# ---------- save figure (PNG raster + SVG for Illustrator) ----------
SVG = FIG.replace(".png", ".svg")
for path in (FIG, f"{BT}/figures/fig_metabolic_heatmap_HSJ019_BT245_concordant.png"):
    fig.savefig(path, dpi=170, bbox_inches="tight")
fig.savefig(SVG, bbox_inches="tight")                       # vector, editable in Illustrator
fig.savefig(f"{BT}/figures/fig_metabolic_heatmap_HSJ019_BT245_concordant.svg",
            bbox_inches="tight")
plt.close(fig)
print("saved", FIG)
print("saved", SVG)
print("saved", SRC)
