#!/usr/bin/env python3
"""fig4 v2 (fatty acid metabolism) rebuilt in matplotlib for Illustrator-editable
vector export. Changes requested: violin plots instead of boxplots; no plot-area
gridlines. Unresolved cells already dropped (v2). Z recomputed on retained nuclei.

Font handling (editable vector text):
  svg.fonttype='none'  -> SVG keeps real <text> elements.
  pdf.fonttype=42      -> PDF embeds DejaVuSans/-Bold as CID TrueType (emb=yes,uni=yes).
"""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "..", ".."))  # publication: repo-relative root
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.gridspec import GridSpec
import pandas as pd, numpy as np

plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"]  = 42
plt.rcParams["font.family"]  = "DejaVu Sans"

RES = REPO_ROOT + "/scRNA_mouse/gsea_analysis/results"
FIG = REPO_ROOT + "/scRNA_mouse/gsea_analysis/figures"
PATH, NICE, NGENES, PCT0 = "Fatty_acid_metabolism", "Fatty acid metabolism", 152, 0.2

df = pd.read_csv(f"{RES}/umap_scores.csv")
mp = {"Malignant_OPC_like":"OPC","OPC":"OPC","Malignant_cycling":"Cycling",
      "Oligodendrocyte":"Oligodendrocyte","Oligodendrocyte_NFOL":"Oligodendrocyte",
      "Astrocyte":"Astrocyte","Neuron":"Neuronal","Microglia":"Immune","Myeloid":"Immune",
      "T_cell":"Immune","Endothelial":"Vascular/stromal","VLMC_fibroblast":"Vascular/stromal",
      "Unresolved_IFN_high":"Unresolved"}
df["ct"] = df["cell_type"].map(mp)
df = df[df["ct"] != "Unresolved"].copy()                 # v2
df["Z"] = (df[PATH] - df[PATH].mean()) / df[PATH].std()  # Z across retained nuclei

COL = {"OPC":"#009E73","Cycling":"#785EF0","Oligodendrocyte":"#D55E00","Astrocyte":"#56B4E9",
       "Neuronal":"#0072B2","Immune":"#CC79A7","Vascular/stromal":"#E69F00"}
SURF = "#fcfcfb"
order = df.groupby("ct")["Z"].mean().sort_values(ascending=False).index.tolist()

fig = plt.figure(figsize=(12, 11))
fig.patch.set_facecolor(SURF)
gs = GridSpec(2, 1, height_ratios=[1.45, 1], hspace=0.16,
              left=0.07, right=0.92, top=0.93, bottom=0.09)

# ---- top: hex-binned mean-Z UMAP (diverging) ----
axU = fig.add_subplot(gs[0]); axU.set_facecolor(SURF)
cmap = LinearSegmentedColormap.from_list("rbdiv", ["#2166AC", "#F0F0F0", "#B2182B"])
hb = axU.hexbin(df["UMAP_1"], df["UMAP_2"], C=df["Z"], reduce_C_function=np.mean,
                gridsize=90, cmap=cmap, vmin=-1, vmax=1, mincnt=1, linewidths=0)
axU.set_aspect("equal"); axU.axis("off"); axU.set_anchor("W")   # left-align UMAP under title
cb = fig.colorbar(hb, ax=axU, fraction=0.03, pad=0.01, shrink=0.55)
cb.set_label("mean Z", fontsize=9); cb.ax.tick_params(labelsize=8); cb.outline.set_visible(False)
# Panel headers as figure-level text: equal-aspect shrinks axU, so axes-anchored
# titles drift right and clip. Figure coords keep both headers at the left margin.
fig.text(0.055, 0.955, f"{NICE} - Z score of per-nucleus UCell",
         fontsize=13, fontweight="bold", color="#1a1a1a", va="top")
fig.text(0.055, 0.935,
         f"diverging ramp, neutral grey midpoint | hex-binned mean, clipped to +/-1 SD | "
         f"{PCT0}% score 0 ({NGENES} genes) | Unresolved dropped",
         fontsize=9, color="#555555", va="top")

# ---- bottom: violins by cell type ----
axV = fig.add_subplot(gs[1]); axV.set_facecolor(SURF)
data = [df.loc[df["ct"] == c, "Z"].values for c in order]
pos  = np.arange(1, len(order) + 1)
axV.axhline(0, color="#cccccc", lw=0.8, zorder=1)
parts = axV.violinplot(data, positions=pos, widths=0.82, showextrema=False,
                       showmedians=True, points=200)
for body, c in zip(parts["bodies"], order):
    body.set_facecolor(COL[c]); body.set_edgecolor("#444444")
    body.set_linewidth(0.4); body.set_alpha(0.9); body.set_zorder(2)
cm = parts["cmedians"]; cm.set_color("#222222"); cm.set_linewidth(1.1); cm.set_zorder(3)
means = [d.mean() for d in data]
axV.scatter(pos, means, marker="D", s=42, facecolor="white", edgecolor="#111111",
            linewidth=0.7, zorder=4)

axV.set_xticks(pos)
axV.set_xticklabels([f"{c}\n(n={len(df[df['ct']==c]):,})" for c in order], fontsize=8.5, color="#333333")
axV.set_ylim(-2.6, 3.2); axV.set_ylabel("Z score", fontsize=11)
axV.tick_params(axis="y", labelsize=9)
axV.grid(False)                                          # no plot-area gridlines
for s in ("top", "right"): axV.spines[s].set_visible(False)
for s in ("left", "bottom"): axV.spines[s].set_color("#999999")
axV_top = axV.get_position().y1
fig.text(0.055, axV_top + 0.030, f"{NICE} score by cell type",
         fontsize=13, fontweight="bold", color="#1a1a1a", va="top")
fig.text(0.055, axV_top + 0.012,
         "per-nucleus UCell, Z-scored across 73,151 retained nuclei | violin = full distribution | "
         "black line = median | white diamond = mean | ordered by mean",
         fontsize=8.5, color="#555555", va="top")

for ext in ("png", "pdf", "svg"):
    fig.savefig(f"{FIG}/fig4_fattyacid_zscore_by_celltype_v2.{ext}",
                dpi=300 if ext == "png" else None, facecolor=SURF)
print("wrote png/pdf/svg. Z by cell type (ordered):")
for c in order:
    z = df.loc[df["ct"] == c, "Z"]
    print(f"  {c:<18} n={len(z):>6}  mean={z.mean():+.3f}  median={z.median():+.3f}")
