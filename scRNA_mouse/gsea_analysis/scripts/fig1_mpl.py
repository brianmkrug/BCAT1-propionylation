#!/usr/bin/env python3
"""fig1 cell-identity UMAP (v2, Unresolved dropped) rebuilt in matplotlib for
Illustrator-editable vector export.

Font handling (the whole point of this rebuild):
  svg.fonttype = 'none'  -> SVG keeps real <text> elements (editable strings, not
                            outlined paths).
  pdf.fonttype = 42      -> PDF embeds TrueType as CID fonts; pdffonts should report
                            DejaVuSans / DejaVuSans-Bold, emb=yes, uni=yes, so text
                            opens as editable type boxes in Illustrator.
  font.family = DejaVu Sans -> the bundled matplotlib face, always available/embeddable.
"""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "..", ".."))  # publication: repo-relative root
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
import pandas as pd, numpy as np

plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42          # TrueType (Type42/CID), not Type3 outlines
plt.rcParams["ps.fonttype"]  = 42
plt.rcParams["font.family"]  = "DejaVu Sans"

RES = REPO_ROOT + "/scRNA_mouse/gsea_analysis/results"
FIG = REPO_ROOT + "/scRNA_mouse/gsea_analysis/figures"

df = pd.read_csv(f"{RES}/umap_scores.csv")
mp = {"Malignant_OPC_like":"OPC","OPC":"OPC","Malignant_cycling":"Cycling",
      "Oligodendrocyte":"Oligodendrocyte","Oligodendrocyte_NFOL":"Oligodendrocyte",
      "Astrocyte":"Astrocyte","Neuron":"Neuronal",
      "Microglia":"Immune","Myeloid":"Immune","T_cell":"Immune",
      "Endothelial":"Vascular/stromal","VLMC_fibroblast":"Vascular/stromal",
      "Unresolved_IFN_high":"Unresolved"}
df["ct"]  = df["cell_type"].map(mp)
df = df[df["ct"] != "Unresolved"].copy()   # v2: drop Unresolved entirely
df = df.sample(frac=1, random_state=42)    # shuffle so no type paints over another

LEV = ["OPC","Cycling","Oligodendrocyte","Astrocyte","Neuronal","Immune","Vascular/stromal"]
COL = {"OPC":"#009E73","Cycling":"#785EF0","Oligodendrocyte":"#D55E00","Astrocyte":"#56B4E9",
       "Neuronal":"#0072B2","Immune":"#CC79A7","Vascular/stromal":"#E69F00"}
LAB = {"OPC":"#008560","Cycling":"#7257F0","Oligodendrocyte":"#BE5400","Astrocyte":"#377BA2",
       "Neuronal":"#0072B2","Immune":"#A35D84","Vascular/stromal":"#9A6A00"}
SURF = "#fcfcfb"
n = df["ct"].value_counts()

def density_peak(sub, bins=60):
    """Anchor a label in the visually densest region of a category (the peak bin
    of a 2D histogram), avoiding the empty gap a plain median can fall into when a
    cluster has two lobes (e.g. oligodendrocyte)."""
    H, xe, ye = np.histogram2d(sub["UMAP_1"], sub["UMAP_2"], bins=bins)
    i, j = np.unravel_index(np.argmax(H), H.shape)
    return (xe[i] + xe[i+1]) / 2, (ye[j] + ye[j+1]) / 2

fig, ax = plt.subplots(figsize=(12, 9.2))
fig.patch.set_facecolor(SURF); ax.set_facecolor(SURF)
fig.subplots_adjust(left=0.02, right=0.80, top=0.88, bottom=0.06)   # room for header/caption

for ctype in LEV:                          # fixed order -> stable colour assignment
    sub = df[df["ct"] == ctype]
    ax.scatter(sub["UMAP_1"], sub["UMAP_2"], s=5.5, c=COL[ctype],   # larger dots for clarity
               linewidths=0, alpha=0.65, rasterized=True,
               label=f"{ctype} ({n[ctype]:,})")

for ctype in LEV:
    x, y = density_peak(df[df["ct"] == ctype])
    # White rounded bbox instead of a stroke halo: path_effects withStroke OUTLINES
    # the text to a <path> (not editable); a bbox keeps the label a real <text>
    # element and renders the halo as a separate, deletable <rect>.
    ax.text(x, y, ctype, color=LAB[ctype], fontsize=13, fontweight="bold",
            ha="center", va="center", zorder=10,
            bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.72))

ax.set_aspect("equal"); ax.axis("off")

# Header/caption as figure-level text (distinct y bands, above the axes -> no overlap)
fig.text(0.02, 0.965, "KPP/KNF mouse H3K27M glioma - cell identity",
         fontsize=15, fontweight="bold", color="#1a1a1a", va="top")
fig.text(0.02, 0.935,
         "73,151 nuclei | harmony-integrated UMAP | 5 samples | labelled by transcriptional cell type\n"
         "Cells named by best-matching identity (no tumour/normal call - GFP too ambient); "
         "OPC includes formerly tumour-flagged OPC-like nuclei",
         fontsize=9.5, color="#555555", va="top", linespacing=1.4)
fig.text(0.02, 0.028,
         "v2: 4,909 IFN-high 'Unresolved' nuclei (clusters 6/18) dropped. "
         "Cycling = proliferative state (best lineage match OPC).",
         fontsize=8.5, color="#777777", va="bottom")

leg = ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False,
                handletextpad=0.5, labelspacing=0.9, fontsize=10, markerscale=2.6)
for t in leg.get_texts(): t.set_color("#333333")

# No bbox_inches='tight' -> keep the manual subplots_adjust / fig.text layout intact.
fig.savefig(f"{FIG}/fig1_cell_identity_umap_v2.png", dpi=300, facecolor=SURF)
fig.savefig(f"{FIG}/fig1_cell_identity_umap_v2.pdf", facecolor=SURF)
fig.savefig(f"{FIG}/fig1_cell_identity_umap_v2.svg", facecolor=SURF)
print("wrote png/pdf/svg; cells per type:")
print(n[LEV].to_string())
