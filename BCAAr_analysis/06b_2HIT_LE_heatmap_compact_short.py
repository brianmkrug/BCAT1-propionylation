#!/usr/bin/env python3
"""2HIT leading-edge heatmap — COMPACT variant (version 2).

Same gene selection and data as 05_2HIT_LE_heatmap.py, but reformatted for a small
footprint:
  * dramatically thinner rows;
  * per-gene names removed;
  * only 2-3 leading-edge genes per category are called out with a leader line;
  * the propionyl-CoA-upstream set (Hsd17b10, Hadha, Bckdha, Acads, Hibch) is bold and
    bracketed as "upstream of propionyl-CoA handling";
  * OPC block calls out Olig1 / Olig2 / Cspg4.
Colours come from BCAT1_figure_palette.json. Saves PNG + SVG (editable text).
"""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))  # publication: repo-relative root
import json
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["svg.fonttype"] = "none"   # SVG: real <text> elements, editable
matplotlib.rcParams["pdf.fonttype"] = 42       # PDF: embed CID TrueType (editable type)
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.transforms import blended_transform_factory
import re as _re
from pydeseq2.dds import DeseqDataSet

ROOT = REPO_ROOT + ""
OUT = f"{ROOT}/BCAAr_analysis/out"
FIG = f"{ROOT}/BCAAr_analysis/figures"

with open(f"{ROOT}/BCAT1_figure_palette.json") as _f:
    PAL = json.load(_f)
sig_colors = {k: v["hex"] for k, v in PAL["pathways"].items()}
NICE = {k: v["label"] for k, v in PAL["pathways"].items()}
COND_COLORS = {k: v["hex"] for k, v in PAL["conditions"].items()}
HM = PAL["heatmap"]; GAP = HM["block_separator"]

# ---------------------------------------------------------------- data / selection
counts = pd.read_csv(f"{OUT}/counts_gene.csv", index_col=0)
meta = pd.read_csv(f"{OUT}/sample_metadata.csv", index_col=0)
m2 = meta[meta.experiment == "2HIT"].copy()
order = (list(m2[m2.condition == "control"].sort_values("rep").index) +
         list(m2[m2.condition == "BCAAr"].sort_values("rep").index))
m2 = m2.loc[order]
c2 = counts[order]; c2 = c2[c2.sum(axis=1) >= 10]

dds = DeseqDataSet(counts=c2.T, metadata=m2[["condition"]], design="~condition", quiet=True)
dds.deseq2(); dds.vst(use_design=False)
vst = pd.DataFrame(dds.layers["vst_counts"], index=order, columns=c2.index)

gsea = pd.read_csv(f"{OUT}/GSEA_all_results_long.csv")
gmet = gsea[(gsea.experiment == "2HIT") &
            (gsea.collection == "Curated_metabolic_panel")].drop_duplicates("Term").set_index("Term")
def lead(term):
    return [g for g in str(gmet.loc[term, "Lead_genes"]).split(";") if g in vst.columns]

rnk = pd.read_csv(f"{OUT}/rank_2HIT.rnk", sep="\t", header=None, index_col=0)[1]
rnk = rnk[~rnk.index.duplicated()]
de = pd.read_csv(f"{OUT}/deseq2_2HIT_BCAArVScontrol.csv", index_col=0)
padj = de["padj"]; l2fc = de["log2FoldChange"]

SIGS = list(PAL["pathways"].keys())
CAPS = {"Glycolysis": 11, "OXPHOS": 15}
FIRST_IDENTITY = next(k for k, v in PAL["pathways"].items() if v.get("section") == "cell_identity")

def padj_of(g):
    p = padj.get(g, 1.0)
    return p if pd.notna(p) else 1.0

groups, seen = [], set()
DROP = {"Rlbp1", "Cpt1a", "Hacd4", "Eci1", "Decr1", "Acad8", "Aldh7a1", "Acss1",
        "mt-Nd1", "Atp6v1c2", "Gck", "Aldh3b2", "Eno4"}
def add(sig, gene):
    if gene not in seen and gene not in DROP:
        groups.append((sig, gene)); seen.add(gene)

for t in ["Glycolysis", "OXPHOS"]:
    for g in sorted(lead(t), key=padj_of)[:CAPS[t]]:
        add(t, g)
for g in sorted(lead("Propanoate_metab"), key=lambda g: -abs(rnk.get(g, 0))):
    if g != "Bckdha":
        add("Propanoate_metab", g)
FORCE_BCAA = ["Bckdha"]
bcaa_le = [g for g in lead("BCAA_degradation") if g not in seen]
pool = [g for g in bcaa_le if g not in FORCE_BCAA]
pool_down = [g for g in pool if l2fc.get(g, 0) < 0]
bcaa_keep = sorted(pool_down, key=padj_of)[:max(1, round(len(pool) / 2))]
for f in FORCE_BCAA:
    if f in bcaa_le and f not in bcaa_keep:
        bcaa_keep.append(f)
for g in sorted(bcaa_keep, key=padj_of):
    add("BCAA_degradation", g)

LIPID_SIGS = ["Unsat_FA_synthesis", "FA_elongation", "FA_oxidation",
              "HALLMARK_Fatty_acid_metab", "HALLMARK_Adipogenesis"]
lipid_membership = {}
for t in LIPID_SIGS:
    if t in gmet.index:
        for g in str(gmet.loc[t, "Lead_genes"]).split(";"):
            if g in vst.columns:
                lipid_membership.setdefault(g, set()).add(t)
lipid_genes = [g for g, s in lipid_membership.items()
               if g not in seen and (len(s) >= 2 or padj_of(g) < 0.25)]
lipid_genes = sorted(lipid_genes, key=lambda g: (-len(lipid_membership[g]), padj_of(g)))
for g in lipid_genes:
    add("Lipid_FA", g)

bx = pd.read_csv(f"{ROOT}/BCAAr_analysis/gsea/Braindex_mouse_2HIT.csv").set_index("Term")
opc_terms = [t for t in bx.index if _re.search(r"OPC", t)
             and pd.to_numeric(bx.loc[t, "NES"]) < 0 and pd.to_numeric(bx.loc[t, "FDR q-val"]) < 0.25]
opc_mem = {}
for t in opc_terms:
    for g in str(bx.loc[t, "Lead_genes"]).split(";"):
        if g in vst.columns:
            opc_mem[g] = opc_mem.get(g, 0) + 1
for g in sorted([g for g, n in opc_mem.items() if n >= 3 and l2fc.get(g, 0) < 0], key=padj_of)[:14]:
    add("OPC", g)
ol_terms = [t for t in bx.index if _re.search(r"MOL|NFOL", t)
            and pd.to_numeric(bx.loc[t, "NES"]) > 0 and pd.to_numeric(bx.loc[t, "FDR q-val"]) < 0.25]
ol_mem = {}
for t in ol_terms:
    for g in str(bx.loc[t, "Lead_genes"]).split(";"):
        if g in vst.columns:
            ol_mem[g] = ol_mem.get(g, 0) + 1
for g in sorted([g for g in ol_mem if l2fc.get(g, 0) > 0], key=padj_of)[:10]:
    add("OL_diff", g)

gdf = pd.DataFrame(groups, columns=["sig", "gene"])

# ---- narrative block order (display only; membership/dedup unchanged) ----
DISPLAY_ORDER = ["OPC", "OL_diff", "BCAA_degradation", "Propanoate_metab",
                 "Lipid_FA", "Glycolysis", "OXPHOS"]
gdf["sig"] = pd.Categorical(gdf["sig"], categories=DISPLAY_ORDER, ordered=True)
gdf = gdf.sort_values("sig", kind="stable").reset_index(drop=True)   # stable: keep within-block order
SIGS = DISPLAY_ORDER

Z = vst[gdf.gene].copy(); Z = (Z - Z.mean(0)) / Z.std(0, ddof=0)
Zm = Z.T.values
nrows, ncols = len(gdf), len(order)
print("rows:", nrows)

# ---------------------------------------------------------------- compact plot
# explicit inch-based geometry: set the heatmap BODY height directly so a 30% cut is exact.
# compact v1 renders a ~5.21 in tall body (81 rows); 30% shorter -> 3.65 in.
BODY_IN = 3.65
ROW_H = BODY_IN / nrows                    # effective in/row (all text sizes stay fixed)
TOP_IN, BOTTOM_IN = 0.55, 0.28             # margins (annotation strip above; x-area below)
fig_h = TOP_IN + BODY_IN + BOTTOM_IN
fig = plt.figure(figsize=(7.5, fig_h))
ax = fig.add_axes([0.14, BOTTOM_IN / fig_h, 0.24, BODY_IN / fig_h])
im = ax.imshow(Zm, cmap=HM["colormap"], vmin=HM["vmin"], vmax=HM["vmax"], aspect="auto")

# top condition annotation
cond_of = m2["condition"].to_dict()
bar_h = 0.15 / ROW_H; bar_top = -0.5 - bar_h   # annotation strip (drawn above axes, fixed inches)
for j, s in enumerate(order):
    ax.add_patch(Rectangle((j - 0.5, bar_top), 1, bar_h,
                           color=COND_COLORS[cond_of[s]], clip_on=False, lw=0))
n_ctrl = sum(cond_of[s] == "control" for s in order)
ax.text((n_ctrl - 1) / 2, bar_top + bar_h / 2, PAL["conditions"]["control"]["label"],
        ha="center", va="center", fontsize=9, fontweight="bold", color="white", clip_on=False)
ax.text(n_ctrl + (ncols - n_ctrl - 1) / 2, bar_top + bar_h / 2, PAL["conditions"]["BCAAr"]["label"],
        ha="center", va="center", fontsize=9, fontweight="bold", color="white", clip_on=False)
ax.axvline(n_ctrl - 0.5, color=GAP["color"], lw=GAP["linewidth"], zorder=5)
ax.set_xticks([]); ax.set_yticks([])

# category strip + abbreviated labels
NICE_SHORT = {"Glycolysis": "Glyc.", "OXPHOS": "OXPHOS", "Propanoate_metab": "Prop.",
              "BCAA_degradation": "BCAA", "Lipid_FA": "Lipid/FA", "OPC": "OPC", "OL_diff": "OL"}
strip_right = -0.6; strip_w = 0.10 * ncols; strip_left = strip_right - strip_w
name_x = strip_left - 0.05 * ncols
start = 0
for t in SIGS:
    n = int((gdf.sig == t).sum())
    if n == 0:
        continue
    ax.add_patch(Rectangle((strip_left, start - 0.5), strip_w, n,
                           color=sig_colors[t], clip_on=False, lw=0))
    lbl = NICE_SHORT[t]
    lbl_fs = max(5.0, min(8.5, n * ROW_H * 110 / len(lbl)))
    ax.text(name_x, start + n / 2 - 0.5, lbl, rotation=90, va="center", ha="center",
            fontsize=lbl_fs, fontweight="bold", color=sig_colors[t], clip_on=False)
    if start > 0:
        ax.axhline(start - 0.5, color=GAP["color"], lw=GAP["linewidth"], zorder=5)
    start += n
ax.set_xlim(name_x - 0.05 * ncols, ncols - 0.5)
ax.set_ylim(nrows - 0.5, -0.5)   # exactly the data rows -> body height == BODY_IN (annotation sits above)
for sp in ax.spines.values():
    sp.set_visible(False)

# ---- callout labels with leader lines (blended: x = axes frac, y = data rows) ----
blend = blended_transform_factory(ax.transAxes, ax.transData)
PROP = {"Acads", "Hibch", "Hadha", "Hsd17b10", "Bckdha"}   # upstream of propionyl-CoA
label_specs = [("Tpi1", False), ("Eno2", False), ("Sdha", False), ("Ndufv1", False),
               ("Acads", True), ("Hibch", True), ("Hadha", True), ("Hsd17b10", True), ("Bckdha", True),
               ("Elovl5", False), ("Elovl6", False),
               ("Olig1", False), ("Olig2", False), ("Cspg4", False),
               ("Tppp", False), ("Lingo2", False)]
pos_of = {g: i for i, g in enumerate(gdf.gene)}
label_specs = [(g, b) for g, b in label_specs if g in pos_of]
label_specs.sort(key=lambda x: pos_of[x[0]])
src = [pos_of[g] for g, _ in label_specs]

MIN_GAP = 7.0 / (0.9 * ROW_H * 72)      # min vertical spacing (rows) so labels don't overlap
ys = []
for p in src:
    y = p if not ys else max(p, ys[-1] + MIN_GAP)
    ys.append(y)
overflow = ys[-1] - (nrows - 1)
if overflow > 0:
    ys = [y - overflow for y in ys]

x_lead0, x_text = 1.005, 1.05
for (g, bold), p, y in zip(label_specs, src, ys):
    ax.plot([x_lead0, x_text - 0.004], [p, y], transform=blend,
            color="0.45", lw=0.6, clip_on=False, zorder=6)
    ax.text(x_text, y, g, transform=blend, ha="left", va="center", fontsize=7,
            fontweight="bold" if bold else "normal", clip_on=False, zorder=7)

# propionyl-CoA bracket + group designation
py = [y for (g, _), y in zip(label_specs, ys) if g in PROP]
if py:
    y0, y1, x_br = min(py), max(py), 1.34
    for seg in ([[x_br, x_br], [y0, y1]], [[x_br - 0.02, x_br], [y0, y0]], [[x_br - 0.02, x_br], [y1, y1]]):
        ax.plot(seg[0], seg[1], transform=blend, color="0.2", lw=1.1, clip_on=False, zorder=6)
    ax.text(x_br + 0.03, (y0 + y1) / 2, "upstream of\npropionyl-CoA\nhandling",
            transform=blend, ha="left", va="center", fontsize=7.5, clip_on=False, zorder=7)

# colourbar
cax = fig.add_axes([0.66, 0.32, 0.02, 0.36])
cb = fig.colorbar(im, cax=cax); cb.set_label("Row Z-score (VST)", fontsize=9)

fig.savefig(f"{FIG}/2HIT_LE_heatmap_compact_short.png", dpi=200, bbox_inches="tight")
fig.savefig(f"{FIG}/2HIT_LE_heatmap_compact_short.svg", bbox_inches="tight")
fig.savefig(f"{FIG}/2HIT_LE_heatmap_compact_short.pdf", bbox_inches="tight")
plt.close(fig)
print("saved", f"{FIG}/2HIT_LE_heatmap_compact_short.png / .svg")
