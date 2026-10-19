#!/usr/bin/env python3
"""2HIT metabolic GSEA dot plot (BCAAr vs control).

Focuses on 2HIT only (4HIT excluded). Curated metabolic panel, showing only pathways
with FDR q < 0.25. Each pathway is one dot:
  * x position + fill colour = NES (RdBu_r gradient, symmetric about 0);
  * dot size = FDR significance, -log10(q)  (bigger = more significant).
Saves PNG + SVG (editable text for Illustrator).
"""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))  # publication: repo-relative root
import json
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["svg.fonttype"] = "none"   # SVG: real <text> elements, editable
matplotlib.rcParams["pdf.fonttype"] = 42       # PDF: embed CID TrueType (editable type)
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = REPO_ROOT + ""
GSEA = f"{ROOT}/BCAAr_analysis/gsea"
FIG = f"{ROOT}/BCAAr_analysis/figures"

with open(f"{ROOT}/BCAT1_figure_palette.json") as f:
    PAL = json.load(f)
CMAP = PAL["heatmap"]["colormap"]        # RdBu_r

Q_CUTOFF = 0.25
NICE = {"Propanoate_metab": "Propionate metab.", "OXPHOS": "OXPHOS",
        "BCAA_degradation": "BCAA degradation", "Glycolysis": "Glycolysis",
        "One-carbon/folate": "One-carbon / folate", "FA_oxidation": "FA oxidation",
        "FA_elongation": "FA elongation", "Unsat_FA_synthesis": "Unsat. FA synthesis",
        "Butanoate_metab": "Butanoate metab.", "Gly/Ser/Thr_metab": "Gly / Ser / Thr metab.",
        "Amino/nucleotide_sugar": "Amino / nucleotide sugar",
        "HALLMARK_Adipogenesis": "Adipogenesis (Hallmark)",
        "HALLMARK_Fatty_acid_metab": "Fatty-acid metab. (Hallmark)"}

d = pd.read_csv(f"{GSEA}/Curated_metabolic_panel_2HIT.csv").set_index("Term")
d["NES"] = pd.to_numeric(d["NES"]); d["FDR q-val"] = pd.to_numeric(d["FDR q-val"])
d = d[d["FDR q-val"] < Q_CUTOFF].copy()
d = d.sort_values("NES")                 # ascending -> most enriched ends up at top of y-axis
d["label"] = [NICE.get(t, t) for t in d.index]
d["neglogq"] = -np.log10(d["FDR q-val"].clip(lower=1e-4))

y = np.arange(len(d))
vmax = np.nanmax(np.abs(d["NES"]))
size = 30 + 260 * d["neglogq"]           # marker area from FDR significance

fig, ax = plt.subplots(figsize=(6.2, max(4.5, 0.42 * len(d))))
ax.axvline(0, color="0.6", lw=0.8, zorder=0)
sc = ax.scatter(d["NES"], y, s=size, c=d["NES"], cmap=CMAP, vmin=-vmax, vmax=vmax,
                edgecolor="0.25", linewidth=0.6, zorder=3)

ax.set_yticks(y); ax.set_yticklabels(d["label"], fontsize=9)
ax.set_xlabel("NES (BCAAr vs control)", fontsize=10)
ax.set_ylim(-0.6, len(d) - 0.4)
ax.margins(x=0.18)
ax.tick_params(axis="y", length=0)
for s in ("top", "right", "left"):
    ax.spines[s].set_visible(False)

# NES colour bar
cb = fig.colorbar(sc, ax=ax, fraction=0.045, pad=0.02)
cb.set_label("NES", fontsize=10)

# FDR size legend
q_ref = [0.2, 0.1, 0.05]
handles = [Line2D([0], [0], marker="o", linestyle="none", markerfacecolor="0.7",
                  markeredgecolor="0.25",
                  markersize=np.sqrt(30 + 260 * (-np.log10(q))) , label=f"{q:g}")
           for q in q_ref]
leg = ax.legend(handles=handles, title="FDR q-val", loc="lower right",
                frameon=False, labelspacing=1.3, borderpad=1.0, handletextpad=1.0,
                fontsize=8.5, title_fontsize=9)

plt.tight_layout()
plt.savefig(f"{FIG}/GSEA_metabolic_2HIT_dotplot.png", dpi=200, bbox_inches="tight")
plt.savefig(f"{FIG}/GSEA_metabolic_2HIT_dotplot.svg", bbox_inches="tight")
plt.savefig(f"{FIG}/GSEA_metabolic_2HIT_dotplot.pdf", bbox_inches="tight")
plt.close()
print(f"{len(d)} pathways (q<{Q_CUTOFF}) plotted")
print("saved", f"{FIG}/GSEA_metabolic_2HIT_dotplot.png / .svg")
