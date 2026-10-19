#!/usr/bin/env python3
"""
Merged IleuR-vs-DMSO dot plot pulling signatures from three GSEA collections
(all computed on the same Wald-stat rankings, analysis_excl_outliers):

  Metabolic / Reactome : SREBF/SREBP gene activation           (GSEA_SREBP2)
  Braindex OPC-P (down): P-p3_OPC-P, P-p6_OPAS-P               (GSEA_braindex)
  Braindex astrocyte   : top-2 by mean positive NES            (GSEA_braindex)
  Hallmark             : MYC Targets V1                         (GSEA_IleuR_metabolic)

Columns: BT245, DIPGXIII. Categories separated by thin divider lines.
Title: "IleuR vs Ctrl".

Usage:
    python3 ileur_merged_dotplot.py --root analysis_excl_outliers
"""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))  # publication: repo-relative root
import argparse, json, re, warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["svg.fonttype"] = "none"      # editable <text> in SVG
matplotlib.rcParams["pdf.fonttype"] = 42           # embed TrueType (editable in AI)
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D

warnings.filterwarnings("ignore")

BASE = Path(REPO_ROOT + "/plasmidsaurusRNAseq")
PAL = json.load(open(REPO_ROOT + "/"
                     "BCAT1_figure_palette.json"))
ALPHA = 0.05
CELL_LINES = ["BT245", "DIPGXIII"]

CAT_COLORS = {
    "Metabolic": "#1a5276",
    "OPC-P": PAL["pathways"]["OPC"]["hex"],
    "Astrocyte": "#2e86c1",
    "Hallmark": "#7b3294",
}


def load_all(root):
    # SREBP (Reactome) — comparison col like "BT245 | IleuR"
    s = pd.read_csv(root / "GSEA_SREBP2" / "GSEA_results_all.csv")
    s = s[s.comparison.str.contains("IleuR") &
          s.gene_set.str.contains("gene activation")]
    s["cell_line"] = s.comparison.str.split(" | ", regex=False).str[0]
    srebp = s.set_index("cell_line")[["NES", "fdr"]]

    # braindex (IleuR)
    b = pd.read_csv(root / "GSEA_braindex" / "GSEA_braindex_all_results.csv")
    b = b[b.contrast == "IleuR"]

    def bx(sig):
        d = b[b.signature == sig].set_index("cell_line")
        return d[["NES", "fdr"]]

    # top-2 astrocytes by mean positive NES
    ast = b[b.signature.str.contains(r"_(ASTR|ASTP)", regex=True)]
    piv = ast.pivot_table(index="signature", columns="cell_line", values="NES")
    piv["mean"] = piv[CELL_LINES].mean(axis=1)
    top_ast = list(piv.sort_values("mean", ascending=False).head(2).index)

    # Hallmark MYC V1 (IleuR)
    h = pd.read_csv(root / "GSEA_IleuR_metabolic" / "GSEA_IleuR_all_results.csv")
    myc = h[h.gene_set == "Hallmark: Myc Targets V1"].set_index("cell_line")[["NES", "fdr"]]

    # assemble row list: (category, display label, per-cell-line NES/fdr frame)
    rows = [("Metabolic", "SREBP gene activation", srebp)]
    for sig in ["P-p3_OPC-P", "P-p6_OPAS-P"]:
        rows.append(("OPC-P", sig, bx(sig)))
    for sig in top_ast:
        rows.append(("Astrocyte", sig, bx(sig)))
    rows.append(("Hallmark", "MYC Targets V1", myc))
    return rows


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--root", required=True)
    a = p.parse_args()
    root = BASE / a.root
    out = root / "IleuR_merged_dotplot"
    out.mkdir(parents=True, exist_ok=True)

    rows = load_all(root)
    labels = [lab for _, lab, _ in rows]
    cats = [c for c, _, _ in rows]
    nes = np.array([[fr.loc[cl, "NES"] for cl in CELL_LINES] for _, _, fr in rows])
    fdr = np.array([[fr.loc[cl, "fdr"] for cl in CELL_LINES] for _, _, fr in rows])
    pd.DataFrame(nes, index=labels, columns=[f"{c}_NES" for c in CELL_LINES]).join(
        pd.DataFrame(fdr, index=labels, columns=[f"{c}_FDR" for c in CELL_LINES])
    ).to_csv(out / "IleuR_merged_values.csv")

    vmax = np.nanmax(np.abs(nes))
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    size = lambda q: 60 + 300 * np.clip(-np.log10(np.clip(q, 1e-4, 1)) / 4, 0, 1)

    xpos = [0.0, 0.62]
    fig, ax = plt.subplots(figsize=(5.0, 0.44 * len(rows) + 1.7))
    for i in range(len(rows)):
        for j in range(len(CELL_LINES)):
            v, q = nes[i, j], fdr[i, j]
            ax.scatter(xpos[j], i, s=size(q), c=[plt.cm.RdBu_r(norm(v))],
                       edgecolor="#111111" if q < ALPHA else "#c4c4c4",
                       linewidth=1.3 if q < ALPHA else 0.6, zorder=3)
    # column headers separated
    ax.set_xticks([])
    y_hdr = -0.78
    ax.text(xpos[0] + 0.18, y_hdr, CELL_LINES[0], ha="right", va="bottom",
            fontsize=10.5, clip_on=False)
    ax.text(xpos[1] - 0.18, y_hdr, CELL_LINES[1], ha="left", va="bottom",
            fontsize=10.5, clip_on=False)
    ax.text(np.mean(xpos), -1.5, "IleuR vs Ctrl", ha="center", va="bottom",
            fontsize=13, fontweight="bold", clip_on=False)

    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(labels, fontsize=9.5)
    for lbl, c in zip(ax.get_yticklabels(), cats):
        lbl.set_color(CAT_COLORS[c])
    ax.invert_yaxis()
    ax.set_xlim(-0.45, xpos[1] + 0.28); ax.set_ylim(len(rows) - 0.5, -0.6)

    # thin divider lines between categories
    for i in range(1, len(rows)):
        if cats[i] != cats[i - 1]:
            ax.axhline(i - 0.5, color="#999999", lw=0.9, zorder=1)
    ax.grid(axis="y", color="#f4f4f4", lw=0.6, zorder=0); ax.set_axisbelow(True)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)

    sm = plt.cm.ScalarMappable(cmap="RdBu_r", norm=norm); sm.set_array([])
    cb = fig.colorbar(sm, ax=ax, fraction=0.05, pad=0.04); cb.set_label("NES", fontsize=9)
    handles = [Line2D([], [], marker="o", ls="", markersize=np.sqrt(size(q)) * 0.5,
                      markerfacecolor="#d9d9d9", markeredgecolor="#666666",
                      label=("FDR ≤ 1e-4" if q <= 1e-4 else f"FDR {q:g}"))
               for q in [0.05, 0.001, 1e-4]]
    ax.legend(handles=handles, frameon=False, fontsize=8, loc="upper left",
              bbox_to_anchor=(1.30, 1.0), labelspacing=1.2)
    fig.tight_layout()
    for ext in ("png", "svg", "pdf"):
        fig.savefig(out / f"IleuR_merged_dotplot.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(pd.DataFrame(nes, index=labels, columns=CELL_LINES).round(2).to_string())
    print("\nwrote:", out / "IleuR_merged_dotplot.png")
