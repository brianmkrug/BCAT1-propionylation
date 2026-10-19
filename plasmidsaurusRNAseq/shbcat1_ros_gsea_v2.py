#!/usr/bin/env python3
"""
Compact ROS dot plot: ROS Pathway [Hallmark] + the two glutathione sets,
BT245 and DIPGXIII shBCAT1, with the FDR printed above each of the 6 dots.
Reads the already-computed GSEA results.

Usage:
    python3 shbcat1_ros_gsea_v2.py --root analysis_excl_outliers
"""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))  # publication: repo-relative root
import argparse
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

BASE = Path(REPO_ROOT + "/plasmidsaurusRNAseq")
ALPHA = 0.05
CELL_LINES = ["BT245", "DIPGXIII"]
ROWS = [
    ("ROS Pathway [Hallmark]", "ROS Pathway [Hallmark]", "#c0392b"),
    ("Glutathione metabolism [KEGG]", "Glutathione metabolism [KEGG]", "#1e8449"),
    ("Glutathione conjugation [Reactome]", "Glutathione conjugation [Reactome]", "#1e8449"),
]


def fdr_txt(q):
    return "<1e-4" if q < 1e-4 else f"{q:.2f}" if q >= 0.01 else f"{q:.1e}"


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--root", required=True)
    a = p.parse_args()
    out = BASE / a.root / "GSEA_shBCAT1_ROS"
    G = pd.read_csv(out / "GSEA_shBCAT1_ROS_results.csv")
    nes = G.pivot(index="gene_set", columns="cell_line", values="NES")
    fdr = G.pivot(index="gene_set", columns="cell_line", values="fdr")
    keys = [k for k, _, _ in ROWS]
    labels = [l for _, l, _ in ROWS]
    lab_c = [c for _, _, c in ROWS]
    NES = nes.loc[keys, CELL_LINES]
    FDR = fdr.loc[keys, CELL_LINES]

    norm = TwoSlopeNorm(vmin=-2, vcenter=0, vmax=2)      # expanded fixed NES scale
    size = lambda q: 90 + 320 * np.clip(-np.log10(np.clip(q, 1e-4, 1)) / 4, 0, 1)

    xpos = [0.0, 1.05]
    fig, ax = plt.subplots(figsize=(6.0, 3.4))
    for i, k in enumerate(keys):
        for j, cl in enumerate(CELL_LINES):
            v, q = NES.loc[k, cl], FDR.loc[k, cl]
            ax.scatter(xpos[j], i, s=size(q), c=[plt.cm.RdBu_r(norm(v))],
                       edgecolor="#111111" if q < ALPHA else "#c4c4c4",
                       linewidth=1.3 if q < ALPHA else 0.6, zorder=3)
            ax.text(xpos[j], i - 0.36, f"FDR {fdr_txt(q)}", ha="center", va="bottom",
                    fontsize=8, color="#333333", zorder=4)
    ax.set_xticks([])
    y_hdr = -0.74
    for xp, cl in zip(xpos, CELL_LINES):
        ax.text(xp, y_hdr, cl, ha="center", va="bottom", fontsize=10.5, clip_on=False)
    ax.text(np.mean(xpos), -1.3, "shBCAT1/shNT", ha="center", va="bottom",
            fontsize=12.5, fontweight="bold", clip_on=False)
    ax.set_yticks(range(len(keys))); ax.set_yticklabels(labels, fontsize=9.5)
    for lbl, c in zip(ax.get_yticklabels(), lab_c):
        lbl.set_color(c)
    ax.invert_yaxis()
    ax.set_xlim(-0.55, xpos[1] + 0.55); ax.set_ylim(len(keys) - 0.4, -0.9)
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
              bbox_to_anchor=(1.30, 1.0), labelspacing=1.3)
    fig.tight_layout()
    for ext in ("png", "svg", "pdf"):
        fig.savefig(out / f"GSEA_shBCAT1_ROS_glutathione_dotplot.{ext}", dpi=300,
                    bbox_inches="tight")
    plt.close(fig)
    print(pd.concat({"NES": NES.round(2), "FDR": FDR.round(3)}).to_string())
    print("\nwrote:", out / "GSEA_shBCAT1_ROS_glutathione_dotplot.png")
