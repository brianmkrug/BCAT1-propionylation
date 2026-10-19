#!/usr/bin/env python3
"""
Second version of the concordant braindex glial dot plot
(earlier inducible DIPG-XIII vs new BT245 shBCAT1).

Changes vs the first version:
  - RGC / OPC / Astrocyte capped to their top-5 signatures by |mean NES|
    (OL kept as-is, only 2 concordant).
  - two columns pulled close together
  - family labels sit tight against the signature list
  - "shBCAT1/shNT" title above the cell-line columns
  - saved as PNG + SVG

Usage:
    python3 compare_inducible_vs_BT245_concordant_top5.py --root analysis_excl_outliers
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
FAMILY_PATTERNS = [("RGC", r"RGC"), ("OPC", r"_(OPC|OPAS)"),
                   ("OL", r"_(NFOL|MOL)"), ("Astrocyte", r"_(ASTR|ASTP)")]
FAMILY_COLOR = {"RGC": "#8a6d3b", "OPC": PAL["pathways"]["OPC"]["hex"],
                "OL": PAL["pathways"]["OL_diff"]["hex"], "Astrocyte": "#2e86c1"}
FAM_ORDER = ["RGC", "OPC", "OL", "Astrocyte"]
TOPN = {"RGC": 5, "OPC": 5, "Astrocyte": 5, "OL": 99}   # OL: keep all concordant


def family(n):
    for f, pat in FAMILY_PATTERNS:
        if re.search(pat, n):
            return f
    return None


def is_prolif(n):
    return bool(re.search(r"-P(\d*)$", n))


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--root", required=True)
    a = p.parse_args()
    root = BASE / a.root
    out = root / "COMPARE_inducible_vs_BT245"

    M = pd.read_csv(out / "braindex_NES_old_vs_BT245.csv")
    M["family"] = M.signature.map(family)
    M["proliferating"] = M.signature.map(is_prolif)
    g = M[M.family.notna()].copy()
    g = g[np.sign(g.NES_old) == np.sign(g.NES_new)]      # sign-concordant only
    g["mean_nes"] = g[["NES_old", "NES_new"]].mean(axis=1)

    # per family: top-N by |mean NES|
    idx = []
    for f in FAM_ORDER:
        s = g[g.family == f].copy()
        s = s.reindex(s.mean_nes.abs().sort_values(ascending=False).index).head(TOPN[f])
        s = s.sort_values(["proliferating", "mean_nes"])   # display order
        idx += list(s.signature)
    g = g.set_index("signature").loc[idx]
    g.reset_index().to_csv(out / "braindex_glial_concordant_top5.csv", index=False)
    print(f"kept {len(idx)} signatures:",
          {f: sum(family(s) == f for s in idx) for f in FAM_ORDER})

    cols = ["DIPGXIII", "BT245"]
    nes = g[["NES_old", "NES_new"]].copy(); nes.columns = cols
    fdr = g[["fdr_old", "fdr_new"]].copy(); fdr.columns = cols

    vmax = np.nanmax(np.abs(nes.values))
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    size = lambda q: 42 + 250 * np.clip(-np.log10(np.clip(q, 1e-4, 1)) / 4, 0, 1)

    xpos = [0.0, 0.82]                                    # columns close but non-touching
    fig, ax = plt.subplots(figsize=(4.8, 0.40 * len(idx) + 1.5))
    for i, sg in enumerate(idx):
        for j, c in enumerate(cols):
            v, q = nes.loc[sg, c], fdr.loc[sg, c]
            ax.scatter(xpos[j], i, s=size(q), c=[plt.cm.RdBu_r(norm(v))],
                       edgecolor="#111111" if q < ALPHA else "#c4c4c4",
                       linewidth=1.3 if q < ALPHA else 0.6, zorder=3)
    ax.set_xticks([])
    # column headers, anchored apart so they never collide
    y_hdr = -0.82
    ax.text(xpos[0] + 0.22, y_hdr, cols[0], ha="right", va="bottom",
            fontsize=10.5, clip_on=False)
    ax.text(xpos[1] - 0.22, y_hdr, cols[1], ha="left", va="bottom",
            fontsize=10.5, clip_on=False)
    ax.text(np.mean(xpos), -1.55, "shBCAT1/shNT", ha="center", va="bottom",
            fontsize=12.5, fontweight="bold", clip_on=False)
    ax.set_yticks(range(len(idx)))
    ax.set_yticklabels([f"{s}  ●" if is_prolif(s) else s for s in idx], fontsize=8.5)
    for lbl, sg in zip(ax.get_yticklabels(), idx):
        lbl.set_color(FAMILY_COLOR[family(sg)])
        if is_prolif(sg):
            lbl.set_fontweight("bold")
    ax.invert_yaxis()
    ax.set_xlim(-0.55, xpos[1] + 0.28); ax.set_ylim(len(idx) - 0.4, -0.6)
    ax.grid(color="#f2f2f2", lw=0.6, zorder=0); ax.set_axisbelow(True)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)

    sm = plt.cm.ScalarMappable(cmap="RdBu_r", norm=norm); sm.set_array([])
    cb = fig.colorbar(sm, ax=ax, fraction=0.05, pad=0.04); cb.set_label("NES", fontsize=9)
    handles = [Line2D([], [], marker="o", ls="", markersize=np.sqrt(size(q)) * 0.55,
                      markerfacecolor="#d9d9d9", markeredgecolor="#666666",
                      label=("FDR ≤ 1e-4" if q <= 1e-4 else f"FDR {q:g}"))
               for q in [0.05, 0.001, 1e-4]]
    handles += [Line2D([], [], marker="o", ls="", markersize=6,
                       markerfacecolor="#555", markeredgecolor="#111",
                       label="● proliferating (-P)")]
    ax.legend(handles=handles, frameon=False, fontsize=8, loc="upper left",
              bbox_to_anchor=(1.32, 1.0), labelspacing=1.1)
    fig.tight_layout()

    # family labels — measured AFTER tight_layout so the axis geometry is final,
    # placed well left of the (variable-length) signature names.
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    inv = ax.transData.inverted()
    left_px = min(t.get_window_extent(renderer).x0 for t in ax.get_yticklabels())
    left_data = inv.transform((left_px, 0))[0]
    x_bar = left_data - 0.05                              # colored bracket, left of names
    x_txt = x_bar - 0.10                                  # family label, left of bracket
    fams = [family(s) for s in idx]
    for f in FAM_ORDER:
        pos = [i for i, ff in enumerate(fams) if ff == f]
        if not pos:
            continue
        lo, hi = min(pos), max(pos)
        ax.plot([x_bar, x_bar], [lo - 0.34, hi + 0.34], color=FAMILY_COLOR[f],
                lw=3.5, clip_on=False, zorder=5, solid_capstyle="round")
        ax.text(x_txt, (lo + hi) / 2, f, rotation=90, va="center", ha="center",
                fontsize=10.5, fontweight="bold", color=FAMILY_COLOR[f], clip_on=False)
        if lo > 0:
            ax.axhline(lo - 0.5, color="#cccccc", lw=0.8, zorder=0)

    for ext in ("png", "svg", "pdf"):
        fig.savefig(out / f"GSEA_braindex_glial_dotplot_concordant_top5.{ext}",
                    dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(g[["NES_old", "NES_new"]].round(2).to_string())
    print("\nwrote:", out / "GSEA_braindex_glial_dotplot_concordant_top5.png")
