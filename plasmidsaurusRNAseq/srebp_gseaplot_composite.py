#!/usr/bin/env python3
"""
Composite of the two SREBP gene-activation enrichment plots stacked vertically
(BT245 top, DIPG-XIII bottom), with minimal panel heights to economise space.
Shared x-axis; running-ES + gene-rank ticks; ~5 leading-edge genes labelled.

Usage:
    python3 srebp_gseaplot_composite.py --root analysis_excl_outliers --nlabel 5
"""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))  # publication: repo-relative root
import argparse, json, warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["pdf.fonttype"] = 42
import matplotlib.pyplot as plt
from matplotlib import gridspec
import gseapy

warnings.filterwarnings("ignore")

BASE = Path(REPO_ROOT + "/plasmidsaurusRNAseq")
OLD = Path(REPO_ROOT + "/shBCAT1_inducible_analysis")
PAL = json.load(open(REPO_ROOT + "/"
                     "BCAT1_figure_palette.json"))
SET = "Activation of Gene Expression by SREBF (SREBP)"
LABEL = "SREBP gene activation (Reactome)"
LINE = "#69b023"
SREBP_C = PAL["metabolic_categories"]["Steroid biosynthesis"]["hex"]


def rank_bt245(root):
    d = pd.read_csv(root / "DE_output" / "DE_BT245_shBCAT1_vs_shNT.csv", index_col=0)
    d = d.dropna(subset=["stat", "gene_name"]); d = d[d.gene_name.astype(str) != ""]
    return d.groupby(d.gene_name.str.upper())["stat"].mean().sort_values(ascending=False)


def rank_dipg_old():
    d = pd.read_csv(OLD / "out" / "deseq2_shBCAT1_vs_shNT.csv", index_col=0)
    d.index = d.index.astype(str).str.upper(); d = d.dropna(subset=["stat"])
    return d.groupby(level=0)["stat"].mean().sort_values(ascending=False)


def panel(fig, gs_es, gs_tick, rnk, genes, cl, nlabel, share, xlabel):
    pre = gseapy.prerank(rnk=rnk.reset_index(), gene_sets={SET: genes},
                         permutation_num=1000, min_size=8, max_size=600,
                         seed=0, threads=4, outdir=None, verbose=False)
    res = pre.results[SET]
    RES = np.asarray(res["RES"]); hits = np.asarray(res["hits"])
    row = pre.res2d.iloc[0]; nes = float(row["NES"]); fdr = float(row["FDR q-val"])
    stars = ("****" if fdr < 1e-4 else "***" if fdr < 1e-3 else
             "**" if fdr < 1e-2 else "*" if fdr < 0.05 else "ns")
    go = list(rnk.index); n = len(RES); peak = int(np.argmax(np.abs(RES))); neg = nes < 0
    le_hits = [h for h in hits if h >= peak] if neg else [h for h in hits if h <= peak]
    le_sorted = sorted([(h, go[h]) for h in le_hits],
                       key=lambda x: -x[0] if neg else x[0])[:nlabel]

    ax = fig.add_subplot(gs_es, sharex=share)
    ax.plot(np.arange(n), RES, color=LINE, lw=1.6, zorder=2)
    ax.axhline(0, color="#888888", lw=0.7)
    ax.axvline(peak, color="#333333", ls="--", lw=0.7)
    ax.plot(peak, RES[peak], "v" if neg else "^", color="#c0392b", ms=7, zorder=5)
    ax.set_ylabel("Enrichment\nScore", fontsize=8.5)
    ax.set_xlim(0, n); ax.set_xticks([]); ax.set_facecolor("white")
    ax.set_title(f"{cl}   {LABEL}  (NES={nes:.2f}, p.adjust={fdr:.2g}, {stars})",
                 fontsize=10, loc="left", pad=3)
    ax.tick_params(labelsize=8)
    lab = sorted(le_sorted, key=lambda x: RES[x[0]], reverse=True)
    for (h, g), yf in zip(lab, np.linspace(0.86, 0.14, len(lab))):
        ax.plot(h, RES[h], "o", color="#111111", ms=3.2, zorder=6, clip_on=False)
        ax.annotate(g, xy=(h, RES[h]), xycoords="data", xytext=(1.015, yf),
                    textcoords=ax.transAxes, fontsize=8, fontweight="bold",
                    color=SREBP_C, ha="left", va="center", zorder=7,
                    annotation_clip=False,
                    arrowprops=dict(arrowstyle="-", color="#999999", lw=0.6,
                                    shrinkA=1, shrinkB=2))
    ax.spines[["top", "right"]].set_visible(False)

    axt = fig.add_subplot(gs_tick, sharex=ax)
    axt.vlines(hits, 0, 1, color="#222222", lw=0.6)
    grad = np.linspace(0, 1, 256).reshape(1, -1)
    axt.imshow(grad, aspect="auto", cmap="RdBu_r", extent=[0, n, -0.55, -0.05], alpha=0.9)
    axt.set_ylim(-0.55, 1); axt.set_yticks([]); axt.set_xlim(0, n)
    if xlabel:
        axt.set_xlabel("Rank in Ordered Dataset", fontsize=9)
        axt.tick_params(labelsize=8)
    else:
        axt.set_xticks([])
    for s in axt.spines.values():
        s.set_visible(False)
    return dict(cell_line=cl, NES=nes, FDR=fdr, labelled=";".join(g for _, g in le_sorted))


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--root", required=True)
    p.add_argument("--nlabel", type=int, default=5)
    a = p.parse_args()
    root = BASE / a.root
    out = root / "COMPARE_inducible_vs_BT245" / "srebp_gseaplots"
    out.mkdir(parents=True, exist_ok=True)
    genes = sorted({g.split(",")[0].strip().upper()
                    for g in gseapy.get_library("Reactome_Pathways_2024", "Human")[SET]})

    # minimal heights: two short ES panels + thin tick strips
    fig = plt.figure(figsize=(7.4, 4.2))
    gs = gridspec.GridSpec(5, 1, height_ratios=[3, 0.35, 0.45, 3, 0.35], hspace=0.08)
    s1 = panel(fig, gs[0], gs[1], rank_bt245(root), genes, "BT245", a.nlabel, None, False)
    s2 = panel(fig, gs[3], gs[4], rank_dipg_old(), genes, "DIPGXIII",
               a.nlabel, None, True)

    for ext in ("png", "svg", "pdf"):
        fig.savefig(out / f"gseaplot_SREBP_composite.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(s1); print(s2)
    print("\nwrote:", out / "gseaplot_SREBP_composite.png")
