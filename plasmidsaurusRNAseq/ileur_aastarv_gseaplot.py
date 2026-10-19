#!/usr/bin/env python3
"""
Simplified GSEA enrichment plots (running ES + gene-rank ticks, no ranked-metric
panel) for GO 'Response to Amino Acid Starvation', IleuR vs DMSO, with BT245 and
DIPGXIII stacked vertically.

Usage:
    python3 ileur_aastarv_gseaplot.py --root analysis_excl_outliers
"""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))  # publication: repo-relative root
import argparse, json, warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["svg.fonttype"] = "none"      # editable <text> in SVG
matplotlib.rcParams["pdf.fonttype"] = 42           # embed TrueType (editable in AI)
import matplotlib.pyplot as plt
from matplotlib import gridspec
import gseapy

warnings.filterwarnings("ignore")

BASE = Path(REPO_ROOT + "/plasmidsaurusRNAseq")
SET = "Response to AA starvation [GO]"
LINE = "#e8862a"                                   # ISR/AA-deprivation accent (orange)


def ranking(root, cl):
    d = pd.read_csv(root / "DE_output" / f"DE_{cl}_IleuR_vs_DMSO.csv", index_col=0)
    d = d.dropna(subset=["stat", "gene_name"])
    d = d[d.gene_name.astype(str) != ""]
    return d.groupby(d.gene_name.str.upper())["stat"].mean().sort_values(ascending=False)


def panel(fig, gs_es, gs_tick, rnk, gmt, cl, share=None):
    pre = gseapy.prerank(rnk=rnk.reset_index(), gene_sets={SET: gmt[SET]},
                         permutation_num=1000, min_size=8, max_size=600,
                         seed=0, threads=4, outdir=None, verbose=False)
    res = pre.results[SET]
    RES = np.asarray(res["RES"]); hits = np.asarray(res["hits"])
    row = pre.res2d.iloc[0]
    nes = float(row["NES"]); fdr = float(row["FDR q-val"])
    stars = ("****" if fdr < 1e-4 else "***" if fdr < 1e-3 else
             "**" if fdr < 1e-2 else "*" if fdr < 0.05 else "ns")
    n = len(RES); peak = int(np.argmax(np.abs(RES)))
    neg = nes < 0

    ax = fig.add_subplot(gs_es, sharex=share)
    ax.plot(np.arange(n), RES, color=LINE, lw=1.8)
    ax.axhline(0, color="#888888", lw=0.8)
    ax.axvline(peak, color="#333333", ls="--", lw=0.8)
    ax.plot(peak, RES[peak], "v" if neg else "^", color="#c0392b", ms=8, zorder=5)
    ax.set_ylabel("Enrichment\nScore", fontsize=9)
    ax.set_xlim(0, n); ax.set_xticks([])
    ax.set_facecolor("#f6efe6")
    ax.set_title(f"{cl}   {SET}  (NES={nes:.2f}, p.adjust={fdr:.2g}, {stars})",
                 fontsize=11, loc="left")
    ax.spines[["top", "right"]].set_visible(False)

    axt = fig.add_subplot(gs_tick, sharex=ax)
    axt.vlines(hits, 0, 1, color="#222222", lw=0.6)
    grad = np.linspace(0, 1, 256).reshape(1, -1)
    axt.imshow(grad, aspect="auto", cmap="RdBu_r", extent=[0, n, -0.55, -0.05], alpha=0.9)
    axt.set_ylim(-0.55, 1); axt.set_yticks([]); axt.set_xlim(0, n)
    for s in axt.spines.values():
        s.set_visible(False)
    return ax, axt, dict(cell_line=cl, NES=nes, FDR=fdr, peak=peak, n_hits=len(hits))


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--root", required=True)
    a = p.parse_args()
    root = BASE / a.root
    out = root / "GSEA_IleuR_ISR_ATF4"
    out.mkdir(parents=True, exist_ok=True)
    gmt = json.load(open(out / "genesets_used.json"))

    fig = plt.figure(figsize=(7.6, 6.6))
    gs = gridspec.GridSpec(5, 1, height_ratios=[4, 0.5, 0.5, 4, 0.5], hspace=0.08)

    stats = []
    ax1, axt1, s1 = panel(fig, gs[0], gs[1], ranking(root, "BT245"), gmt, "BT245")
    stats.append(s1)
    ax2, axt2, s2 = panel(fig, gs[3], gs[4], ranking(root, "DIPGXIII"), gmt,
                          "DIPGXIII", share=ax1)
    stats.append(s2)
    axt2.set_xlabel("Rank in Ordered Dataset  (IleuR-up → IleuR-down)", fontsize=10)

    fig.tight_layout()
    for ext in ("png", "svg", "pdf"):
        fig.savefig(out / f"gseaplot_AAstarvation_stacked.{ext}", dpi=300,
                    bbox_inches="tight")
    plt.close(fig)
    print(pd.DataFrame(stats).to_string(index=False))
    print("\nwrote:", out / "gseaplot_AAstarvation_stacked.png")
