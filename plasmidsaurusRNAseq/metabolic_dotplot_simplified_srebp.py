#!/usr/bin/env python3
"""
Simplified metabolic dot plot (inducible DIPG-XIII vs new BT245 shBCAT1) with the
Reactome 'SREBF/SREBP gene activation' signature added as the top row.

KEGG rows are read from the precomputed metabolic GSEA; the SREBP row is computed
here on the same Wald-stat rankings (BT245 = new Plasmidsaurus shBCAT1,
DIPGXIII = OLD inducible shBCAT1) so it is directly comparable.

Usage:
    python3 metabolic_dotplot_simplified_srebp.py --root analysis_excl_outliers
"""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))  # publication: repo-relative root
import argparse, warnings
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
import gseapy

warnings.filterwarnings("ignore")

BASE = Path(REPO_ROOT + "/plasmidsaurusRNAseq")
OLD = Path(REPO_ROOT + "/shBCAT1_inducible_analysis")
ALPHA = 0.05
KEGG_C, REAC_C = "#1a5276", "#7b3294"
SREBP_TERM = "Activation of Gene Expression by SREBF (SREBP)"

# (source key, display label, colour)  — SREBP added at top
ROWS = [
    ("SREBP",                                  "SREBP gene activation",           REAC_C),
    ("KEGG: Steroid biosynthesis",             "Steroid biosynthesis",            KEGG_C),
    ("KEGG: Pantothenate and CoA biosynthesis","Pantothenate and CoA biosynthesis", KEGG_C),
    ("KEGG: Arginine and proline metabolism",  "Arginine and proline metabolism", KEGG_C),
    ("KEGG: Terpenoid backbone biosynthesis",  "Terpenoid backbone biosynthesis", KEGG_C),
    ("KEGG: Fatty acid degradation",           "Fatty acid metabolism",           KEGG_C),
    ("KEGG: Glycosaminoglycan degradation",    "Glycosaminoglycan degradation",   KEGG_C),
    ("KEGG: Other glycan degradation",         "Other glycan degradation",        KEGG_C),
]
DIVIDER_AFTER = 5   # rows 0..5 are the down block; up block below


def rank_bt245(root):
    d = pd.read_csv(root / "DE_output" / "DE_BT245_shBCAT1_vs_shNT.csv", index_col=0)
    d = d.dropna(subset=["stat", "gene_name"])
    d = d[d.gene_name.astype(str) != ""]
    return d.groupby(d.gene_name.str.upper())["stat"].mean().sort_values(ascending=False)


def rank_dipg_old():
    d = pd.read_csv(OLD / "out" / "deseq2_shBCAT1_vs_shNT.csv", index_col=0)
    d.index = d.index.astype(str).str.upper()
    d = d.dropna(subset=["stat"])
    return d.groupby(level=0)["stat"].mean().sort_values(ascending=False)


def srebp_nes(rnk, genes):
    pre = gseapy.prerank(rnk=rnk.reset_index(), gene_sets={"SREBP": genes},
                         permutation_num=1000, min_size=8, max_size=600,
                         seed=0, threads=4, outdir=None, verbose=False)
    r = pre.res2d.iloc[0]
    return float(r["NES"]), float(r["FDR q-val"])


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--root", required=True)
    a = p.parse_args()
    root = BASE / a.root
    out = root / "COMPARE_inducible_vs_BT245"
    cols = ["DIPGXIII", "BT245"]

    G = pd.read_csv(out / "GSEA_metabolic_DIPGold_vs_BT245_all.csv")
    nes_k = G.pivot(index="gene_set", columns="col", values="NES")
    fdr_k = G.pivot(index="gene_set", columns="col", values="fdr")

    srebp = {g.split(",")[0].strip().upper()
             for g in gseapy.get_library("Reactome_Pathways_2024", "Human")[SREBP_TERM]}
    srebp = sorted(srebp)
    n_bt, f_bt = srebp_nes(rank_bt245(root), srebp)
    n_dp, f_dp = srebp_nes(rank_dipg_old(), srebp)
    srebp_row = {"DIPGXIII": (n_dp, f_dp), "BT245": (n_bt, f_bt)}
    print(f"SREBP gene activation: DIPGXIII NES={n_dp:.2f} (FDR {f_dp:.3g}); "
          f"BT245 NES={n_bt:.2f} (FDR {f_bt:.3g})")

    labels = [l for _, l, _ in ROWS]
    lab_c = [c for _, _, c in ROWS]
    NES = np.zeros((len(ROWS), 2)); FDR = np.zeros((len(ROWS), 2))
    for i, (key, _, _) in enumerate(ROWS):
        for j, c in enumerate(cols):
            if key == "SREBP":
                NES[i, j], FDR[i, j] = srebp_row[c]
            else:
                NES[i, j], FDR[i, j] = nes_k.loc[key, c], fdr_k.loc[key, c]
    pd.DataFrame(NES, index=labels, columns=[f"{c}_NES" for c in cols]).join(
        pd.DataFrame(FDR, index=labels, columns=[f"{c}_FDR" for c in cols])
    ).to_csv(out / "GSEA_metabolic_simplified_srebp_values.csv")

    vmax = np.nanmax(np.abs(NES))
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    size = lambda q: 55 + 210 * np.clip(-np.log10(np.clip(q, 1e-4, 1)) / 4, 0, 1)

    xpos = [0.0, 0.8]
    fig, ax = plt.subplots(figsize=(6.2, 4.1))
    for i in range(len(ROWS)):
        for j in range(2):
            v, q = NES[i, j], FDR[i, j]
            ax.scatter(xpos[j], i, s=size(q), c=[plt.cm.RdBu_r(norm(v))],
                       edgecolor="#111111" if q < ALPHA else "#c4c4c4",
                       linewidth=1.3 if q < ALPHA else 0.6, zorder=3)
    ax.set_xticks(xpos); ax.set_xticklabels(cols, fontsize=11)
    ax.xaxis.tick_top()
    ax.text(np.mean(xpos), -1.45, "shBCAT1/shNT", ha="center", va="bottom",
            fontsize=12, fontweight="bold", clip_on=False)
    ax.set_yticks(range(len(ROWS))); ax.set_yticklabels(labels, fontsize=9.5)
    for lbl, c in zip(ax.get_yticklabels(), lab_c):
        lbl.set_color(c)
    ax.invert_yaxis()
    ax.set_xlim(-0.7, 1.5); ax.set_ylim(len(ROWS) - 0.5, -0.65)
    ax.axhline(DIVIDER_AFTER + 0.5, color="#cccccc", lw=0.9, zorder=0)
    ax.grid(axis="y", color="#f2f2f2", lw=0.6, zorder=0); ax.set_axisbelow(True)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)

    sm = plt.cm.ScalarMappable(cmap="RdBu_r", norm=norm); sm.set_array([])
    cb = fig.colorbar(sm, ax=ax, fraction=0.04, pad=0.10); cb.set_label("NES", fontsize=9)
    handles = [Line2D([], [], marker="o", ls="", markersize=np.sqrt(size(q)) * 0.5,
                      markerfacecolor="#d9d9d9", markeredgecolor="#666666",
                      label=("FDR ≤ 1e-4" if q <= 1e-4 else f"FDR {q:g}"))
               for q in [0.05, 0.001, 1e-4]]
    ax.legend(handles=handles, frameon=False, fontsize=8, loc="upper left",
              bbox_to_anchor=(1.30, 1.0), labelspacing=1.2)
    fig.tight_layout()
    for ext in ("png", "pdf", "svg"):
        fig.savefig(out / f"GSEA_metabolic_dotplot_simplified_srebp.{ext}", dpi=300,
                    bbox_inches="tight")
    plt.close(fig)
    print(pd.DataFrame(NES, index=labels, columns=cols).round(2).to_string())
    print("\nwrote:", out / "GSEA_metabolic_dotplot_simplified_srebp.png")
