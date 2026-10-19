#!/usr/bin/env python3
"""
Does BCAT1 knockdown induce ROS / oxidative-stress signatures?
GSEA of curated ROS/NRF2/glutathione/oxidative-stress gene sets on the shBCAT1
Wald-stat rankings for:
   BT245     = new Plasmidsaurus shBCAT1 vs shNT (outlier-excluded)
   DIPGXIII  = OLD inducible shBCAT1 vs shNT

Both columns on one dot plot (NES colour, FDR size), title "shBCAT1/shNT".
Positive NES = up in knockdown → would indicate ROS-stress induction.

Usage:
    python3 shbcat1_ros_gsea.py --root analysis_excl_outliers
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
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D
import gseapy

warnings.filterwarnings("ignore")

BASE = Path(REPO_ROOT + "/plasmidsaurusRNAseq")
OLD = Path(REPO_ROOT + "/shBCAT1_inducible_analysis")
ALPHA = 0.05

CAT_COLORS = {
    "Global ROS": "#c0392b",
    "NRF2 / detox": "#e8862a",
    "Glutathione / redox": "#1e8449",
    "GO oxidative response": "#7b3294",
}


def load_sets():
    hall = gseapy.get_library("MSigDB_Hallmark_2020", "Human")
    reac = gseapy.get_library("Reactome_Pathways_2024", "Human")
    kegg = gseapy.get_library("KEGG_2021_Human", "Human")
    wp = gseapy.get_library("WikiPathway_2023_Human", "Human")
    go = gseapy.get_library("GO_Biological_Process_2023", "Human")
    up = lambda lst: sorted({g.split(",")[0].strip().upper() for g in lst})
    rows = [
        ("Global ROS", "ROS Pathway [Hallmark]", up(hall["Reactive Oxygen Species Pathway"])),
        ("NRF2 / detox", "KEAP1-NFE2L2 pathway [Reactome]", up(reac["KEAP1-NFE2L2 Pathway"])),
        ("NRF2 / detox", "Detoxification of ROS [Reactome]",
         up(reac["Detoxification of Reactive Oxygen Species"])),
        ("NRF2 / detox", "NRF2 pathway [WikiPathways]", up(wp["NRF2 Pathway WP2884"])),
        ("NRF2 / detox", "Oxidative stress response [WikiPathways]",
         up(wp["Oxidative Stress Response WP408"])),
        ("NRF2 / detox", "FOXO transcription of oxidative-stress genes [Reactome]",
         up(reac["FOXO-mediated Transcription of Oxidative Stress, Metabolic and Neuronal Genes"])),
        ("Glutathione / redox", "Glutathione metabolism [KEGG]", up(kegg["Glutathione metabolism"])),
        ("Glutathione / redox", "Glutathione conjugation [Reactome]",
         up(reac["Glutathione Conjugation"])),
        ("Glutathione / redox", "Cell redox homeostasis [Reactome]",
         up(reac["Cell Redox Homeostasis"])),
        ("GO oxidative response", "Cellular response to oxidative stress [GO]",
         up(go["Cellular Response To Oxidative Stress (GO:0034599)"])),
        ("GO oxidative response", "Cellular response to ROS [GO]",
         up(go["Cellular Response To Reactive Oxygen Species (GO:0034614)"])),
        ("GO oxidative response", "Cellular response to H2O2 [GO]",
         up(go["Cellular Response To Hydrogen Peroxide (GO:0070301)"])),
    ]
    return rows


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


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--root", required=True)
    p.add_argument("--seed", type=int, default=0)
    a = p.parse_args()
    root = BASE / a.root
    out = root / "GSEA_shBCAT1_ROS"
    out.mkdir(parents=True, exist_ok=True)

    setrows = load_sets()
    gmt = {lab: g for _, lab, g in setrows}
    cats = {lab: c for c, lab, _ in setrows}
    json.dump({lab: g for _, lab, g in setrows}, open(out / "genesets_used.json", "w"), indent=2)

    ranks = {"BT245": rank_bt245(root), "DIPGXIII": rank_dipg_old()}
    res = []
    for cl, rnk in ranks.items():
        pre = gseapy.prerank(rnk=rnk.reset_index(), gene_sets=gmt, permutation_num=1000,
                             min_size=8, max_size=600, seed=a.seed, threads=4,
                             outdir=None, verbose=False)
        r = pre.res2d.copy(); r["cell_line"] = cl
        res.append(r)
        print(f"  {cl}: {len(rnk)} ranked, {len(r)} sets tested")
    G = pd.concat(res, ignore_index=True).rename(
        columns={"Term": "gene_set", "FDR q-val": "fdr", "NOM p-val": "pval"})
    for c in ["NES", "fdr", "pval"]:
        G[c] = pd.to_numeric(G[c], errors="coerce")
    G["category"] = G.gene_set.map(cats)
    G.sort_values(["cell_line", "category"]).to_csv(out / "GSEA_shBCAT1_ROS_results.csv",
                                                    index=False)

    CELL_LINES = ["BT245", "DIPGXIII"]
    order = [lab for _, lab, _ in setrows if lab in set(G.gene_set)]
    cat_of = [cats[l] for l in order]
    nes = G.pivot(index="gene_set", columns="cell_line", values="NES").reindex(order)[CELL_LINES]
    fdr = G.pivot(index="gene_set", columns="cell_line", values="fdr").reindex(order)[CELL_LINES]

    norm = TwoSlopeNorm(vmin=-2, vcenter=0, vmax=2)      # expanded fixed NES scale
    size = lambda q: 55 + 300 * np.clip(-np.log10(np.clip(q, 1e-4, 1)) / 4, 0, 1)

    xpos = [0.0, 0.62]
    fig, ax = plt.subplots(figsize=(6.8, 0.44 * len(order) + 1.6))
    for i, gs in enumerate(order):
        for j, cl in enumerate(CELL_LINES):
            v, q = nes.loc[gs, cl], fdr.loc[gs, cl]
            if np.isnan(v):
                continue
            ax.scatter(xpos[j], i, s=size(q), c=[plt.cm.RdBu_r(norm(v))],
                       edgecolor="#111111" if q < ALPHA else "#c4c4c4",
                       linewidth=1.3 if q < ALPHA else 0.6, zorder=3)
    ax.set_xticks([])
    y_hdr = -0.85
    ax.text(xpos[0] + 0.18, y_hdr, CELL_LINES[0], ha="right", va="bottom",
            fontsize=10.5, clip_on=False)
    ax.text(xpos[1] - 0.18, y_hdr, CELL_LINES[1], ha="left", va="bottom",
            fontsize=10.5, clip_on=False)
    ax.text(np.mean(xpos), -1.7, "shBCAT1/shNT", ha="center", va="bottom",
            fontsize=13, fontweight="bold", clip_on=False)
    ax.set_yticks(range(len(order))); ax.set_yticklabels(order, fontsize=8.5)
    for lbl, c in zip(ax.get_yticklabels(), cat_of):
        lbl.set_color(CAT_COLORS[c])
    ax.invert_yaxis()
    ax.set_xlim(-0.5, xpos[1] + 0.3); ax.set_ylim(len(order) - 0.5, -0.7)
    for i in range(1, len(order)):
        if cat_of[i] != cat_of[i - 1]:
            ax.axhline(i - 0.5, color="#999999", lw=0.9)
    ax.grid(axis="y", color="#f4f4f4", lw=0.6, zorder=0); ax.set_axisbelow(True)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)

    sm = plt.cm.ScalarMappable(cmap="RdBu_r", norm=norm); sm.set_array([])
    cb = fig.colorbar(sm, ax=ax, fraction=0.045, pad=0.03); cb.set_label("NES", fontsize=9)
    handles = [Line2D([], [], marker="o", ls="", markersize=np.sqrt(size(q)) * 0.5,
                      markerfacecolor="#d9d9d9", markeredgecolor="#666666",
                      label=("FDR ≤ 1e-4" if q <= 1e-4 else f"FDR {q:g}"))
               for q in [0.05, 0.001, 1e-4]]
    handles += [Line2D([], [], marker="o", ls="", markersize=8, markerfacecolor="#d9d9d9",
                       markeredgecolor="#111", markeredgewidth=1.3, label="FDR < 0.05")]
    ax.legend(handles=handles, frameon=False, fontsize=8, loc="upper left",
              bbox_to_anchor=(1.28, 1.0), labelspacing=1.2)
    fig.tight_layout()
    for ext in ("png", "svg", "pdf"):
        fig.savefig(out / f"GSEA_shBCAT1_ROS_dotplot.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)

    show = G.pivot(index="gene_set", columns="cell_line", values="NES").reindex(order)
    showf = G.pivot(index="gene_set", columns="cell_line", values="fdr").reindex(order)
    print("\nNES:\n", show.round(2).to_string())
    print("\nFDR:\n", showf.round(3).to_string())
    print("\nwrote:", out)
