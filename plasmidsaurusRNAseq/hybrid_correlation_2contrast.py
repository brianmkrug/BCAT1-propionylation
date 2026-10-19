#!/usr/bin/env python3
"""
Second heatmap set for the hybrid data (DIPGXIII shBCAT1 = OLD inducible),
restricted to the shBCAT1 and IleuR contrasts only (4 columns), rendered with
the BrBG diverging palette. Spearman, Pearson and directional-concordance
panels, matching the earlier layout.

Usage:
    python3 hybrid_correlation_2contrast.py --root analysis_excl_outliers
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

warnings.filterwarnings("ignore")

BASE = Path(REPO_ROOT + "/plasmidsaurusRNAseq")
OLD = Path(REPO_ROOT + "/shBCAT1_inducible_analysis")
ALPHA = 0.05
CMAP = "BrBG"                                              # requested palette change
# column order: keep the two lines grouped
COLS = ["BT245\nshBCAT1", "BT245\nIleuR",
        "DIPGXIII\nshBCAT1", "DIPGXIII\nIleuR"]


def new_de(root, cl, test, base):
    d = pd.read_csv(root / "DE_output" / f"DE_{cl}_{test}_vs_{base}.csv", index_col=0)
    d = d.dropna(subset=["gene_name"])
    d.index = d.gene_name.astype(str).str.upper()
    return d[~d.index.duplicated()]


def old_de():
    d = pd.read_csv(OLD / "out" / "deseq2_shBCAT1_vs_shNT.csv", index_col=0)
    d.index = d.index.astype(str).str.upper()
    return d[~d.index.duplicated()]


def build(root):
    spec = {
        "BT245\nshBCAT1":          new_de(root, "BT245", "shBCAT1", "shNT"),
        "BT245\nIleuR":            new_de(root, "BT245", "IleuR", "DMSO"),
        "DIPGXIII\nshBCAT1":       old_de(),
        "DIPGXIII\nIleuR":         new_de(root, "DIPGXIII", "IleuR", "DMSO"),
    }
    lfc = pd.DataFrame({k: v["log2FoldChange"] for k, v in spec.items()})
    padj = pd.DataFrame({k: v["padj"] for k, v in spec.items()})
    return lfc[COLS], padj[COLS]


def heatmap(ax, M, title, cmap, norm, fmt="{:.2f}", annot2=None):
    im = ax.imshow(M.values, cmap=cmap, norm=norm)
    n = len(M)
    ax.set_xticks(range(n)); ax.set_yticks(range(n))
    ax.set_xticklabels(M.columns, rotation=45, ha="right", fontsize=8.5)
    ax.set_yticklabels(M.index, fontsize=8.5)
    for i in range(n):
        for j in range(n):
            v = M.values[i, j]
            if np.isnan(v):
                continue
            lum = abs(norm(v) - 0.5) * 2
            txt = fmt.format(v)
            if annot2 is not None:
                txt += f"\nn={int(annot2.values[i, j])}"
            ax.text(j, i, txt, ha="center", va="center", fontsize=7.5,
                    color="white" if lum > 0.6 else "#222222")
    for k in range(1, n):
        if M.index[k].split("\n")[0] != M.index[k - 1].split("\n")[0]:
            ax.axhline(k - .5, color="white", lw=2.6)
            ax.axvline(k - .5, color="white", lw=2.6)
    ax.set_title(title, fontsize=11, pad=8)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    return im


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--root", required=True)
    a = p.parse_args()
    root = BASE / a.root
    out = root / "COMPARE_inducible" / "hybrid_oldDIPG_shBCAT1"
    out.mkdir(parents=True, exist_ok=True)

    lfc, padj = build(root)
    common = lfc.dropna().index
    lfc_c, padj_c = lfc.loc[common], padj.loc[common]
    union = (padj_c < ALPHA).any(axis=1)
    print(f"common tested in all 4: {len(common)}   union DEGs: {int(union.sum())}")

    sp = lfc_c[union].corr(method="spearman")
    pe = lfc_c.corr(method="pearson")
    cols = lfc_c.columns
    conc = pd.DataFrame(np.nan, index=cols, columns=cols)
    nboth = pd.DataFrame(0, index=cols, columns=cols)
    for x in cols:
        for y in cols:
            m = (padj_c[x] < ALPHA) & (padj_c[y] < ALPHA)
            nboth.loc[x, y] = int(m.sum())
            if m.sum():
                conc.loc[x, y] = 100 * (np.sign(lfc_c.loc[m, x]) ==
                                        np.sign(lfc_c.loc[m, y])).mean()
    sp.to_csv(out / "IleuR_shBCAT1_spearman.csv")
    pe.to_csv(out / "IleuR_shBCAT1_pearson.csv")
    conc.to_csv(out / "IleuR_shBCAT1_concordance.csv")

    norm_c = TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)
    norm_k = TwoSlopeNorm(vmin=0, vcenter=50, vmax=100)
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 4.7))
    im0 = heatmap(axes[0], sp, f"Spearman of log2FC\nunion DEGs (n={int(union.sum())})",
                  CMAP, norm_c)
    im2 = heatmap(axes[1], conc, "Directional concordance (%)\ngenes sig in both",
                  CMAP, norm_k, fmt="{:.0f}%", annot2=nboth)
    for ax, im, lab in zip(axes, [im0, im2],
                           ["Spearman rho", "% same direction"]):
        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
        cb.set_label(lab, fontsize=8); cb.ax.tick_params(labelsize=7)
    fig.tight_layout()
    for ext in ("png", "svg", "pdf"):
        fig.savefig(out / f"hybrid_IleuR_shBCAT1_heatmaps_BrBG_v2.{ext}", dpi=300,
                    bbox_inches="tight")
    plt.close(fig)
    print("\nSpearman:\n", sp.round(2).to_string())
    print("\nwrote:", out / "hybrid_IleuR_shBCAT1_heatmaps_BrBG_v2.png")
