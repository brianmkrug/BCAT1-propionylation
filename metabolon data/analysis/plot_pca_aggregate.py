#!/usr/bin/env python3
"""Aggregate PCA figure: cell culture and xenograft side by side, PC1 vs PC2.
Xenograft outlier HSJ019_K27M_2 excluded. Preprocess: >=50% detection filter,
half-min impute, log10, autoscale per metabolite. Colour = genotype, shape = cell line."""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42

from matplotlib.lines import Line2D
from sklearn.decomposition import PCA

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(HERE)
XLSX = os.path.join(BASE, "norm_cell_data copy.xlsx")
FIG = os.path.join(HERE, "figures")
META = ["Cell line", "Genotype", "Source"]
DROP_LABEL = ("HSJ019", "K27M", 2)
GENO_C = {"K27M": "#d1495b", "KO": "#2e86ab", "WT": "#5f8f4e"}
LINE_M = {"HSJ019": "o", "DIPGXIII": "s", "BT245": "^", "hfNSC": "D"}


def load(sheet, drop_outlier=False):
    df = pd.read_excel(XLSX, sheet_name=sheet)
    if drop_outlier:
        cnt, keep = {}, []
        for _, r in df.iterrows():
            k = (r["Cell line"], r["Genotype"]); cnt[k] = cnt.get(k, 0) + 1
            keep.append(not (r["Cell line"] == DROP_LABEL[0] and
                             r["Genotype"] == DROP_LABEL[1] and cnt[k] == DROP_LABEL[2]))
        df = df.loc[keep].reset_index(drop=True)
    meta = df[META].copy()
    X = df.drop(columns=META).apply(pd.to_numeric, errors="coerce")
    return meta, X


def preprocess(X):
    Xz = X.replace(0, np.nan)
    Xz = Xz.loc[:, Xz.notna().mean(axis=0) >= 0.5]
    half = Xz.min(axis=0, skipna=True) / 2.0
    Xl = np.log10(Xz.fillna(half))
    mu, sd = Xl.mean(0), Xl.std(0, ddof=1).replace(0, np.nan)
    return ((Xl - mu) / sd).dropna(axis=1, how="any")


def pca_panel(ax, meta, X, title):
    Xs = preprocess(X)
    p = PCA(n_components=2)
    sc = p.fit_transform(Xs.values)
    ev = p.explained_variance_ratio_ * 100
    for i in range(len(meta)):
        g, cl = meta.iloc[i]["Genotype"], meta.iloc[i]["Cell line"]
        ax.scatter(sc[i, 0], sc[i, 1], c=GENO_C.get(g, "gray"),
                   marker=LINE_M.get(cl, "o"), s=100, edgecolor="k",
                   linewidth=.7, alpha=.9, zorder=3)
    ax.set_xlabel(f"PC1 ({ev[0]:.1f}%)", fontsize=18)
    ax.set_ylabel(f"PC2 ({ev[1]:.1f}%)", fontsize=18)
    ax.set_title(title, fontweight="bold", fontsize=20)
    ax.tick_params(labelsize=15)
    for s in ax.spines.values():          # full box around plot region
        s.set_visible(True)


def main():
    mc, Xc = load("cell culture")
    mx, Xx = load("xenograft", drop_outlier=True)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.4))
    pca_panel(axes[0], mc, Xc, "Cell culture")
    pca_panel(axes[1], mx, Xx, "Xenograft")
    genos = set(mc["Genotype"]) | set(mx["Genotype"])
    lines = set(mc["Cell line"]) | set(mx["Cell line"])
    glabel = {"K27M": "K27M", "KO": "KO", "WT": "hfNSC (WT)"}
    gleg = [Line2D([0], [0], marker="o", color="w", markerfacecolor=GENO_C[g],
                   markeredgecolor="k", markersize=13, label=glabel[g])
            for g in ["K27M", "KO", "WT"] if g in genos]
    lleg = [Line2D([0], [0], marker=LINE_M[l], color="w", markerfacecolor="gray",
                   markeredgecolor="k", markersize=13, label=l)
            for l in ["BT245", "DIPGXIII", "HSJ019", "hfNSC"] if l in lines]
    fig.legend(handles=gleg + lleg, loc="center right", frameon=False,
               fontsize=15, bbox_to_anchor=(1.0, 0.5))
    fig.suptitle("PCA (log10 + autoscale) — PC1 vs PC2", fontweight="bold", fontsize=20)
    fig.tight_layout(rect=[0, 0, 0.88, 0.96])
    for ext in ("png", "svg", "pdf"):
        out = os.path.join(FIG, f"pca_aggregate_PC1PC2.{ext}")
        fig.savefig(out, dpi=200, bbox_inches="tight")
        print("wrote", out)
    plt.close(fig)


if __name__ == "__main__":
    main()
