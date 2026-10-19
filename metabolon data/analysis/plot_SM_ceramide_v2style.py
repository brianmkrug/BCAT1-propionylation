#!/usr/bin/env python3
"""Stearoyl SM (d18:1/18:0) + best K27M-depleted ceramide, cell culture (3 lines),
styled to match cellculture_BCAA_BCKA_bars_v2 (full box, open-circle points, compact
p-labels, sci y-offset, bottom K27M/KO legend). Figure sized to a single row of that
3-panel-per-row grid. Exports png/svg/pdf."""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42

from matplotlib.patches import Patch
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(HERE)
XLSX = os.path.join(BASE, "norm_cell_data copy.xlsx")
FIG = os.path.join(HERE, "figures")
GENO_C = {"K27M": "#d1495b", "KO": "#2e86ab"}
LINES = ["BT245", "DIPGXIII", "HSJ019"]

PANELS = [
    ("Stearoyl sphingomyelin (d18:1/18:0)", "stearoyl sphingomyelin (d18:1/18:0)"),
    ("N-stearoyl-sphingadienine (d18:2/18:0)", "N-stearoyl-sphingadienine (d18:2/18:0)*"),
]


def load():
    df = pd.read_excel(XLSX, sheet_name="cell culture")
    return df[df["Cell line"].isin(LINES)].reset_index(drop=True)


def pfmt(p):
    return f"p={p:.0e}" if p < 0.01 else f"p={p:.2f}"


def panel(ax, df, title, met, ylab):
    vals = pd.to_numeric(df[met], errors="coerce")
    xpos, w = np.arange(len(LINES)), 0.40
    for gi, geno in enumerate(["K27M", "KO"]):
        for li, line in enumerate(LINES):
            m = (df["Cell line"] == line) & (df["Genotype"] == geno)
            v = vals[m.values].dropna().values
            x = xpos[li] + (gi - 0.5) * w
            if len(v):
                ax.bar(x, v.mean(), width=w, color=GENO_C[geno], edgecolor="k",
                       linewidth=1.0, zorder=2)
                ax.errorbar(x, v.mean(),
                            yerr=(v.std(ddof=1) / np.sqrt(len(v)) if len(v) > 1 else 0),
                            color="k", elinewidth=1.2, capsize=4, capthick=1.2, zorder=4)
                ax.scatter(np.full(len(v), x) + np.random.uniform(-.07, .07, len(v)),
                           v, s=42, facecolors="white", edgecolors="k",
                           linewidths=1.0, zorder=5)
    ymax = np.nanmax(vals.values)
    nz = vals[vals > 0]; hm = (nz.min() / 2.0) if len(nz) else np.nan
    imp = np.log10(vals.replace(0, np.nan).fillna(hm))
    for li, line in enumerate(LINES):
        k = imp[((df["Cell line"] == line) & (df["Genotype"] == "K27M")).values].dropna()
        o = imp[((df["Cell line"] == line) & (df["Genotype"] == "KO")).values].dropna()
        if len(k) >= 2 and len(o) >= 2:
            _, p = stats.ttest_ind(k, o, equal_var=False)
            yb = ymax * 1.10
            ax.plot([xpos[li] - w/2, xpos[li] + w/2], [yb, yb], color="k", lw=1.0)
            ax.text(xpos[li], yb * 1.015, pfmt(p), ha="center", va="bottom", fontsize=11)
    ax.set_xticks(xpos); ax.set_xticklabels(LINES, fontsize=12)
    ax.set_title(title, fontsize=12, fontweight="bold")
    if ylab:
        ax.set_ylabel("norm. abund.", fontsize=12)
    ax.set_ylim(0, ymax * 1.28)
    ax.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))
    ax.yaxis.get_offset_text().set_fontsize(11)
    ax.tick_params(labelsize=11)
    for s in ax.spines.values():
        s.set_visible(True); s.set_linewidth(1.0)


def main():
    np.random.seed(7)
    df = load()
    # single-row footprint of the v2 grid: full width, one-row height
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 4.7))
    for i, (ax, (title, met)) in enumerate(zip(axes, PANELS)):
        panel(ax, df, title, met, ylab=(i == 0))
    fig.legend(handles=[Patch(facecolor=GENO_C["K27M"], edgecolor="k", label="K27M"),
                        Patch(facecolor=GENO_C["KO"], edgecolor="k", label="KO")],
               loc="lower center", ncol=2, frameon=False, fontsize=15,
               handletextpad=0.5, columnspacing=1.6, bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=[0, 0.09, 1, 1])
    for ext in ("png", "svg", "pdf"):
        out = os.path.join(FIG, f"cellculture_SM_ceramide_bars_v2.{ext}")
        fig.savefig(out, dpi=200, bbox_inches="tight")
        print("wrote", out)
    plt.close(fig)


if __name__ == "__main__":
    main()
