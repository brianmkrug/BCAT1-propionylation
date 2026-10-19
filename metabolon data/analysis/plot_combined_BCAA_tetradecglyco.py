#!/usr/bin/env python3
"""Two-block xenograft panel (14 metabolites, equal-width sub-plots):
 LEFT  block 2x3 : BCAAs (top) with each matching BCKA directly beneath
                   (Valine/Ketovaline, Leucine/Ketoleucine, Isoleucine/Ketoisoleucine)
 RIGHT block 2x4 : 4 lipids (top) over 4 glycolysis-related metabolites (bottom)
Bars have no edge outlines; y-axis sci-notation offset sits left of the axis;
single shared K27M/KO legend. Titles / axis labels / legend share one size."""
import os, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__)); BASE = os.path.dirname(HERE)
XLSX = os.path.join(BASE, "norm_cell_data copy.xlsx"); FIG = os.path.join(HERE, "figures")
GENO_C = {"K27M": "#d1495b", "KO": "#2e86ab"}
LINES = ["HSJ019", "DIPGXIII"]; DROP_LABEL = ("HSJ019", "K27M", 2)
FST = 6   # metabolite-title size (== axis labels & legend)

# LEFT block, row-major 2x3: BCAA row then matching BCKA row
PANELS_L = [("Valine", "valine"), ("Leucine", "leucine"), ("Isoleucine", "isoleucine"),
            ("Ketovaline", "3-methyl-2-oxobutyrate"),
            ("Ketoleucine", "4-methyl-2-oxopentanoate"),
            ("Ketoisoleucine", "3-methyl-2-oxovalerate")]
# RIGHT block, row-major 2x4: lipids row then glycolysis row
PANELS_R = [("Tetradec.\n(14:2)", "tetradecadienoate (14:2)*"),
            ("Oleoyl-GPI\n(18:1)", "1-oleoyl-GPI (18:1)"),
            ("Stearoyl-GPI\n(18:0)", "1-stearoyl-GPI (18:0)"),
            ("Linoleoyl-GPI\n(18:2)", "1-linoleoyl-GPI* (18:2)*"),
            ("DHAP", "dihydroxyacetone phosphate (DHAP)"),
            ("Hexose\ndiphosphates", "Isobar: hexose diphosphates"),
            ("PEP", "phosphoenolpyruvate (PEP)"),
            ("3-phospho-\nglycerate", "3-phosphoglycerate")]

plt.rcParams.update({"axes.facecolor": "white", "figure.facecolor": "white",
                     "axes.edgecolor": "black", "axes.linewidth": 0.35, "axes.grid": False})


def load():
    df = pd.read_excel(XLSX, sheet_name="xenograft"); cnt, keep = {}, []
    for _, r in df.iterrows():
        k = (r["Cell line"], r["Genotype"]); cnt[k] = cnt.get(k, 0) + 1
        keep.append(not (r["Cell line"] == DROP_LABEL[0] and r["Genotype"] == DROP_LABEL[1] and cnt[k] == DROP_LABEL[2]))
    return df.loc[keep].reset_index(drop=True)


def pfmt(p):
    return f"p={p:.0e}" if p < 1e-2 else f"p={p:.2f}"


def panel(ax, df, title, met, show_x, show_ylab):
    vals = pd.to_numeric(df[met], errors="coerce"); xpos, width = np.arange(len(LINES)), 0.36
    for gi, geno in enumerate(["K27M", "KO"]):
        for li, line in enumerate(LINES):
            m = (df["Cell line"] == line) & (df["Genotype"] == geno); v = vals[m.values].dropna().values
            x = xpos[li] + (gi - 0.5) * width
            if len(v):
                ax.bar(x, v.mean(), width=width, color=GENO_C[geno], alpha=.85, edgecolor="none",
                       yerr=(v.std(ddof=1) / np.sqrt(len(v)) if len(v) > 1 else 0),
                       capsize=1.4, error_kw=dict(lw=.5))
                ax.scatter(np.full(len(v), x) + np.random.uniform(-.05, .05, len(v)),
                           v, s=3.5, c="white", edgecolor="k", lw=.35, zorder=5)
    ymax = np.nanmax(vals.values); nz = vals[vals > 0]; hm = (nz.min() / 2.0) if len(nz) else np.nan
    imp = np.log10(vals.replace(0, np.nan).fillna(hm))
    for li, line in enumerate(LINES):
        k = imp[((df["Cell line"] == line) & (df["Genotype"] == "K27M")).values].dropna()
        o = imp[((df["Cell line"] == line) & (df["Genotype"] == "KO")).values].dropna()
        if len(k) >= 2 and len(o) >= 2:
            _, p = stats.ttest_ind(k, o, equal_var=False); yb = ymax * 1.14
            ax.plot([xpos[li] - width/2, xpos[li] + width/2], [yb, yb], color="k", lw=.4)
            ax.text(xpos[li], yb * 1.006, pfmt(p), ha="center", va="bottom", fontsize=FST - 2.5)
    ax.set_xticks(xpos)
    ax.set_xticklabels(LINES if show_x else [], rotation=35, ha="right", fontsize=FST)
    ax.set_title(title, fontsize=FST, fontweight="normal", pad=1.5)
    ax.set_ylim(0, ymax * 1.28); ax.tick_params(length=1.3, width=0.35, labelsize=FST)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=3))
    ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    if show_ylab:
        ax.set_ylabel("norm. abund.", fontsize=FST)


def relocate_offset(ax):
    ot = ax.yaxis.get_offset_text(); s = ot.get_text(); ot.set_visible(False)
    if s:
        ax.text(-0.02, 1.0, s, transform=ax.transAxes, ha="right", va="bottom", fontsize=FST - 2.5)


def main():
    np.random.seed(7); df = load()
    fig = plt.figure(figsize=(16.8 / 2.54, 6.6 / 2.54))
    outer = fig.add_gridspec(1, 2, left=0.055, right=0.997, top=0.92, bottom=0.19,
                             wspace=0.07, width_ratios=[3, 4])
    gsL = outer[0].subgridspec(2, 3, hspace=0.5, wspace=0.55)
    gsR = outer[1].subgridspec(2, 4, hspace=0.5, wspace=0.55)
    axesL = [fig.add_subplot(gsL[i, j]) for i in range(2) for j in range(3)]
    axesR = [fig.add_subplot(gsR[i, j]) for i in range(2) for j in range(4)]

    for idx, (ax, (t, m)) in enumerate(zip(axesL, PANELS_L)):
        panel(ax, df, t, m, show_x=(idx >= 3), show_ylab=(idx % 3 == 0))
    for idx, (ax, (t, m)) in enumerate(zip(axesR, PANELS_R)):
        panel(ax, df, t, m, show_x=(idx >= 4), show_ylab=False)  # left block carries the y-label

    fig.canvas.draw()
    for ax in axesL + axesR:
        relocate_offset(ax)

    fig.legend(handles=[Patch(facecolor=GENO_C["K27M"], label="K27M"),
                        Patch(facecolor=GENO_C["KO"], label="KO")],
               loc="lower center", ncol=2, frameon=False, fontsize=FST,
               handlelength=1.1, handleheight=1.1, columnspacing=1.6,
               bbox_to_anchor=(0.5, 0.0))

    out = os.path.join(FIG, "xenograft_BCAA_tetradecglyco_combined")
    for ext in ("png", "svg", "pdf"):
        fig.savefig(f"{out}.{ext}", dpi=600, bbox_inches="tight")
    print("wrote", out)
    plt.close(fig)


if __name__ == "__main__":
    main()
