#!/usr/bin/env python3
"""Composite figure (16.8 x 13.2 cm):
 TOP half    = 14-panel xenograft BCAA/BCKA (left 2x3) + lipids/glycolysis (right 2x4)
 BOTTOM-left = cell-culture stearoyl-SM + ceramide (2 panels), styled to match the top
 BOTTOM-right= Myelin_composition_nested_pie (vector rasterised at 600 dpi)
Shared K27M/KO styling: red/blue bars, no edge outlines, open-circle points, FST=6,
sci-notation offset moved left of the y-axis, angled x-tick labels."""
import os, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__)); BASE = os.path.dirname(HERE)
XLSX = os.path.join(BASE, "norm_cell_data copy.xlsx"); FIG = os.path.join(HERE, "figures")
MYELIN = os.path.join(FIG, "myelin_raster.png")
GENO_C = {"K27M": "#d1495b", "KO": "#2e86ab"}
XLINES = ["HSJ019", "DIPGXIII"]; CLINES = ["BT245", "DIPGXIII", "HSJ019"]
DROP_LABEL = ("HSJ019", "K27M", 2)
FST = 6
plt.rcParams.update({"axes.facecolor": "white", "figure.facecolor": "white",
                     "axes.edgecolor": "black", "axes.linewidth": 0.35, "axes.grid": False})

# --- top-block panels ---
BCAA = [("Valine", "valine"), ("Leucine", "leucine"), ("Isoleucine", "isoleucine"),
        ("Ketovaline", "3-methyl-2-oxobutyrate"),
        ("Ketoleucine", "4-methyl-2-oxopentanoate"),
        ("Ketoisoleucine", "3-methyl-2-oxovalerate")]
GLYC = [("Tetradec.\n(14:2)", "tetradecadienoate (14:2)*"),
        ("Oleoyl-GPI\n(18:1)", "1-oleoyl-GPI (18:1)"),
        ("Stearoyl-GPI\n(18:0)", "1-stearoyl-GPI (18:0)"),
        ("Linoleoyl-GPI\n(18:2)", "1-linoleoyl-GPI* (18:2)*"),
        ("DHAP", "dihydroxyacetone phosphate (DHAP)"),
        ("Hexose\ndiphosphates", "Isobar: hexose diphosphates"),
        ("PEP", "phosphoenolpyruvate (PEP)"),
        ("3-phospho-\nglycerate", "3-phosphoglycerate")]
# bottom-left cell-culture panels
SMCER = [("Stearoyl SM\n(d18:1/18:0)", "stearoyl sphingomyelin (d18:1/18:0)"),
         ("N-stearoyl-sphingadienine\n(ceramide)", "N-stearoyl-sphingadienine (d18:2/18:0)*")]


def load_xe():
    df = pd.read_excel(XLSX, sheet_name="xenograft"); cnt, keep = {}, []
    for _, r in df.iterrows():
        k = (r["Cell line"], r["Genotype"]); cnt[k] = cnt.get(k, 0) + 1
        keep.append(not (r["Cell line"] == DROP_LABEL[0] and r["Genotype"] == DROP_LABEL[1] and cnt[k] == DROP_LABEL[2]))
    return df.loc[keep].reset_index(drop=True)


def load_cc():
    df = pd.read_excel(XLSX, sheet_name="cell culture")
    return df[df["Cell line"].isin(CLINES)].reset_index(drop=True)


def pfmt(p):
    return f"p={p:.0e}" if p < 1e-2 else f"p={p:.2f}"


def panel(ax, df, title, met, lines, show_x, show_ylab, ylim_mult=1.28, pbrk=1.14, pfs=FST - 2.5, xrot=35):
    vals = pd.to_numeric(df[met], errors="coerce"); xpos, width = np.arange(len(lines)), 0.36
    for gi, geno in enumerate(["K27M", "KO"]):
        for li, line in enumerate(lines):
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
    for li, line in enumerate(lines):
        k = imp[((df["Cell line"] == line) & (df["Genotype"] == "K27M")).values].dropna()
        o = imp[((df["Cell line"] == line) & (df["Genotype"] == "KO")).values].dropna()
        if len(k) >= 2 and len(o) >= 2:
            _, p = stats.ttest_ind(k, o, equal_var=False); yb = ymax * pbrk
            ax.plot([xpos[li] - width/2, xpos[li] + width/2], [yb, yb], color="k", lw=.4)
            ax.text(xpos[li], yb * 1.004, pfmt(p), ha="center", va="bottom", fontsize=pfs)
    ax.set_xticks(xpos)
    ax.set_xticklabels(lines if show_x else [], rotation=xrot,
                       ha=("center" if xrot == 0 else "right"), fontsize=FST)
    ax.set_title(title, fontsize=FST, fontweight="normal", pad=1.5)
    ax.set_ylim(0, ymax * ylim_mult); ax.tick_params(length=1.3, width=0.35, labelsize=FST)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=3))
    ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    if show_ylab:
        ax.set_ylabel("norm. abund.", fontsize=FST)


def relocate_offset(ax):
    ot = ax.yaxis.get_offset_text(); s = ot.get_text(); ot.set_visible(False)
    if s:
        ax.text(-0.02, 1.0, s, transform=ax.transAxes, ha="right", va="bottom", fontsize=FST - 2.5)


def main():
    np.random.seed(7); xe = load_xe(); cc = load_cc()
    fig = plt.figure(figsize=(16.8 / 2.54, 13.2 / 2.54))

    # ---- TOP: 14 panels in ONE evenly-spaced 2x7 grid ----
    # cols 0-2 = BCAA (top) / BCKA (bottom); cols 3-6 = lipids (top) / glycolysis (bottom)
    gtop = fig.add_gridspec(2, 7, left=0.052, right=0.996, top=0.955, bottom=0.585, hspace=0.5, wspace=0.55)
    top_rows = [[BCAA[0], BCAA[1], BCAA[2], GLYC[0], GLYC[1], GLYC[2], GLYC[3]],
                [BCAA[3], BCAA[4], BCAA[5], GLYC[4], GLYC[5], GLYC[6], GLYC[7]]]
    axTop = []
    for i, rowitems in enumerate(top_rows):
        for j, (t, m) in enumerate(rowitems):
            ax = fig.add_subplot(gtop[i, j])
            panel(ax, xe, t, m, XLINES, show_x=(i == 1), show_ylab=(j == 0),
                  ylim_mult=1.34, pbrk=1.14, pfs=5)
            axTop.append(ax)

    # ---- BOTTOM-LEFT: SM + ceramide (same panel height as top; tight y, p at top) ----
    gbL = fig.add_gridspec(1, 2, left=0.052, right=0.62, top=0.42, bottom=0.272, wspace=0.4)
    axB = [fig.add_subplot(gbL[0, j]) for j in range(2)]
    for idx, (ax, (t, m)) in enumerate(zip(axB, SMCER)):
        panel(ax, cc, t, m, CLINES, show_x=True, show_ylab=(idx == 0),
              ylim_mult=1.16, pbrk=1.03, pfs=5, xrot=0)

    fig.canvas.draw()
    for ax in axTop + axB:
        relocate_offset(ax)

    # ---- BOTTOM-RIGHT: myelin nested pie (rasterised vector; ~50% smaller) ----
    img = mpimg.imread(MYELIN); ar = img.shape[1] / img.shape[0]     # w/h
    mw = 0.26                                                        # ~half previous width
    mh = (mw * (16.8 / 13.2)) / ar                                   # keep image aspect
    mx0 = 0.70; my0 = 0.346 - mh / 2                                # centred on the SM/ceramide row band
    axmy = fig.add_axes([mx0, my0, mw, mh]); axmy.imshow(img); axmy.axis("off")

    # ---- one centred K27M/KO legend per data row ----
    handles = [Patch(facecolor=GENO_C["K27M"], label="K27M"), Patch(facecolor=GENO_C["KO"], label="KO")]
    for xc, yc in [(0.524, 0.505), (0.336, 0.19)]:                  # under top row / under SM-ceramide row
        lg = fig.legend(handles=handles, loc="center", ncol=2, frameon=False, fontsize=FST + 1,
                        handlelength=1.2, handleheight=1.2, columnspacing=1.6, bbox_to_anchor=(xc, yc))
        fig.add_artist(lg)

    out = os.path.join(FIG, "composite_xeno_cc_myelin")
    for ext in ("png", "svg", "pdf"):
        fig.savefig(f"{out}.{ext}", dpi=600, bbox_inches="tight")
    print("wrote", out)
    plt.close(fig)


if __name__ == "__main__":
    main()
