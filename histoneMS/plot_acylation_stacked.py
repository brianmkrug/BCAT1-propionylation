#!/usr/bin/env python3
"""Figures 5B and 5C -- histone acyl-PTM stoichiometry by lysine site.

Input:  histoneMS/DIPGXIII PTM.xlsx
          sheet 'shNT means'  Peptide, Site, Modification, DIPGXIII mean   -> Figure 5B
          sheet 'mouse K9'    Peptide, Site, Modification, Mouse brain     -> Figure 5C
Output: DIPGXIII_acylation_stacked_Set2.{png,svg,pdf}
        mouseK9_acylation_stacked_Set2.{png,svg,pdf}
        acylation_stacked_values.csv

Values are relative abundance (fraction of the peptide carrying each modification) from
quantitative histone mass spectrometry. Only the eight ACYL modifications are plotted;
the methyl marks present in the sheet are excluded. Bars are stacked in a fixed acyl
order, coloured with matplotlib Set2. Figure 5B uses a broken x axis because H3K23 and
the H4 tail run an order of magnitude above every other site.

NOTE: the sheet supplies a single mean per site x modification, so no error bars or
replicate points can be drawn -- the underlying per-replicate values are not in this
workbook.

Run from anywhere:  python histoneMS/plot_acylation_stacked.py
"""
import os
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def repo_root(start=None):
    d = os.path.abspath(start or os.path.dirname(__file__))
    while not os.path.exists(os.path.join(d, "BCAT1_figure_palette.json")):
        p = os.path.dirname(d)
        if p == d:
            raise RuntimeError("repository root not found")
        d = p
    return d


ROOT = repo_root(); HERE = os.path.join(ROOT, "histoneMS")
XL = os.path.join(HERE, "DIPGXIII PTM.xlsx")
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans"],
                     "font.size": 8.5, "axes.linewidth": 0.9, "svg.fonttype": "none",
                     "pdf.fonttype": 42, "ps.fonttype": 42})

ACYL = ["acetyl", "propionyl", "butyryl", "hydroxyisobutyryl",
        "succinyl", "crotonyl", "malonyl", "glutaryl"]
COL = dict(zip(ACYL, plt.get_cmap("Set2").colors))
SITES = ["H3K4", "H3K9", "H3K18", "H3K23", "H3K27", "H3K36", "H3K56", "H3K64",
         "H3K79", "H3K122", "H4K5-16 (1ac)", "H4K5-16 (2ac)", "H4K20", "H4K31",
         "H4K79", "H4K91"]


def matrix(sheet, valcol, sites):
    df = pd.read_excel(XL, sheet_name=sheet).dropna(how="all")
    df[valcol] = pd.to_numeric(df[valcol], errors="coerce").fillna(0.0)
    df = df[df.Modification.isin(ACYL)]
    m = (df.pivot_table(index="Site", columns="Modification", values=valcol,
                        aggfunc="sum")
           .reindex(index=sites, columns=ACYL).fillna(0.0))
    return m


def legend(fig, anchor):
    h = [plt.Rectangle((0, 0), 1, 1, fc=COL[a], ec="none") for a in ACYL]
    fig.legend(h, ACYL, loc="center right", bbox_to_anchor=anchor, ncol=2,
               frameon=False, fontsize=7.5, handlelength=1.1, columnspacing=1.0,
               title="PTM", title_fontsize=8)


# ---------------- Figure 5B : DIPGXIII, broken x axis ----------------
m = matrix("shNT means", "DIPGXIII mean", SITES)
BREAK_LO, BREAK_HI = 0.075, 0.245          # x-range hidden by the axis break
fig, (a1, a2) = plt.subplots(1, 2, figsize=(5.6, 4.2), sharey=True,
                             gridspec_kw=dict(width_ratios=[2.4, 1.0], wspace=0.06))
y = np.arange(len(SITES))
for ax in (a1, a2):
    left = np.zeros(len(SITES))
    for a in ACYL:
        ax.barh(y, m[a].values, left=left, height=0.62, color=COL[a],
                edgecolor="none", zorder=2)
        left += m[a].values
a1.set_xlim(0, BREAK_LO); a2.set_xlim(BREAK_HI, 0.405)
a1.spines[["top", "right"]].set_visible(False)
a2.spines[["top", "left"]].set_visible(False)
a2.tick_params(axis="y", length=0)
a1.set_yticks(y); a1.set_yticklabels(SITES); a1.invert_yaxis()
a1.set_ylabel("Histone lysine site")
for ax in (a1, a2):
    ax.tick_params(width=0.9, length=3)
d = 0.012
for ax, xs in ((a1, [1]), (a2, [0])):
    for xx in xs:
        ax.plot([xx - d, xx + d], [-d, +d], transform=ax.transAxes, color="k",
                lw=0.9, clip_on=False)
fig.supxlabel("Relative abundance", fontsize=9, y=0.02)
legend(fig, (1.22, 0.32))
for e, kw in (("png", dict(dpi=400)), ("svg", {}), ("pdf", {})):
    fig.savefig(os.path.join(HERE, f"DIPGXIII_acylation_stacked_Set2.{e}"),
                bbox_inches="tight", **kw)
plt.close(fig); print("wrote DIPGXIII_acylation_stacked_Set2.{png,svg,pdf}")

# ---------------- Figure 5C : mouse brain, H3K9 only ----------------
mk = matrix("mouse K9", "Mouse brain", ["H3K9"])
fig, ax = plt.subplots(figsize=(3.6, 1.35))
left = 0.0
for a in ACYL:
    v = float(mk[a].iloc[0])
    ax.barh(0, v, left=left, height=0.5, color=COL[a], edgecolor="none", zorder=2)
    left += v
ax.set_yticks([0]); ax.set_yticklabels(["H3K9"])
ax.set_xlabel("Relative abundance")
ax.set_title("Mouse brain acylation", fontsize=9, pad=4)
ax.set_xlim(0, left * 1.05)
ax.spines[["top", "right"]].set_visible(False)
ax.tick_params(width=0.9, length=3)
legend(fig, (1.46, 0.5))
for e, kw in (("png", dict(dpi=400)), ("svg", {}), ("pdf", {})):
    fig.savefig(os.path.join(HERE, f"mouseK9_acylation_stacked_Set2.{e}"),
                bbox_inches="tight", **kw)
plt.close(fig); print("wrote mouseK9_acylation_stacked_Set2.{png,svg,pdf}")

out = m.copy(); out.insert(0, "dataset", "DIPGXIII (shNT mean)")
mk2 = mk.copy(); mk2.insert(0, "dataset", "mouse brain")
pd.concat([out, mk2]).to_csv(os.path.join(HERE, "acylation_stacked_values.csv"))
print("\nH3K9  DIPGXIII: acetyl %.4f  propionyl %.4f  (pr/ac = %.1fx)"
      % (m.loc["H3K9", "acetyl"], m.loc["H3K9", "propionyl"],
         m.loc["H3K9", "propionyl"] / m.loc["H3K9", "acetyl"]))
print("H3K9  mouse   : acetyl %.4f  propionyl %.4f  (pr/ac = %.1fx)"
      % (mk.loc["H3K9", "acetyl"], mk.loc["H3K9", "propionyl"],
         mk.loc["H3K9", "propionyl"] / mk.loc["H3K9", "acetyl"]))
print("total H3 propionyl at H3K9 (fraction): %.4f" % m.loc["H3K9", "propionyl"])
