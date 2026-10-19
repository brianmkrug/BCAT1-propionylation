#!/usr/bin/env python3
"""Figure 7G -- BCAA concentrations in DMEM/F-12 versus human plasma.

Input:  media/Media concentration BCAA.xlsx  (Amino acid, Media, Concentration (uM))
Output: BCAA_DMEMF12_vs_plasma.{png,svg,pdf}, BCAA_DMEMF12_vs_plasma_values.csv

Horizontal bars, one pair per BCAA, annotated with the concentration in uM. The
fold-excess of DMEM/F-12 over plasma is printed to stdout and written to the CSV.
Human plasma values are reference concentrations, not measurements from this study.

Run from anywhere:  python media/plot_media_BCAA.py
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


ROOT = repo_root(); HERE = os.path.join(ROOT, "media")
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans"],
                     "font.size": 9, "axes.linewidth": 0.9, "svg.fonttype": "none",
                     "pdf.fonttype": 42, "ps.fonttype": 42})

df = pd.read_excel(os.path.join(HERE, "Media concentration BCAA.xlsx")).dropna(how="all")
V = [c for c in df.columns if "µM" in str(c) or "uM" in str(c)][0]
df[V] = pd.to_numeric(df[V])

ORDER = ["Valine", "Isoleucine", "Leucine"]     # bottom-to-top as in the figure
MEDIA = ["DMEM-F12", "Human plasma"]
C = {"DMEM-F12": "#7b5ea7", "Human plasma": "#f0a020"}

fig, ax = plt.subplots(figsize=(3.4, 2.2))
for i, aa in enumerate(ORDER):
    for off, med in ((0.19, "DMEM-F12"), (-0.19, "Human plasma")):
        v = df[(df["Amino acid"] == aa) & (df["Media"] == med)][V]
        if not len(v):
            continue
        val = float(v.iloc[0])
        ax.barh(i + off, val, height=0.34, color=C[med], edgecolor="none", zorder=2)
        ax.text(val + 8, i + off, f"{val:.0f}", va="center", fontsize=7.5)

ax.set_yticks(range(len(ORDER))); ax.set_yticklabels(ORDER)
ax.set_xlabel("Concentration (µM)")
ax.set_xlim(0, df[V].max() * 1.22)
ax.spines[["top", "right"]].set_visible(False)
ax.tick_params(width=0.9, length=3)
h = [plt.Rectangle((0, 0), 1, 1, fc=C[m], ec="none") for m in MEDIA]
ax.legend(h, MEDIA, frameon=False, fontsize=8, loc="lower right", handlelength=1.1)

for e, kw in (("png", dict(dpi=400)), ("svg", {}), ("pdf", {})):
    fig.savefig(os.path.join(HERE, f"BCAA_DMEMF12_vs_plasma.{e}"),
                bbox_inches="tight", **kw)
print("wrote BCAA_DMEMF12_vs_plasma.{png,svg,pdf}")

w = df.pivot_table(index="Amino acid", columns="Media", values=V)
w["fold_DMEMF12_over_plasma"] = w["DMEM-F12"] / w["Human plasma"]
w.to_csv(os.path.join(HERE, "BCAA_DMEMF12_vs_plasma_values.csv"))
print(); print(w.round(2).to_string())
