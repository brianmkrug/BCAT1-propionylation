#!/usr/bin/env python3
"""Figure S5B -- absolute abundance (iBAQ) of BCAA-catabolic and acyl-CoA enzymes.

Input:  histoneMS/iBAQ.xlsx  (Protein, Description, r1, r2, r3, Group)
Output: iBAQ_BCAA_AcylCoA_barplot.{png,svg,pdf}, iBAQ_BCAA_AcylCoA_values.csv

Horizontal bars of the mean iBAQ across three replicates, points overlaid, sorted
within each group. The two groups ('BCAA catabolism', 'Acyl-CoA metabolism') are drawn
as separate blocks with a gap between them; BCAT1 is bolded.

Run from anywhere:  python histoneMS/plot_iBAQ_BCAA_AcylCoA.py
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
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans"],
                     "font.size": 8.5, "axes.linewidth": 0.9, "svg.fonttype": "none",
                     "pdf.fonttype": 42, "ps.fonttype": 42})

REPS = ["r1", "r2", "r3"]
GCOL = {"BCAA catabolism": "#3fa48c", "Acyl-CoA metabolism": "#7b6bb0"}

df = pd.read_excel(os.path.join(HERE, "iBAQ.xlsx")).dropna(how="all")
for r in REPS:
    df[r] = pd.to_numeric(df[r], errors="coerce")
df["mean"] = df[REPS].mean(axis=1)
df["sem"] = df[REPS].std(axis=1, ddof=1) / np.sqrt(len(REPS))

blocks, ypos, y = [], [], 0.0
for g in ["Acyl-CoA metabolism", "BCAA catabolism"]:      # bottom block first
    sub = df[df.Group == g].sort_values("mean")
    for _, r in sub.iterrows():
        ypos.append((y, r)); y += 1
    blocks.append((g, y)); y += 1.4                        # gap between blocks

fig, ax = plt.subplots(figsize=(3.6, 0.19 * len(ypos) + 1.2))
rng = np.random.default_rng(5)
for yy, r in ypos:
    ax.barh(yy, r["mean"], height=0.68, color=GCOL[r["Group"]], edgecolor="none",
            xerr=r["sem"], error_kw=dict(elinewidth=0.8, capsize=2, capthick=0.8),
            zorder=2)
    ax.scatter([r[c] for c in REPS], yy + rng.uniform(-0.14, 0.14, len(REPS)),
               s=9, facecolor="white", edgecolor="black", linewidth=0.5, zorder=4)

ax.set_yticks([yy for yy, _ in ypos])
ax.set_yticklabels([r["Protein"] for _, r in ypos], fontsize=7)
for t in ax.get_yticklabels():
    if t.get_text() == "BCAT1":
        t.set_fontweight("bold")
ax.set_xlabel("iBAQ")
ax.set_ylim(-1, y - 1.0)
ax.spines[["top", "right"]].set_visible(False)
ax.tick_params(width=0.9, length=3)
for g, top in blocks:
    idx = [yy for yy, r in ypos if r["Group"] == g]
    ax.text(1.02, (min(idx) + max(idx)) / 2, g, transform=ax.get_yaxis_transform(),
            rotation=270, va="center", ha="left", fontsize=8, color=GCOL[g])

for e, kw in (("png", dict(dpi=400)), ("svg", {}), ("pdf", {})):
    fig.savefig(os.path.join(HERE, f"iBAQ_BCAA_AcylCoA_barplot.{e}"),
                bbox_inches="tight", **kw)
print("wrote iBAQ_BCAA_AcylCoA_barplot.{png,svg,pdf}")
df[["Protein", "Description", "Group"] + REPS + ["mean", "sem"]] \
  .sort_values(["Group", "mean"], ascending=[True, False]) \
  .to_csv(os.path.join(HERE, "iBAQ_BCAA_AcylCoA_values.csv"), index=False)
print("\nn proteins: %d  (%s)" % (len(df), dict(df.Group.value_counts())))
print(df.nlargest(6, "mean")[["Protein", "Group", "mean"]].to_string(index=False))
print("BCAT1 rank within BCAA catabolism: %d of %d"
      % ((df[df.Group == "BCAA catabolism"]["mean"] >
          float(df[df.Protein == "BCAT1"]["mean"].iloc[0])).sum() + 1,
         (df.Group == "BCAA catabolism").sum()))
