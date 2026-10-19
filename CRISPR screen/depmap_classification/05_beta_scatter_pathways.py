#!/usr/bin/env python3
"""
K27M vs KO beta scatter, v2 — Terpenoid backbone biosynthesis and Fatty acid
biosynthesis members highlighted in two colours over the region background.

Same region definitions as 04_beta_scatter_K27M_vs_KO.py.
Outputs fig_beta_scatter_pathways.{png,pdf,svg} + pathway_region_membership.csv
"""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import json
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

matplotlib.rcParams["svg.fonttype"] = "none"

SCR = REPO_ROOT + "/CRISPR screen"
OUT = f"{SCR}/depmap_classification"
PAL = json.load(open(REPO_ROOT + "/BCAT1_figure_palette.json"))

DEP_CUT, DIAG_BAND = -0.50, 0.50
SETS = {
    "Terpenoid backbone biosynthesis": PAL["metabolic_categories"]["Steroid biosynthesis"]["hex"],   # #7b3294
    "Fatty acid biosynthesis":         PAL["pathways"]["Lipid_FA"]["hex"],                            # #1e8449
}
C_K27M = PAL["conditions"]["K27M"]["hex"]
C_KO = PAL["conditions"]["KO"]["hex"]
C_BG = "#dedede"
C_SHARED_BG = "#a8a8a8"

d = pd.read_excel(f"{SCR}/beta_correlation.xlsx")
gmt = pd.read_pickle(f"{SCR}/gsea_out/gmt.pkl")
x, y = d["K27M.beta"], d["KO.beta"]
delta = x - y

shared = (x <= DEP_CUT) & (y <= DEP_CUT) & (delta.abs() <= DIAG_BAND)
k27m_pref = (x <= DEP_CUT) & (delta <= -DIAG_BAND)
ko_pref = (y <= DEP_CUT) & (delta >= DIAG_BAND)
d["region"] = np.select([shared, k27m_pref, ko_pref],
                        ["shared dependency", "K27M-preferential", "KO-preferential"],
                        default="neither")
d["pathway"] = "other"
for name in SETS:
    d.loc[d.Gene.isin(gmt[name]), "pathway"] = name

d[d.pathway != "other"].sort_values(["pathway", "K27M.beta"]).to_csv(
    f"{OUT}/pathway_region_membership.csv", index=False)

fig, ax = plt.subplots(figsize=(8.0, 7.6))
lo, hi = -2.45, 1.65

# region wedges
ax.fill([lo, DEP_CUT, DEP_CUT, lo], [lo + DIAG_BAND, DEP_CUT + DIAG_BAND, hi, hi],
        color=C_K27M, alpha=0.055, lw=0, zorder=0)
ax.fill([lo + DIAG_BAND, hi, hi, DEP_CUT + DIAG_BAND], [lo, lo, DEP_CUT, DEP_CUT],
        color=C_KO, alpha=0.055, lw=0, zorder=0)
ax.plot([lo, hi], [lo, hi], ls="--", color="#8c8c8c", lw=1.0, zorder=1)
for off in (-DIAG_BAND, DIAG_BAND):
    ax.plot([lo, hi], [lo + off, hi + off], ls=":", color="#b8b8b8", lw=0.9, zorder=1)
for v in (0, DEP_CUT):
    ax.axhline(v, color="#cfcfcf", lw=0.7, ls="-" if v == 0 else (0, (6, 4)), zorder=0)
    ax.axvline(v, color="#cfcfcf", lw=0.7, ls="-" if v == 0 else (0, (6, 4)), zorder=0)

# background points
bg = d.pathway == "other"
ax.scatter(x[bg & ~shared], y[bg & ~shared], s=8, c=C_BG, edgecolor="none",
           alpha=0.7, zorder=2)
ax.scatter(x[bg & shared], y[bg & shared], s=17, c=C_SHARED_BG, edgecolor="none",
           alpha=0.9, zorder=2.5)

# pathway points
for name, col in SETS.items():
    m = d.pathway == name
    ax.scatter(x[m], y[m], s=62, c=col, edgecolor="white", linewidth=0.7, zorder=5)

# ---- labels with a simple greedy de-overlap -------------------------------
placed = []          # display-coord boxes already used
CAND = [(9, 6), (9, -12), (-9, 6), (-9, -12), (9, 18), (-9, 18), (9, -24), (-9, -24),
        (22, 6), (-22, 6), (0, 20), (0, -26)]


def collide(bx, boxes):
    return any(not (bx[2] < b[0] or bx[0] > b[2] or bx[3] < b[1] or bx[1] > b[3])
               for b in boxes)


def smart_label(gx, gy, text, colour, bold):
    px, py = ax.transData.transform((gx, gy))
    w, h = 6.4 * len(text) + 6, 15
    for i, (dx, dy) in enumerate(CAND):
        x0, y0 = px + dx - (w if dx < 0 else 0), py + dy - h / 2
        box = (x0, y0, x0 + w, y0 + h)
        if not collide(box, placed):
            placed.append(box)
            arrow = dict(arrowstyle="-", color=colour, lw=0.6, alpha=0.65,
                         shrinkA=0, shrinkB=5) if i > 0 else None
            ax.annotate(text, (gx, gy), xytext=(dx, dy), textcoords="offset points",
                        fontsize=8.0, color=colour, ha="right" if dx < 0 else "left",
                        va="center", fontweight="bold" if bold else "normal",
                        arrowprops=arrow, zorder=6)
            return
    placed.append((px, py, px + w, py + h))
    ax.annotate(text, (gx, gy), xytext=(9, 6), textcoords="offset points",
                fontsize=8.0, color=colour, zorder=6,
                fontweight="bold" if bold else "normal")


lab = d[(d.pathway != "other")].copy()
lab["prio"] = lab.region.map({"K27M-preferential": 0, "KO-preferential": 0,
                              "shared dependency": 1, "neither": 2})
for _, r in lab.sort_values(["prio", "K27M.beta"]).iterrows():
    smart_label(r["K27M.beta"], r["KO.beta"], r.Gene, SETS[r.pathway],
                bold=r.prio == 0)

# BCAT1 reference marker
b = d[d.Gene == "BCAT1"].iloc[0]
ax.scatter([b["K27M.beta"]], [b["KO.beta"]], s=95, facecolor="none",
           edgecolor=C_K27M, linewidth=1.8, zorder=6)
ax.annotate("BCAT1", (b["K27M.beta"], b["KO.beta"]), xytext=(-11, 13),
            textcoords="offset points", fontsize=9.5, fontweight="bold",
            color=C_K27M, ha="right", zorder=6)

# ---- captions -------------------------------------------------------------
ax.text(-2.34, 1.52, f"K27M-preferential  n = {int(k27m_pref.sum())}",
        fontsize=9, color=C_K27M, fontweight="bold", ha="left", va="center")
ax.text(1.56, -2.32, f"KO-preferential  n = {int(ko_pref.sum())}",
        fontsize=9, color=C_KO, fontweight="bold", ha="right", va="center")
ax.text(-2.34, -2.26, f"shared\ndependency\nn = {int(shared.sum())}",
        fontsize=9, color="#5a5a5a", fontweight="bold", ha="left", va="center",
        linespacing=1.35)

r_p = stats.pearsonr(x, y)[0]
ax.text(0.985, 0.985, f"n = {len(d):,} genes\nPearson r = {r_p:.2f}",
        transform=ax.transAxes, ha="right", va="top", fontsize=8.4, color="#444444")

ax.set_xlim(lo, hi); ax.set_ylim(lo, hi); ax.set_aspect("equal")
ax.set_xlabel("K27M MAGeCK MLE $\\beta$", fontsize=11)
ax.set_ylabel("H3-WT KO MAGeCK MLE $\\beta$", fontsize=11)
ax.set_title("Terpenoid backbone vs fatty acid biosynthesis\nin the DIPGXIII K27M / KO screen",
             fontsize=12.5, fontweight="bold", pad=10)
ax.spines[["top", "right"]].set_visible(False)
ax.tick_params(labelsize=9)

counts = {n: d[(d.pathway == n)].region.value_counts() for n in SETS}
handles = [Line2D([], [], marker="o", ls="", ms=8, mfc=c, mec="white",
                  label=f"{n}  ({(d.pathway == n).sum()} genes)")
           for n, c in SETS.items()] + [
    Line2D([], [], marker="o", ls="", ms=5.5, mfc=C_SHARED_BG, mec="none",
           label=f"shared dependency, other genes"),
    Line2D([], [], ls="--", color="#8c8c8c", label="1:1 (dotted = ±%.1f band)" % DIAG_BAND)]
ax.legend(handles=handles, fontsize=8.3, loc="upper center", frameon=False,
          bbox_to_anchor=(0.5, -0.095), ncol=2, handlelength=1.4, columnspacing=1.6)

fig.text(0.01, -0.05,
         "Bold gene labels = member falls in the K27M- or KO-preferential wedge.  "
         "Δ = β$_{K27M}$ − β$_{KO}$; region cut-offs β ≤ −0.5, |Δ| ≶ 0.5.",
         fontsize=7, color="#666666")

fig.tight_layout()
for ext in ("png", "pdf", "svg"):
    fig.savefig(f"{OUT}/fig_beta_scatter_pathways.{ext}",
                dpi=300 if ext == "png" else None, bbox_inches="tight")

for n in SETS:
    s = d[d.pathway == n]
    print(f"\n=== {n} ({len(s)}/{len(gmt[n])} in library)")
    for reg in ("K27M-preferential", "KO-preferential", "shared dependency"):
        sel = s[s.region == reg].sort_values("K27M.beta")
        if len(sel):
            print(f"  {reg}: " + ", ".join(
                f"{r.Gene} (β_K27M {r['K27M.beta']:+.2f}, β_KO {r['KO.beta']:+.2f})"
                for _, r in sel.iterrows()))
