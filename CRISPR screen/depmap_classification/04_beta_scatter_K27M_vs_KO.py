#!/usr/bin/env python3
"""
K27M vs KO beta scatter for the DIPGXIII metabolic CRISPR screen.

Regions (thresholds at the top — change them in one place):
  shared dependency  : both betas <= DEP_CUT and |K27M - KO| <= DIAG_BAND   -> bold, dark
  K27M-preferential  : K27M <= DEP_CUT and (K27M - KO) <= -DIAG_BAND        -> red wedge, BCAT1 here
  KO-preferential    : KO   <= DEP_CUT and (KO - K27M) <= -DIAG_BAND        -> blue wedge
  everything else                                                            -> faint grey

Outputs fig_beta_scatter_K27M_vs_KO.{png,pdf,svg} + beta_scatter_regions.csv
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

matplotlib.rcParams["svg.fonttype"] = "none"   # keep SVG text editable in Illustrator

SCR = REPO_ROOT + "/CRISPR screen"
OUT = f"{SCR}/depmap_classification"
PAL = json.load(open(REPO_ROOT + "/BCAT1_figure_palette.json"))

DEP_CUT = -0.50        # a beta this negative or lower counts as depleted
DIAG_BAND = 0.50       # half-width of the "close to 1:1" band, in beta units
N_LABEL = 8            # genes to label in each preferential wedge

C_K27M = PAL["conditions"]["K27M"]["hex"]     # #c0392b
C_KO = PAL["conditions"]["KO"]["hex"]         # #2166ac
C_SHARED = "#2f2f2f"
C_BG = "#d9d9d9"

d = pd.read_excel(f"{SCR}/beta_correlation.xlsx")
x, y = d["K27M.beta"], d["KO.beta"]
delta = x - y

shared = (x <= DEP_CUT) & (y <= DEP_CUT) & (delta.abs() <= DIAG_BAND)
k27m_pref = (x <= DEP_CUT) & (delta <= -DIAG_BAND)
ko_pref = (y <= DEP_CUT) & (delta >= DIAG_BAND)
other = ~(shared | k27m_pref | ko_pref)

d["region"] = np.select([shared, k27m_pref, ko_pref],
                        ["shared dependency", "K27M-preferential", "KO-preferential"],
                        default="neither")
d.sort_values("K27M.beta").to_csv(f"{OUT}/beta_scatter_regions.csv", index=False)

r_p = stats.pearsonr(x, y)
r_s = stats.spearmanr(x, y)

fig, ax = plt.subplots(figsize=(7.6, 7.2))
lo, hi = -2.45, 1.65

# ---- region shading -------------------------------------------------------
# K27M-preferential: x <= DEP_CUT and y >= x + DIAG_BAND  (upper-left wedge)
ax.fill([lo, DEP_CUT, DEP_CUT, lo],
        [lo + DIAG_BAND, DEP_CUT + DIAG_BAND, hi, hi],
        color=C_K27M, alpha=0.06, lw=0, zorder=0)
# KO-preferential: y <= DEP_CUT and x >= y + DIAG_BAND  (lower-right wedge)
ax.fill([lo + DIAG_BAND, hi, hi, DEP_CUT + DIAG_BAND],
        [lo, lo, DEP_CUT, DEP_CUT],
        color=C_KO, alpha=0.06, lw=0, zorder=0)

# 1:1 line and band edges
ax.plot([lo, hi], [lo, hi], ls="--", color="#8c8c8c", lw=1.0, zorder=1)
for off in (-DIAG_BAND, DIAG_BAND):
    ax.plot([lo, hi], [lo + off, hi + off], ls=":", color="#b0b0b0", lw=0.9, zorder=1)
ax.axhline(0, color="#cccccc", lw=0.7, zorder=0)
ax.axvline(0, color="#cccccc", lw=0.7, zorder=0)
ax.axhline(DEP_CUT, color="#cccccc", lw=0.7, ls=(0, (6, 4)), zorder=0)
ax.axvline(DEP_CUT, color="#cccccc", lw=0.7, ls=(0, (6, 4)), zorder=0)

# ---- points ---------------------------------------------------------------
ax.scatter(x[other], y[other], s=9, c=C_BG, edgecolor="none", alpha=0.75, zorder=2)
ax.scatter(x[ko_pref], y[ko_pref], s=22, c=C_KO, edgecolor="none", alpha=0.85, zorder=3)
ax.scatter(x[k27m_pref], y[k27m_pref], s=22, c=C_K27M, edgecolor="none", alpha=0.9, zorder=3)
ax.scatter(x[shared], y[shared], s=34, c=C_SHARED, edgecolor="white", linewidth=0.35, zorder=4)

# ---- labels ---------------------------------------------------------------
def label(sel, sort_col, colour, dx=6, dy=4, n=N_LABEL, weight="normal"):
    sub = d[sel].nsmallest(n, sort_col)
    for _, r in sub.iterrows():
        if r.Gene == "BCAT1":
            continue
        ax.annotate(r.Gene, (r["K27M.beta"], r["KO.beta"]), xytext=(dx, dy),
                    textcoords="offset points", fontsize=7.2, color=colour, fontweight=weight)

label(k27m_pref, "difference", "#8a2620", n=6)
label(ko_pref, "KO.beta", "#17405f", n=5, dx=7, dy=7)
label(shared, "K27M.beta", "#2f2f2f", n=5)

# BCAT1
b = d[d.Gene == "BCAT1"].iloc[0]
ax.scatter([b["K27M.beta"]], [b["KO.beta"]], s=130, facecolor="none",
           edgecolor=C_K27M, linewidth=2.1, zorder=6)
ax.annotate(f"BCAT1\n$\\beta_{{K27M}}$ = {b['K27M.beta']:.2f},  $\\beta_{{KO}}$ = {b['KO.beta']:+.2f}",
            (b["K27M.beta"], b["KO.beta"]), xytext=(52, 40), textcoords="offset points",
            fontsize=10.5, fontweight="bold", color=C_K27M, ha="left", va="center",
            arrowprops=dict(arrowstyle="-", color=C_K27M, lw=1.1, shrinkA=0, shrinkB=9),
            zorder=7)

# ---- region captions ------------------------------------------------------
ax.text(-2.34, -2.26, f"shared\ndependency\nn = {int(shared.sum())}",
        fontsize=9, color=C_SHARED, fontweight="bold", ha="left", va="center",
        linespacing=1.35)
ax.text(-2.34, 1.52, f"K27M-preferential  n = {int(k27m_pref.sum())}",
        fontsize=9, color=C_K27M, fontweight="bold", ha="left", va="center")
ax.text(1.56, -2.32, f"KO-preferential  n = {int(ko_pref.sum())}",
        fontsize=9, color=C_KO, fontweight="bold", ha="right", va="center")

ax.text(0.985, 0.985,
        f"n = {len(d):,} genes\nPearson r = {r_p[0]:.2f}\nSpearman ρ = {r_s[0]:.2f}",
        transform=ax.transAxes, ha="right", va="top", fontsize=8.4, color="#444444")

ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
ax.set_aspect("equal")
ax.set_xlabel("K27M MAGeCK MLE $\\beta$", fontsize=11)
ax.set_ylabel("H3-WT KO MAGeCK MLE $\\beta$", fontsize=11)
ax.set_title("DIPGXIII metabolic CRISPR screen: K27M vs isogenic KO",
             fontsize=12.5, fontweight="bold", pad=10)
ax.spines[["top", "right"]].set_visible(False)
ax.tick_params(labelsize=9)

handles = [
    Line2D([], [], marker="o", ls="", ms=7, mfc=C_SHARED, mec="white",
           label=f"shared dependency (both β ≤ {DEP_CUT:g}, |Δ| ≤ {DIAG_BAND:g})"),
    Line2D([], [], marker="o", ls="", ms=6, mfc=C_K27M, mec="none",
           label=f"K27M-preferential (β$_{{K27M}}$ ≤ {DEP_CUT:g}, Δ ≤ −{DIAG_BAND:g})"),
    Line2D([], [], marker="o", ls="", ms=6, mfc=C_KO, mec="none",
           label=f"KO-preferential (β$_{{KO}}$ ≤ {DEP_CUT:g}, Δ ≥ +{DIAG_BAND:g})"),
    Line2D([], [], ls="--", color="#8c8c8c", label="1:1 (dotted = ±%.1f band)" % DIAG_BAND),
]
ax.legend(handles=handles, fontsize=8.2, loc="upper center", frameon=False,
          bbox_to_anchor=(0.5, -0.10), ncol=2, handlelength=1.4, columnspacing=1.6)

fig.text(0.01, -0.055, "Δ = β$_{K27M}$ − β$_{KO}$.  Thresholds set in 04_beta_scatter_K27M_vs_KO.py.",
         fontsize=7, color="#666666")

fig.tight_layout()
for ext in ("png", "pdf", "svg"):
    fig.savefig(f"{OUT}/fig_beta_scatter_K27M_vs_KO.{ext}",
                dpi=300 if ext == "png" else None, bbox_inches="tight")

print(d.region.value_counts().to_string())
print("\nBCAT1:", b["K27M.beta"], b["KO.beta"], "->", b.region if "region" in b else
      d.loc[d.Gene == "BCAT1", "region"].iloc[0])
