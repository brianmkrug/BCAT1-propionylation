#!/usr/bin/env python3
"""
figB5 v2 -- compact BCAT1/BCAT2 vs MYC/MYCN panel along the OPC->oligodendrocyte
pseudotime, sized for print at 6 cm (H) x 10 cm (W)  [width doubled 2026-08].

Changes from figB5_MYC.png:
  - dropped the S-phase and G2/M annotation bars (sample age bar retained)
  - dropped PDGFRA / CSPG4 / PLP1 / MBP from the heatmap (MYC, MYCN, BCAT1,
    BCAT2 retained)
  - laid out for a 6 x 5 cm final panel, text scaled accordingly

=============================================================================
METHOD PROVENANCE  (reconstructed 2026-08-05)
=============================================================================
No script was ever committed for figB5_MYC.png / figB6_BCAT_simple.pdf / figB3 --
they were produced in Cowork sessions whose code was not saved. The three
quantifications below were therefore recovered by reverse-engineering the
committed intermediate CSVs and were each confirmed against them numerically.
Do not change them without re-running the checks noted.

(1) EXPRESSION -- pseudobulk CPM per pseudotime bin
        CPM_bin(g) = sum(counts_g over cells in bin)
                     / sum(library sizes of cells in bin) * 1e6
    i.e. cells are POOLED per bin, then normalised. This is NOT the mean of
    per-cell CPM -- that alternative is wrong by up to 8.6 CPM on BCAT1.
    Check: reproduces BCAT1_CPM / BCAT2_CPM in BCAT1_BCAT2_trajectory.csv
           to max abs diff 0.04 CPM.

(2) HEATMAP Z -- z-score of LOG1P(pseudobulk CPM) across the 20 bins
        z = (log1p(CPM_bin) - mean_bins) / sd_bins
    The log1p step matters: z-scoring raw CPM instead is off by up to 0.30 z.
    Check: reproduces BCAT1_Z / BCAT2_Z in BCAT1_BCAT2_trajectory.csv to
           max abs diff 0.005 (vs 0.295 without the log).

(3) SAMPLE AGE -- mean donor age in YEARS per bin, from the Velmeshev
    `development_stage` field (saved as `dev` in human_meta.npz):
        "<n>-year-old stage"   -> n
        "<n>-month-old stage"  -> n/12
        "<x> LMP month stage"  -> (x - 9)/12       prenatal, NEGATIVE years
        "newborn stage"        -> 14/365
        "infant stage"         -> 0.5
    Prenatal samples are negative years relative to birth, which is why bin 0
    has mean_age = -0.056 yr rather than 0.
    Check: reproduces mean_age in bcat_annotation_bins.csv to 0.0085 yr and
           pct_prenatal exactly. (Residual is the newborn/infant midpoint
           convention, which the original did not record.)

Layout follows figB6/figB3: age bar on top, Z heatmap, then CPM line plot.
Binning is 20 equal-frequency quantile bins of the pseudotime `pt`, matching
BCAT1_BCAT2_trajectory.csv and bcat_annotation_bins.csv.
=============================================================================

Inputs (all in this directory):
    hs_pseudotime.npz        pt, genes, celltype   (from hs_pseudotime.py)
    human_counts.npz         raw count matrix, cells x genes
    bcat_annotation_bins.csv per-bin mean sample age (see (3) to regenerate)
Outputs:
    figB5v2_bcat_myc_bins.csv   the plotted per-bin values
    figB5v2_bcat_myc.{svg,pdf,png}
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp

matplotlib.rcParams["svg.fonttype"] = "none"   # editable text in Illustrator
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]

HERE = Path(__file__).resolve().parent
CM = 1 / 2.54

GENES = ["MYC", "MYCN", "BCAT1", "BCAT2"]
COLOR = {                      # BCAT pair = project-canonical Dark2 (see BCAT1_figure_palette.json)
    "BCAT1": "#1b9e77",
    "BCAT2": "#7570b3",
    "MYC":   "#d62728",
    "MYCN":  "#ff7f0e",
}
NBIN_HEAT = 20      # matches the original annotation/heatmap binning
NBIN_LINE = 60      # finer bins for the CPM curves
FS = 5.0            # base font size (pt) for a 6 x 5 cm panel
LW = 0.8

# ----------------------------------------------------------------------------
# Data
# ----------------------------------------------------------------------------
d = np.load(HERE / "hs_pseudotime.npz", allow_pickle=True)
pt = d["pt"]
genes = d["genes"].astype(str)
gidx = {g: i for i, g in enumerate(genes)}

M = sp.load_npz(HERE / "human_counts.npz").tocsc()
lib = np.asarray(M.sum(1)).ravel()
counts = {g: np.asarray(M[:, gidx[g]].todense()).ravel() for g in GENES}


def bin_assign(nbin):
    edges = np.quantile(pt, np.linspace(0, 1, nbin + 1))
    edges[-1] += 1e-9
    return np.clip(np.digitize(pt, edges) - 1, 0, nbin - 1)


def pseudobulk_cpm(gene, b, nbin):
    c = counts[gene]
    return np.array([c[b == i].sum() / lib[b == i].sum() * 1e6 for i in range(nbin)])


b20 = bin_assign(NBIN_HEAT)
b60 = bin_assign(NBIN_LINE)

heat_cpm = {g: pseudobulk_cpm(g, b20, NBIN_HEAT) for g in GENES}
# NOTE: z is taken on log1p(CPM), not raw CPM -- see METHOD PROVENANCE (2).
heat_z = {}
for g, v in heat_cpm.items():
    lv = np.log1p(v)
    heat_z[g] = (lv - lv.mean()) / lv.std()
line_cpm = {g: pseudobulk_cpm(g, b60, NBIN_LINE) for g in GENES}

# mean donor age (years) per bin -- see METHOD PROVENANCE (3)
age = pd.read_csv(HERE / "bcat_annotation_bins.csv")["mean_age"].values

# light smoothing of the fine-binned curves (moving average, width 5 bins)
def smooth(y, w=5):
    k = np.ones(w) / w
    return np.convolve(np.pad(y, w // 2, mode="edge"), k, mode="valid")[:len(y)]


line_s = {g: smooth(v) for g, v in line_cpm.items()}

pd.DataFrame({"bin": np.arange(NBIN_LINE),
              **{f"{g}_CPM": line_cpm[g] for g in GENES},
              **{f"{g}_CPM_smooth": line_s[g] for g in GENES}}
             ).to_csv(HERE / "figB5v2_bcat_myc_bins.csv", index=False)

# self-check: the 20-bin BCAT1/BCAT2 values must still match the committed CSV
_ref = pd.read_csv(HERE / "BCAT1_BCAT2_trajectory.csv")
for _g in ("BCAT1", "BCAT2"):
    _dc = np.abs(heat_cpm[_g] - _ref[f"{_g}_CPM"].values).max()
    _dz = np.abs(heat_z[_g] - _ref[f"{_g}_Z"].values).max()
    assert _dc < 0.1 and _dz < 0.02, f"{_g} drifted from reference: dCPM={_dc}, dZ={_dz}"
    print(f"  check {_g}: max |dCPM|={_dc:.3f}, max |dZ|={_dz:.4f}")

# ----------------------------------------------------------------------------
# Figure -- 6 cm high x 10 cm wide  (width doubled from the original 5 cm)
# ----------------------------------------------------------------------------
fig = plt.figure(figsize=(10 * CM, 6 * CM))

# Explicit axes rectangles rather than a gridspec. With the width doubled,
# L..R spans ~5.2 cm over 20 bins, so each heatmap cell is ~2.60 mm wide;
# ROW_H stays at 1.60 mm (rows are unchanged), giving wide (~1:0.6) cells.
# Font sizes are in points, so text keeps its physical size on the wider panel.
L, R = 0.255, 0.775
ROW_H = 0.16 / 6.0                       # 1.6 mm expressed as a figure fraction
TOP = 0.975
GAP_AH = 0.024                           # age bar -> heatmap
GAP_HL = 0.072                           # heatmap -> line plot

age_y = TOP - ROW_H
heat_y = age_y - GAP_AH - len(GENES) * ROW_H

# ---- sample age bar --------------------------------------------------------
ax_age = fig.add_axes([L, age_y, R - L, ROW_H])
im_age = ax_age.imshow(age[None, :], aspect="auto", cmap="viridis",
                       extent=[0, NBIN_HEAT, 0, 1])
ax_age.set_xticks([])
ax_age.set_yticks([0.5])
ax_age.set_yticklabels(["Age"], fontsize=FS - 0.5)
ax_age.tick_params(length=0, pad=1.5)
for s in ax_age.spines.values():
    s.set_linewidth(0.4)

# ---- 4-gene z-score heatmap ------------------------------------------------
ax_h = fig.add_axes([L, heat_y, R - L, len(GENES) * ROW_H])
Z = np.vstack([heat_z[g] for g in GENES])
im_z = ax_h.imshow(Z, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2,
                   extent=[0, NBIN_HEAT, len(GENES), 0])
ax_h.set_xticks([])
ax_h.set_yticks(np.arange(len(GENES)) + 0.5)
ax_h.set_yticklabels(GENES, fontsize=FS - 0.5)
for t, g in zip(ax_h.get_yticklabels(), GENES):
    t.set_color(COLOR[g])
    t.set_fontweight("bold")
ax_h.tick_params(length=0, pad=1.5)
ax_h.axhline(2, color="black", lw=0.5)          # MYC/MYCN | BCAT1/BCAT2 divider
for s in ax_h.spines.values():
    s.set_linewidth(0.4)

# ---- colour scales for the two annotation blocks ---------------------------
CB_X, CB_W = 0.800, 0.030
CB_FS = FS - 1.2


def add_cbar(mappable, y, h, ticks, ticklabels, title):
    cax = fig.add_axes([CB_X, y, CB_W, h])
    cb = fig.colorbar(mappable, cax=cax, orientation="vertical")
    cb.set_ticks(ticks)
    cb.set_ticklabels(ticklabels)
    cb.ax.tick_params(labelsize=CB_FS, length=1.0, width=0.3, pad=1.0)
    cb.outline.set_linewidth(0.3)
    cax.set_title(title, fontsize=CB_FS, pad=1.6)
    return cb


# stacked in the right margin; positions are explicit so the two never collide
add_cbar(im_age, 0.930, 0.045, [0, 30], ["0", "30"], "yr")
add_cbar(im_z, 0.815, 0.080, [-2, 0, 2], ["−2", "0", "2"], "z")

# ---- CPM curves, BCAT on the left axis, MYC/MYCN on the right --------------
ax = fig.add_axes([L, 0.115, R - L, heat_y - GAP_HL - 0.115])
x = np.linspace(0, 1, NBIN_LINE)

for g in ("BCAT1", "BCAT2"):
    ax.plot(x, line_s[g], color=COLOR[g], lw=LW, solid_capstyle="round")
ax.set_ylabel("BCAT1/2 CPM", fontsize=FS, labelpad=1.5)
ax.set_ylim(0, max(line_s["BCAT1"].max(), line_s["BCAT2"].max()) * 1.30)

axr = ax.twinx()
for g in ("MYC", "MYCN"):
    axr.plot(x, line_s[g], color=COLOR[g], lw=LW, ls="--", dashes=(2.4, 1.6))
axr.set_ylabel("MYC/MYCN CPM", fontsize=FS, labelpad=2.0, color="#555555")
axr.set_ylim(0, max(line_s["MYC"].max(), line_s["MYCN"].max()) * 1.30)

for a in (ax, axr):
    a.tick_params(labelsize=FS - 0.4, length=1.6, width=0.4, pad=1.2)
    for s in a.spines.values():
        s.set_linewidth(0.4)
    a.spines["top"].set_visible(False)
axr.tick_params(colors="#555555")
axr.spines["right"].set_color("#555555")

ax.set_xlim(0, 1)
ax.set_xticks([0, 1])
ax.set_xticklabels(["OPC", "oligo"], fontsize=FS)
ax.get_xticklabels()[0].set_ha("left")
ax.get_xticklabels()[1].set_ha("right")
ax.set_xlabel("pseudotime", fontsize=FS, labelpad=1.5)

handles = [plt.Line2D([], [], color=COLOR[g], lw=LW,
                      ls="--" if g in ("MYC", "MYCN") else "-", label=g)
           for g in ("BCAT1", "BCAT2", "MYC", "MYCN")]
ax.legend(handles=handles, fontsize=FS - 0.5, ncol=2, frameon=False,
          loc="upper right", handlelength=1.3, columnspacing=0.8,
          handletextpad=0.4, labelspacing=0.18, borderaxespad=0.15)

for ext in ("svg", "pdf", "png"):
    p = HERE / f"figB5v2_bcat_myc.{ext}"
    fig.savefig(p, dpi=600 if ext == "png" else None, facecolor="white")
    print("wrote", p)
