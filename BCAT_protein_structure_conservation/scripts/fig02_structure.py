#!/usr/bin/env python3
"""Figure 2 - AlphaFold models of BCAT1/BCAT2 with a 2D per-residue rendering
of secondary structure, model confidence and predicted disorder."""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.patches import Rectangle, FancyArrow, Patch
from bcat_common import *

P = "/root/bc/out/panels"
aln, sf, cons = load_all()


def crop(path, pad=12):
    im = mpimg.imread(path)
    rgb = im[..., :3]
    mask = (rgb.min(axis=2) < 0.985)
    ys, xs = np.where(mask)
    y0, y1 = max(ys.min() - pad, 0), min(ys.max() + pad, im.shape[0])
    x0, x1 = max(xs.min() - pad, 0), min(xs.max() + pad, im.shape[1])
    return im[y0:y1, x0:x1]


def ss_track(ax, ss, y, h, colour):
    """Draw helices as blocks, strands as arrows, coil as a line."""
    n = len(ss)
    ax.plot([1, n], [y, y], color="#8d8d96", lw=0.8, solid_capstyle="butt", zorder=2)
    i = 0
    while i < n:
        c = ss[i]
        j = i
        while j < n and ss[j] == c:
            j += 1
        if c == "H":
            ax.add_patch(Rectangle((i + 1, y - h / 2), j - i, h, fc=colour,
                                   ec="none", zorder=3))
        elif c == "E":
            L = j - i
            head = min(L * 0.45, 6)
            ax.add_patch(FancyArrow(i + 1, y, L - head, 0, width=h * 0.62,
                                    head_width=h * 1.06, head_length=head,
                                    length_includes_head=False, fc="#6f6f7a",
                                    ec="none", zorder=3))
        i = j


fig = plt.figure(figsize=(7.2, 7.4))
gs = fig.add_gridspec(3, 3, height_ratios=[1.45, 1.0, 1.0],
                      hspace=0.42, wspace=0.06,
                      left=0.085, right=0.985, top=0.955, bottom=0.055)

# ------------------------------------------------------------- 3D renders
titles = [("3d_bcat1_plddt.png", "BCAT1 (AF-P54687-F1)", "a"),
          ("3d_bcat2_plddt.png", "BCAT2 (AF-O15382-F1)", "b"),
          ("3d_superposition.png", "superposition", "c")]
for k, (f, t, lab) in enumerate(titles):
    ax = fig.add_subplot(gs[0, k])
    ax.imshow(crop(f"{P}/{f}"))
    ax.axis("off")
    ax.set_title(lab, loc="left", fontweight="bold", y=0.98)
    bb = ax.get_position()
    fig.text((bb.x0 + bb.x1) / 2, 0.652, t, ha="center", va="top", fontsize=7)
    if lab == "c":
        ax.text(0.90, 0.955, "BCAT1\nN-term", transform=ax.transAxes, ha="center",
                va="top", fontsize=6.2, color=C1, fontweight="bold")
        ax.text(0.10, 0.80, "BCAT2\nN-term", transform=ax.transAxes, ha="center",
                va="top", fontsize=6.2, color=C2, fontweight="bold")

# pLDDT legend under a/b
ax = fig.axes[0]
hs = [Patch(fc=c, ec="none") for _, _, c, _ in PLDDT_BANDS]
ls = [l for _, _, _, l in PLDDT_BANDS]
fig.legend(hs, ls, ncol=4, loc="upper center", bbox_to_anchor=(0.38, 0.639),
           handlelength=1.0, handleheight=0.8, columnspacing=1.1,
           title="AlphaFold pLDDT", title_fontsize=6.5)
fig.text(0.86, 0.637, "core Cα RMSD 0.35 Å (235 residues)",
         ha="center", va="top", fontsize=6, color=GREY)

# ------------------------------------------------------------- 2D tracks
for r, g in enumerate(["BCAT1", "BCAT2"]):
    d = sf[g]
    n = len(d["seq"])
    x = np.arange(1, n + 1)
    colour = C1 if g == "BCAT1" else C2
    nt = d["idrs"][0][1]

    ax = fig.add_subplot(gs[r + 1, :])
    ax.set_title("de"[r], loc="left", fontweight="bold")

    # disordered / low-confidence N-terminal segment
    ax.axvspan(1, nt, color=colour, alpha=0.12, lw=0, zorder=0)
    ax.text(nt + 3, 1.02, f"predicted IDR 1–{nt}", fontsize=6, color=colour,
            va="bottom", ha="left")

    # pLDDT as coloured bars + line
    pl = np.array(d["plddt"]) / 100
    for lo, hi, col, _ in PLDDT_BANDS:
        m = (np.array(d["plddt"]) > lo) & (np.array(d["plddt"]) <= hi)
        ax.bar(x[m], pl[m], width=1.0, bottom=0, color=col, lw=0, alpha=0.6, zorder=1)
    ax.plot(x, pl, color="#2b2b33", lw=0.6, zorder=4)

    # disorder scores
    ax.plot(x, d["metapredict"], color="#c0392b", lw=1.0, zorder=5,
            label="metapredict V3 disorder")
    ax.axhline(0.5, color=GREY, lw=0.5, ls=":", zorder=3)
    fi = np.array(d["foldindex"])
    ax.plot(x, fi + 0.62, color="#2e86c1", lw=0.95, ls=(0, (3, 1.5)), zorder=6)
    ax.axhline(0.62, color="#2e86c1", lw=0.45, ls=":", alpha=0.8, zorder=3)
    axf = ax.secondary_yaxis("right", functions=(lambda v: v - 0.62,
                                                 lambda f: f + 0.62))
    axf.set_ylabel("FoldIndex", color="#2e86c1", labelpad=1)
    axf.set_yticks([-0.5, -0.25, 0.0, 0.25])
    axf.tick_params(axis="y", colors="#2e86c1")
    axf.spines["right"].set_color("#2e86c1")

    # secondary structure track below the axes
    ss_track(ax, d["ss"], -0.165, 0.085, colour)
    ax.text(-0.012, -0.165, "SS", transform=ax.get_yaxis_transform(), ha="right",
            va="center", fontsize=6, color=GREY)

    # key residues
    plp = SITES[g]["plp"]; cx = SITES[g]["cxxc"]
    ax.plot([plp], [1.08], marker="v", ms=3.4, color="#111", clip_on=False, zorder=6)
    ax.text(plp, 1.12, f"PLP–K{plp}", ha="center", va="bottom", fontsize=6, clip_on=False)
    ax.plot([cx[0], cx[1]], [1.08] * 2, color=ACCENT, lw=1.6, clip_on=False, zorder=6)
    ax.text(np.mean(cx), 1.12, f"CXXC {cx[0]}–{cx[1]}", ha="center", va="bottom",
            fontsize=6, color=ACCENT, clip_on=False)
    if g == "BCAT2":
        ax.plot([1, 27], [-0.055] * 2, color=C2, lw=1.8, clip_on=False, zorder=6)
        ax.text(30, -0.055, "transit peptide (UniProt 1–27)", fontsize=6, color=C2,
                va="center", ha="left")

    ax.set_xlim(0, 400); ax.set_ylim(-0.235, 1.0)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_yticklabels(["0", "0.25", "0.5", "0.75", "1"])
    ax.set_ylabel("pLDDT / 100,\ndisorder score")
    ax.set_xlabel(f"{g} residue")
    ax.spines["bottom"].set_bounds(0, n)
    ax.spines["left"].set_bounds(0, 1)
    ax.text(0.998, 1.04, g, transform=ax.transAxes, ha="right", va="bottom",
            fontsize=8, fontweight="bold", color=colour)
    if r == 0:
        from matplotlib.lines import Line2D
        hh = [Line2D([], [], color="#c0392b", lw=1.0),
              Line2D([], [], color="#2e86c1", lw=0.9, ls=(0, (3, 1.5))),
              Line2D([], [], color="#2b2b33", lw=0.6)]
        ax.legend(hh, ["metapredict V3 disorder", "FoldIndex (<0 = disordered)",
                       "pLDDT / 100"], ncol=3, loc="lower center",
                  bbox_to_anchor=(0.54, 0.055), handlelength=1.6, fontsize=6)

fig.suptitle("AlphaFold DB models and per-residue confidence / disorder for BCAT1 and BCAT2",
             y=0.995, fontsize=8.5)
save(fig, "Fig2_BCAT_structure_disorder")
