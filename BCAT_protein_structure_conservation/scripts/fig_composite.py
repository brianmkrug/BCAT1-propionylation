#!/usr/bin/env python3
"""Composite publication figure:
  A  residue-level BCAT1/BCAT2 alignment            (was Fig1b)
  B  AlphaFold superposition                        (was Fig2c)
  C  BCAT1 phastCons over the complete CDS          (was Fig3b)
  D  BCAT2 phastCons over the complete CDS          (was Fig3d)
Saved as PDF and SVG (vector, editable text) plus a PNG proof.
"""
import json, textwrap
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.patches import Rectangle
from Bio.Align import substitution_matrices
from bcat_common import *

PANELS = "/root/bc/out/panels"
aln, sf, cons = load_all()
a1, a2 = aln["a1"], aln["a2"]
m1, m2 = col_maps(a1, a2)
L = len(a1)
M = substitution_matrices.load("BLOSUM62")
NT = {"BCAT1": 25, "BCAT2": 36}
COL = {"BCAT1": C1, "BCAT2": C2}
FILL_ID, FILL_SIM = "#3c3c44", "#d3d3da"
MONO = {"family": "DejaVu Sans Mono", "size": 6.3}
PER, NBLK = 60, int(np.ceil(L / 60))
BLOCK_H = 2.62                     # vertical units per alignment block


# --------------------------------------------------------------- A: alignment
def draw_alignment(ax):
    ax.set_xlim(-9.5, PER + 7.0)
    ax.set_ylim(-0.9, NBLK * BLOCK_H + 1.9)
    ax.invert_yaxis(); ax.axis("off")
    cx1 = set(range(*[SITES["BCAT1"]["cxxc"][0], SITES["BCAT1"]["cxxc"][1] + 1]))
    cx2 = set(range(*[SITES["BCAT2"]["cxxc"][0], SITES["BCAT2"]["cxxc"][1] + 1]))

    for b in range(NBLK):
        c0, c1 = b * PER, min((b + 1) * PER, L)
        y = b * BLOCK_H
        n1a = next((m1[c] for c in range(c0, c1) if m1[c]), None)
        n1b = next((m1[c] for c in range(c1 - 1, c0 - 1, -1) if m1[c]), None)
        n2a = next((m2[c] for c in range(c0, c1) if m2[c]), None)
        n2b = next((m2[c] for c in range(c1 - 1, c0 - 1, -1) if m2[c]), None)

        for k, c in enumerate(range(c0, c1)):
            x, u, v = k, a1[c], a2[c]
            same = (u == v and u != "-")
            sim = (u != "-" and v != "-" and not same and M[u, v] > 0)
            if same:
                ax.add_patch(Rectangle((x - .5, y - .42), 1, 1.84, fc=FILL_ID, ec="none"))
            elif sim:
                ax.add_patch(Rectangle((x - .5, y - .42), 1, 1.84, fc=FILL_SIM, ec="none"))
            if m1[c] and m1[c] <= NT["BCAT1"] and not same and not sim:
                ax.add_patch(Rectangle((x - .5, y - .42), 1, .88, fc=C1, alpha=.25, ec="none"))
            if m2[c] and m2[c] <= NT["BCAT2"] and not same and not sim:
                ax.add_patch(Rectangle((x - .5, y + .54), 1, .88, fc=C2, alpha=.25, ec="none"))
            ax.text(x, y, u, ha="center", va="center",
                    color="w" if same else "#16161a", fontdict=MONO)
            ax.text(x, y + 1.0, v, ha="center", va="center",
                    color="w" if same else "#16161a", fontdict=MONO)
            if m1[c] == SITES["BCAT1"]["plp"]:
                ax.text(x, y - 0.70, "★", ha="center", va="center", fontsize=5.6, color="#111")
            if m2[c] == SITES["BCAT2"]["plp"]:
                ax.text(x, y + 1.72, "★", ha="center", va="center", fontsize=5.6, color="#111")
            if m1[c] in cx1:
                ax.add_patch(Rectangle((x - .45, y - .66), .9, .10, fc=ACCENT, ec="none"))
            if m2[c] in cx2:
                ax.add_patch(Rectangle((x - .45, y + 1.58), .9, .10, fc=ACCENT, ec="none"))

        ax.text(-1.8, y, f"{n1a or ''}", ha="right", va="center", fontsize=5.8, color=C1)
        ax.text(-1.8, y + 1.0, f"{n2a or ''}", ha="right", va="center", fontsize=5.8, color=C2)
        ax.text(PER + 1.0, y, f"{n1b or ''}", ha="left", va="center", fontsize=5.8, color=C1)
        ax.text(PER + 1.0, y + 1.0, f"{n2b or ''}", ha="left", va="center", fontsize=5.8, color=C2)
        ax.text(-9.2, y, "BCAT1", ha="left", va="center", fontsize=6.1, color=C1, fontweight="bold")
        ax.text(-9.2, y + 1.0, "BCAT2", ha="left", va="center", fontsize=6.1, color=C2, fontweight="bold")

    # legend, two rows
    CW = 0.92
    ly = NBLK * BLOCK_H - 0.75
    rows = [[(FILL_ID, 1.0, "identical"), (FILL_SIM, 1.0, "similar (BLOSUM62 > 0)")],
            [(C1, .25, "BCAT1 N-terminal extension (1–25)"),
             (C2, .25, "BCAT2 transit-peptide region (1–36)")]]
    for r, row in enumerate(rows):
        xx, yy = 0.0, ly + r * 0.95
        for col, a, lab in row:
            ax.add_patch(Rectangle((xx, yy - .32), 1.05, .64, fc=col, alpha=a, ec="none"))
            ax.text(xx + 1.65, yy, lab, ha="left", va="center", fontsize=5.9)
            xx += 1.65 + len(lab) * CW + 2.2
    yy = ly + 1.95
    ax.text(0, yy, "★", ha="left", va="center", fontsize=5.9, color="#111")
    ax.text(1.4, yy, "PLP-binding lysine", ha="left", va="center", fontsize=5.9)
    x2 = 1.4 + len("PLP-binding lysine") * CW + 2.4
    ax.add_patch(Rectangle((x2, yy - .09), 1.15, .18, fc=ACCENT, ec="none"))
    ax.text(x2 + 1.75, yy, "CXXC redox motif", ha="left", va="center", fontsize=5.9)


# ------------------------------------------------------------- B: structure
def draw_structure(ax):
    im = mpimg.imread(f"{PANELS}/3d_superposition.png")
    mask = im[..., :3].min(axis=2) < 0.985
    ys, xs = np.where(mask)
    pad = 10
    im = im[max(ys.min() - pad, 0):ys.max() + pad, max(xs.min() - pad, 0):xs.max() + pad]
    ax.imshow(im); ax.axis("off")
    ax.text(0.92, 0.97, "BCAT1\nN-term", transform=ax.transAxes, ha="center", va="top",
            fontsize=6.0, color=C1, fontweight="bold")
    ax.text(0.09, 0.80, "BCAT2\nN-term", transform=ax.transAxes, ha="center", va="top",
            fontsize=6.0, color=C2, fontweight="bold")
    ax.text(0.5, 0.015, "Cα RMSD 0.35 Å / 235 residues", transform=ax.transAxes,
            ha="center", va="bottom", fontsize=5.9, color=GREY)


# ------------------------------------------------------------- C,D: phastCons
def draw_phastcons(ax, g, xlabel=True):
    c = cons[g]
    n = len(c["protein"])
    pc = np.array(c["phastCons"], float)[:n * 3].reshape(n, 3).mean(axis=1)
    x = np.arange(1, n + 1)
    col = COL[g]
    sizes = [b - a for a, b in c["blocks"]]
    if c["strand"] == "-":
        sizes = sizes[::-1]
    edges = np.concatenate([[0], np.cumsum(sizes) / 3.0])

    ax.axvspan(1, NT[g], color=col, alpha=0.14, lw=0, zorder=0)
    ax.fill_between(x, 0, pc, color="#3c3c44", lw=0, alpha=0.85, step="mid", zorder=2)
    for e in edges[1:-1]:
        ax.axvline(e, color="#b9b9c2", lw=0.4, ls=(0, (2, 2)), zorder=1)
    for i in range(len(sizes)):
        ax.add_patch(Rectangle((edges[i], 1.22), edges[i + 1] - edges[i], 0.13,
                               fc=col if i % 2 == 0 else "#c9ccd1", ec="w", lw=0.4,
                               clip_on=False, zorder=5))
    ax.text(-0.014, 1.285, "CDS\nexons", transform=ax.get_yaxis_transform(),
            ha="right", va="center", fontsize=5.6, color=GREY)

    plp = SITES[g]["plp"]; cx = SITES[g]["cxxc"]
    ax.plot([plp], [1.06], marker="v", ms=3.0, color="#111", clip_on=False, zorder=6)
    ax.text(plp, 1.10, f"K{plp}", ha="center", va="bottom", fontsize=5.6, clip_on=False)
    ax.add_patch(Rectangle((cx[0], 1.045), cx[1] - cx[0], 0.035, fc=ACCENT, ec="none",
                           clip_on=False, zorder=6))
    ax.text(np.mean(cx), 1.10, f"CXXC", ha="center", va="bottom", fontsize=5.6,
            color=ACCENT, clip_on=False)

    mn, mc = pc[:NT[g]].mean(), pc[NT[g]:].mean()
    ax.annotate(f"N-term 1–{NT[g]}, mean {mn:.2f}", xy=(NT[g] / 2, 0.30),
                xytext=(NT[g] + 30, 0.30), fontsize=5.8, color=col, va="center",
                ha="left", arrowprops=dict(arrowstyle="->", lw=0.6, color=col),
                bbox=dict(fc="w", ec="none", alpha=0.85, pad=0.8), zorder=8)
    ax.text(0.995, 0.86, f"core mean {mc:.2f}", transform=ax.transAxes, fontsize=5.8,
            color=GREY, ha="right", va="top",
            bbox=dict(fc="w", ec="none", alpha=0.85, pad=1.0), zorder=8)
    ax.set_xlim(0, 400); ax.set_ylim(0, 1.04)
    ax.set_yticks([0, 0.5, 1.0])
    ax.set_ylabel("phastCons", fontsize=6.6)
    ax.spines["bottom"].set_bounds(0, n)
    if xlabel:
        ax.set_xlabel("codon position along the complete CDS (residue number)", fontsize=6.6)
    ax.text(0.995, 1.40, f"{g}  ({c['transcript']}, {c['chrom']}, {c['strand']} strand, "
            f"{n * 3 + 3} coding bases incl. stop)", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=6.0, color=col)


# ------------------------------------------------------------------- layout
FIG_W, FIG_H = 7.2, 10.3
fig = plt.figure(figsize=(FIG_W, FIG_H))

A_TOP, A_BOT = 0.982, 0.475          # alignment panel
axA = fig.add_axes([0.075, A_BOT, 0.915, A_TOP - A_BOT]); draw_alignment(axA)

axB = fig.add_axes([0.032, 0.172, 0.295, 0.225]); draw_structure(axB)
axC = fig.add_axes([0.445, 0.298, 0.535, 0.076]); draw_phastcons(axC, "BCAT1", xlabel=False)
axD = fig.add_axes([0.445, 0.160, 0.535, 0.076]); draw_phastcons(axD, "BCAT2", xlabel=True)

for ax, lab, dx, dy in [(axA, "A", -0.062, 0.004), (axB, "B", -0.010, 0.015),
                        (axC, "C", -0.093, 0.40), (axD, "D", -0.093, 0.40)]:
    ax.text(dx, 1.0 + dy, lab, transform=ax.transAxes, ha="left", va="bottom",
            fontsize=10, fontweight="bold")

# --------------------------------------------------------------- captions
CAP = dict(fontsize=5.9, color="#1a1a1a", va="top", ha="left", linespacing=1.45)

capA = ("**A** Global pairwise alignment of human BCAT1 (UniProt P54687, 386 aa) and BCAT2 "
        "(O15382, 392 aa), BCAT1 as reference. Needleman–Wunsch with BLOSUM62, gap open −11 / "
        "extend −1 and free end gaps: 217 of 377 aligned positions are identical (57.6%), 73.7% "
        "similar. Dark boxes, identical residues; light boxes, similar (BLOSUM62 > 0). Tinted "
        "residues mark the BCAT1 N-terminal extension (1–25) and the BCAT2 transit-peptide "
        "region (1–36); neither has a counterpart in the other paralogue, and windowed identity "
        "does not exceed 40% until alignment column 57. The catalytic PLP-binding lysine "
        "(BCAT1 K222 / BCAT2 K229) and the CXXC redox motif (BCAT1 C335–C338 / BCAT2 C342–C345) "
        "are aligned.")

capB = ("**B** AlphaFold DB models AF-P54687-F1 and AF-O15382-F1 (v6) superimposed on "
        "structurally equivalent Cα atoms (iterative outlier removal): RMSD 0.35 Å over 235 "
        "residues. Folded cores, grey; BCAT1 1–25, green; BCAT2 1–30, purple. Both N-terminal "
        "tails are low-confidence predictions (pLDDT < 70) and are called intrinsically "
        "disordered (metapredict V3: BCAT1 1–25, BCAT2 1–36).")

capCD = ("**C, D** phastCons conservation (UCSC hg38, 100-way vertebrate alignment) across the "
         "complete coding sequence of each gene, shown as the per-codon mean of the three base "
         "scores against residue number. Dashed lines and the block track give CDS exon "
         "boundaries; shading marks the divergent N-terminal region. Mean phastCons in the "
         "N-terminal region versus the folded core is 0.31 vs 0.80 for BCAT1 and 0.64 vs 0.69 "
         "for BCAT2; the corresponding phyloP difference is significant for both "
         "(Mann–Whitney p = 1.3 × 10⁻¹¹ and 3.6 × 10⁻⁶). Spliced CDS translates exactly to the "
         "UniProt canonical sequence for both transcripts; phastCons covers 100% of coding bases.")


def emit(x, y, text, width, **kw):
    """Render **bold** lead-ins by splitting the first token off."""
    wrapped = textwrap.fill(text.replace("**", ""), width)
    fig.text(x, y, wrapped, **kw)
    lead = text.split("**")[1]
    fig.text(x, y, lead, fontweight="bold", **kw)


emit(0.075, 0.466, capA, 150, **CAP)
emit(0.032, 0.154, capB, 46, **CAP)
emit(0.445, 0.116, capCD, 88, **CAP)

fig.savefig("/root/bc/out/Fig_composite_BCAT1_BCAT2.pdf")
fig.savefig("/root/bc/out/Fig_composite_BCAT1_BCAT2.svg")
fig.savefig("/root/bc/out/Fig_composite_BCAT1_BCAT2.png", dpi=450)
print("wrote Fig_composite_BCAT1_BCAT2 .pdf / .svg / .png")
plt.close(fig)
