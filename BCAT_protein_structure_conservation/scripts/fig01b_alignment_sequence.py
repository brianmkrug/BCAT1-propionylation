#!/usr/bin/env python3
"""Figure 1b - residue-level BCAT1/BCAT2 alignment, 60 columns per block."""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from Bio.Align import substitution_matrices
from bcat_common import *

aln, sf, cons = load_all()
a1, a2 = aln["a1"], aln["a2"]
m1, m2 = col_maps(a1, a2)
L = len(a1)
M = substitution_matrices.load("BLOSUM62")
NT1, NT2 = 25, 36
PER = 60
nblk = int(np.ceil(L / PER))

FILL_ID, FILL_SIM = "#3c3c44", "#d3d3da"
MONO = {"family": "DejaVu Sans Mono", "size": 6.6}

fig_h = 0.92 * nblk + 1.15
fig, ax = plt.subplots(figsize=(7.2, fig_h))
ax.set_xlim(-9, PER + 6.5); ax.set_ylim(-0.6, nblk * 3.0 + 0.9)
ax.invert_yaxis(); ax.axis("off")

cxxc1 = set(range(SITES["BCAT1"]["cxxc"][0], SITES["BCAT1"]["cxxc"][1] + 1))
cxxc2 = set(range(SITES["BCAT2"]["cxxc"][0], SITES["BCAT2"]["cxxc"][1] + 1))

for b in range(nblk):
    c0, c1 = b * PER, min((b + 1) * PER, L)
    y = b * 3.0
    n1a = next((m1[c] for c in range(c0, c1) if m1[c]), None)
    n1b = next((m1[c] for c in range(c1 - 1, c0 - 1, -1) if m1[c]), None)
    n2a = next((m2[c] for c in range(c0, c1) if m2[c]), None)
    n2b = next((m2[c] for c in range(c1 - 1, c0 - 1, -1) if m2[c]), None)

    for k, c in enumerate(range(c0, c1)):
        x, u, v = k, a1[c], a2[c]
        same = (u == v and u != "-")
        sim = (u != "-" and v != "-" and not same and M[u, v] > 0)
        if same:
            ax.add_patch(Rectangle((x - 0.5, y - 0.45), 1, 1.9, fc=FILL_ID, ec="none"))
        elif sim:
            ax.add_patch(Rectangle((x - 0.5, y - 0.45), 1, 1.9, fc=FILL_SIM, ec="none"))
        # BCAT1-unique N-terminal extension / BCAT2 transit peptide tint
        if m1[c] and m1[c] <= NT1 and not same and not sim:
            ax.add_patch(Rectangle((x - 0.5, y - 0.45), 1, 0.9, fc=C1, alpha=0.25, ec="none"))
        if m2[c] and m2[c] <= NT2 and not same and not sim:
            ax.add_patch(Rectangle((x - 0.5, y + 0.55), 1, 0.9, fc=C2, alpha=0.25, ec="none"))
        ax.text(x, y, u, ha="center", va="center", color="w" if same else "#16161a",
                fontdict=MONO)
        ax.text(x, y + 1.0, v, ha="center", va="center", color="w" if same else "#16161a",
                fontdict=MONO)
        # site markers
        if m1[c] == SITES["BCAT1"]["plp"]:
            ax.text(x, y - 0.72, "★", ha="center", va="center", fontsize=6, color="#111")
        if m2[c] == SITES["BCAT2"]["plp"]:
            ax.text(x, y + 1.78, "★", ha="center", va="center", fontsize=6, color="#111")
        if m1[c] in cxxc1:
            ax.plot([x - 0.45, x + 0.45], [y - 0.62] * 2, color=ACCENT, lw=1.4,
                    solid_capstyle="butt")
        if m2[c] in cxxc2:
            ax.plot([x - 0.45, x + 0.45], [y + 1.68] * 2, color=ACCENT, lw=1.4,
                    solid_capstyle="butt")

    ax.text(-1.6, y, f"{n1a or ''}", ha="right", va="center", fontsize=6, color=C1)
    ax.text(-1.6, y + 1.0, f"{n2a or ''}", ha="right", va="center", fontsize=6, color=C2)
    ax.text(PER - 0.2 + 1.2, y, f"{n1b or ''}", ha="left", va="center", fontsize=6, color=C1)
    ax.text(PER - 0.2 + 1.2, y + 1.0, f"{n2b or ''}", ha="left", va="center", fontsize=6, color=C2)
    ax.text(-8.6, y, "BCAT1", ha="left", va="center", fontsize=6.4, color=C1, fontweight="bold")
    ax.text(-8.6, y + 1.0, "BCAT2", ha="left", va="center", fontsize=6.4, color=C2, fontweight="bold")

# legend (two rows; ~0.92 data units per character at 6 pt)
CW = 0.92
ly = nblk * 3.0 - 0.55
rows = [[(FILL_ID, 1.0, "identical"), (FILL_SIM, 1.0, "similar (BLOSUM62 > 0)")],
        [(C1, 0.25, "BCAT1 N-terminal extension (1–25)"),
         (C2, 0.25, "BCAT2 transit peptide region (1–36)")]]
for r, row in enumerate(rows):
    xx = 0.0
    yy = ly + r * 1.05
    for col, a, lab in row:
        ax.add_patch(Rectangle((xx, yy - 0.34), 1.1, 0.68, fc=col, alpha=a, ec="none"))
        ax.text(xx + 1.7, yy, lab, ha="left", va="center", fontsize=6)
        xx += 1.7 + len(lab) * CW + 2.2
yy = ly + 2.10
ax.text(0, yy, "★", ha="left", va="center", fontsize=6, color="#111")
ax.text(1.4, yy, "PLP-binding lysine", ha="left", va="center", fontsize=6)
x2 = 1.4 + len("PLP-binding lysine") * CW + 2.4
ax.add_patch(Rectangle((x2, yy - 0.09), 1.2, 0.18, fc=ACCENT, ec="none"))
ax.text(x2 + 1.7, yy, "CXXC redox motif", ha="left", va="center", fontsize=6)

ax.set_title("BCAT1 (P54687) vs BCAT2 (O15382) — global pairwise alignment, BCAT1 as reference\n"
             "57.6% identity / 73.7% similarity over 377 aligned positions (BLOSUM62, free end gaps)",
             fontsize=8, pad=14)
save(fig, "Fig1b_BCAT1_BCAT2_sequence_alignment")
