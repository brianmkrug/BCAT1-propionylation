#!/usr/bin/env python3
"""Figure 1 - BCAT1 vs BCAT2 pairwise alignment overview (BCAT1 as reference)."""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrow
from matplotlib.lines import Line2D
from bcat_common import *

aln, sf, cons = load_all()
a1, a2 = aln["a1"], aln["a2"]
s1, s2 = aln["BCAT1"], aln["BCAT2"]
m1, m2 = col_maps(a1, a2)
L = len(a1)
NT1, NT2 = sf["BCAT1"]["idrs"][0][1], sf["BCAT2"]["idrs"][0][1]   # 25, 36

ident = np.array([x == y and x != "-" for x, y in zip(a1, a2)])
from Bio.Align import substitution_matrices
M = substitution_matrices.load("BLOSUM62")
simil = np.array([(x != "-" and y != "-" and M[x, y] > 0 and x != y)
                  for x, y in zip(a1, a2)])
both = np.array([x != "-" and y != "-" for x, y in zip(a1, a2)])
pid = ident.sum() / both.sum()
psim = (ident.sum() + simil.sum()) / both.sum()
xcross = divergent_nterm(a1, a2)
print(f"{both.sum()} aligned columns, {ident.sum()} identical ({pid:.1%}), "
      f"+{simil.sum()} similar ({psim:.1%} total); identity window crosses 40% at col {xcross}")

fig = plt.figure(figsize=(7.2, 5.0))
gs = fig.add_gridspec(3, 1, height_ratios=[1.5, 0.55, 1.0], hspace=0.55,
                      left=0.075, right=0.985, top=0.94, bottom=0.10)

# ---------------------------------------------------------------- A schematic
ax = fig.add_subplot(gs[0]); ax.set_title("a", loc="left", fontweight="bold")
ax.set_xlim(-6, L + 6); ax.set_ylim(-0.10, 1.42); ax.axis("off")
ROWS = {"BCAT1": (0.92, C1, a1, m1, NT1), "BCAT2": (0.28, C2, a2, m2, NT2)}
H = 0.20

# grey band linking the conserved core
core0 = next(i for i, v in enumerate(m1) if v and v > NT1)
ax.add_patch(Rectangle((core0, 0.28), L - core0, 0.84, fc="#f2f2f4",
                       ec="none", zorder=0))
ax.text((core0 + L) / 2, 1.30, "conserved aminotransferase core (fold type IV)",
        ha="center", va="bottom", fontsize=6.8, color=GREY)
ax.annotate("", xy=(core0, 1.265), xytext=(L, 1.265),
            arrowprops=dict(arrowstyle="|-|,widthA=0.25,widthB=0.25", lw=0.7, color=GREY))
ax.annotate("", xy=(core0, 1.27), xytext=(L, 1.27), arrowprops=dict(arrowstyle="-", lw=0))

for g, (y, col, aseq, mp, nt) in ROWS.items():
    present = np.array([c != "-" for c in aseq])
    # backbone
    for start in range(L):
        pass
    ax.add_patch(Rectangle((0, y - H / 2), L, H, fc="none", ec="none"))
    # draw contiguous present blocks
    i = 0
    while i < L:
        if present[i]:
            j = i
            while j < L and present[j]:
                j += 1
            isnt = (mp[i] is not None and mp[i] <= nt)
            ax.add_patch(Rectangle((i, y - H / 2), j - i, H, fc=LIGHT,
                                   ec="#9a9aa2", lw=0.5, zorder=2))
            i = j
        else:
            i += 1
    # divergent N-terminal segment in gene colour
    ntcols = [i for i in range(L) if mp[i] is not None and mp[i] <= nt]
    if ntcols:
        ax.add_patch(Rectangle((ntcols[0], y - H / 2), ntcols[-1] - ntcols[0] + 1, H,
                               fc=col, ec=col, lw=0.5, zorder=3))
    # identity ticks
    for i in np.where(ident)[0]:
        ax.plot([i + 0.5, i + 0.5], [y - H / 2, y + H / 2], color="#7a7a84",
                lw=0.25, alpha=0.55, zorder=4, solid_capstyle="butt")
    ax.text(-4, y, g, ha="right", va="center", fontsize=8, fontweight="bold", color=col)
    ax.text(L + 4, y, f"{SITES[g]['n']} aa", ha="left", va="center", fontsize=6.5, color=GREY)

# gap in BCAT1 opposite the BCAT2 MTS
gap = [i for i in range(L) if a1[i] == "-" and i < 60]
if gap:
    ax.plot([gap[0], gap[-1] + 1], [0.92, 0.92], color="#9a9aa2", lw=0.7, zorder=2)

# site markers
def mark(g, resi, label, colour, dy):
    y = ROWS[g][0]; mp = ROWS[g][3]
    col = next(i for i in range(L) if mp[i] == resi)
    ax.plot([col + 0.5], [y + H / 2 + 0.015], marker="v", ms=3.2, color=colour,
            zorder=6, clip_on=False)
    ax.text(col + 0.5, y + H / 2 + dy, label, ha="center", va="bottom",
            fontsize=6, color=colour)

mark("BCAT1", SITES["BCAT1"]["plp"], "PLP–K222", "#111111", 0.075)
mark("BCAT2", SITES["BCAT2"]["plp"], "PLP–K229", "#111111", 0.075)
mark("BCAT1", SITES["BCAT1"]["cxxc"][0], "CXXC 335–338", ACCENT, 0.075)
c2s = next(i for i in range(L) if m2[i] == SITES["BCAT2"]["cxxc"][0])
ax.plot([c2s + 0.5], [ROWS["BCAT2"][0] - H / 2 - 0.015], marker="^", ms=3.2,
        color=ACCENT, zorder=6)
ax.text(c2s + 0.5, ROWS["BCAT2"][0] - H / 2 - 0.10, "CXXC 342–345", ha="center",
        va="top", fontsize=6, color=ACCENT)
# BCAT2 mitochondrial transit peptide
mtcols = [i for i in range(L) if m2[i] is not None and m2[i] <= 27]
ax.plot([mtcols[0], mtcols[-1] + 1], [0.28 - H / 2 - 0.05] * 2, color=C2, lw=1.4,
        solid_capstyle="butt")
ax.text(mtcols[-1] + 6, 0.28 - H / 2 - 0.05, "mitochondrial transit peptide 1–27 (UniProt)",
        ha="left", va="center", fontsize=6.2, color=C2)
nt1cols = [i for i in range(L) if m1[i] is not None and m1[i] <= NT1]
ax.plot([nt1cols[0], nt1cols[-1] + 1], [0.92 + H / 2 + 0.055] * 2, color=C1, lw=1.4,
        solid_capstyle="butt")
ax.text(nt1cols[-1] + 6, 0.92 + H / 2 + 0.055, "N-terminal extension 1–25  (no BCAT2 counterpart)",
        ha="left", va="center", fontsize=6.2, color=C1)

# ------------------------------------------------------------- B identity bar
ax = fig.add_subplot(gs[1]); ax.set_title("b", loc="left", fontweight="bold")
x = np.arange(L)
ax.bar(x[ident], 1, width=1.0, color="#3c3c44", lw=0)
ax.bar(x[simil], 1, width=1.0, color="#b9b9c2", lw=0)
ax.set_xlim(-6, L + 6); ax.set_ylim(0, 1); ax.set_yticks([])
ax.set_xlabel("alignment column")
ax.spines["left"].set_visible(False)
h = [Rectangle((0, 0), 1, 1, fc="#3c3c44"), Rectangle((0, 0), 1, 1, fc="#b9b9c2")]
ax.legend(h, [f"identical ({pid:.0%})", "similar (BLOSUM62 > 0)"],
          ncol=2, loc="lower right", bbox_to_anchor=(1.0, 1.02), handlelength=1.0,
          handleheight=0.8, columnspacing=1.2)

# ------------------------------------------------ C windowed identity profile
ax = fig.add_subplot(gs[2]); ax.set_title("c", loc="left", fontweight="bold")
prof = identity_profile(a1, a2, 15) * 100
ax.fill_between(x, 0, prof, color="#9a9aa2", lw=0, alpha=0.45)
ax.plot(x, prof, color="#3c3c44", lw=0.8)
ntcols = [i for i in range(L) if m1[i] is not None and m1[i] <= NT1]
mtc = [i for i in range(L) if m2[i] is not None and m2[i] <= 36]
span = (0, max(ntcols[-1], mtc[-1]) + 1)
ax.axvspan(*span, color="#f0c9a8", alpha=0.45, lw=0, zorder=0)
ax.axhline(pid * 100, color=GREY, lw=0.6, ls=(0, (3, 2)))
ax.text(L, pid * 100 - 3, f"whole-protein {pid:.0%}", ha="right", va="top",
        fontsize=6, color=GREY)
ax.text(np.mean(span), 99, "divergent\nN-terminal region", ha="center", va="top",
        fontsize=6.2, color="#a85a12")
ax.set_xlim(-6, L + 6); ax.set_ylim(0, 100)
ax.set_xlabel("alignment column")
ax.set_ylabel("identity\n(15-column window, %)")

fig.suptitle("BCAT1 (P54687) vs BCAT2 (O15382) — pairwise global alignment, BCAT1 as reference",
             y=1.005, fontsize=8.5)
save(fig, "Fig1_BCAT1_BCAT2_alignment_overview")
