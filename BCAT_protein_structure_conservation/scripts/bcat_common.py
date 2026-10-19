"""Shared style + data loading for the BCAT1/BCAT2 sequence-structure figures."""
import json, os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

D = "/root/bc/data"
OUT = "/root/bc/out"
os.makedirs(OUT, exist_ok=True)

# project palette (BCAT1_figure_palette.json v1.8, "genes" block)
C1, C2 = "#1b9e77", "#7570b3"          # BCAT1, BCAT2
ACCENT = "#d95f02"                      # CXXC redox motif
GREY = "#4d4d4d"
LIGHT = "#e8e8ea"

# AlphaFold pLDDT bands
PLDDT_BANDS = [(90, 100, "#0053D6", "Very high (>90)"),
               (70, 90,  "#65CBF3", "Confident (70–90)"),
               (50, 70,  "#FFDB13", "Low (50–70)"),
               (0,  50,  "#FF7D45", "Very low (<50)")]

mpl.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "font.size": 7, "axes.labelsize": 7.5, "axes.titlesize": 8.5,
    "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6.5,
    "axes.linewidth": 0.7, "xtick.major.width": 0.7, "ytick.major.width": 0.7,
    "xtick.major.size": 2.6, "ytick.major.size": 2.6,
    "axes.spines.top": False, "axes.spines.right": False,
    "legend.frameon": False, "figure.dpi": 150,
    "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
})

# ---------------------------------------------------------------- annotations
SITES = {
    "BCAT1": dict(plp=222, cxxc=(335, 338), uniprot="P54687", n=386,
                  tx="NM_005504.7", chrom="chr12", strand="-"),
    "BCAT2": dict(plp=229, cxxc=(342, 345), uniprot="O15382", n=392,
                  mts=(1, 27), tx="NM_001190.4", chrom="chr19", strand="-"),
}


def load_all():
    aln = json.load(open(f"{D}/aln.json"))
    sf = json.load(open(f"{D}/struct_features.json"))
    cons = json.load(open(f"{D}/conservation.json"))
    return aln, sf, cons


def col_maps(a1, a2):
    """alignment column -> 1-based residue index (or None) for each sequence."""
    m1, m2, i, j = [], [], 0, 0
    for x, y in zip(a1, a2):
        if x != "-":
            i += 1; m1.append(i)
        else:
            m1.append(None)
        if y != "-":
            j += 1; m2.append(j)
        else:
            m2.append(None)
    return m1, m2


def identity_profile(a1, a2, win=15):
    """Fraction identical over a sliding window of alignment columns."""
    eq = np.array([1.0 if (x == y and x != "-") else 0.0 for x, y in zip(a1, a2)])
    k = np.ones(win) / win
    return np.convolve(eq, k, mode="same")


def divergent_nterm(a1, a2, win=15, thresh=0.40):
    """First alignment column from which windowed identity stays >= thresh."""
    p = identity_profile(a1, a2, win)
    for c in range(len(p)):
        if p[c] >= thresh:
            return c
    return 0


def per_codon(values, n_aa):
    """Collapse per-base CDS values (transcript order) to per-codon means."""
    v = np.array(values[:n_aa * 3], dtype=float).reshape(n_aa, 3)
    return v.mean(axis=1), v


def save(fig, stem, also_png=True):
    fig.savefig(f"{OUT}/{stem}.pdf")
    if also_png:
        fig.savefig(f"{OUT}/{stem}.png", dpi=450)
    print("wrote", stem + ".pdf/.png")
    plt.close(fig)
