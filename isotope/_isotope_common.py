"""Shared loader and style for the isotope/ acyl-CoA tracing figures.

All panels derive from one workbook:
    Copy of 240911_LuChao_Brian_BCAA_acylCoA_Trace_Peak Area Summary copy.xlsx

The workbook contains no formulas: the raw instrument export (`Sheet1_1`) and the
tidy, plot-ready sheets used here were both produced upstream, and the derivation
between them (natural-abundance correction, molar-enrichment formula) is not
recorded. These scripts therefore reproduce the PLOTTING AND STATISTICS ONLY,
starting from the tidy sheets. See PROVENANCE_AUDIT.md.
"""
import json, os
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

WORKBOOK = ("Copy of 240911_LuChao_Brian_BCAA_acylCoA_Trace_"
            "Peak Area Summary copy.xlsx")


def repo_root(start=None):
    """Walk up from this file until BCAT1_figure_palette.json is found."""
    d = os.path.abspath(start or os.path.dirname(__file__))
    while True:
        if os.path.exists(os.path.join(d, "BCAT1_figure_palette.json")):
            return d
        parent = os.path.dirname(d)
        if parent == d:
            raise RuntimeError("repository root not found "
                               "(no BCAT1_figure_palette.json above %s)" % start)
        d = parent


ROOT = repo_root()
HERE = os.path.join(ROOT, "isotope")
PALETTE = json.load(open(os.path.join(ROOT, "BCAT1_figure_palette.json")))

# Stacked isotopologue panels are black (13C-labelled) over grey (unlabelled).
# This is a figure-specific convention, not a palette entry.
C_LABELED = "#000000"
C_UNLABELED = "#d9d9d9"


def load(sheet):
    """Read one tidy sheet and drop fully-empty rows."""
    df = pd.read_excel(os.path.join(HERE, WORKBOOK), sheet_name=sheet)
    return df.dropna(how="all").reset_index(drop=True)


def set_style():
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans"],
        "font.size": 9, "axes.linewidth": 0.9, "svg.fonttype": "none",
        "pdf.fonttype": 42, "ps.fonttype": 42,
    })


def pfmt(p):
    return f"p={p:.4f}" if p >= 1e-4 else f"p={p:.1e}"


def save(fig, stem):
    out = os.path.join(HERE, stem)
    for ext, kw in (("png", dict(dpi=400)), ("svg", {}), ("pdf", {})):
        fig.savefig(f"{out}.{ext}", bbox_inches="tight", **kw)
    print("wrote", out + ".{png,svg,pdf}")
