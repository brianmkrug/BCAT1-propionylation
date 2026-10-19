"""Shared loader/style for the MDT patient tissue metabolomics figures.

Source: 'PFA_ tissue edited copy.xlsx' (Metabolon tissue panel, 31 samples).
Row 2 of each data sheet is GROUP_DESCRIPTION. The three groups used in the manuscript
map onto the raw labels as:
    Normal        -> Normal brain    (n=5)
    IDH1 MUT GBM  -> IDH-Mut HGG     (n=3)
    WT GBM        -> WT HGG          (n=3)
PF and ST (ependymoma) samples are present in the workbook but excluded from these
comparisons.
"""
import os, json
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

XLSX = "PFA_ tissue edited copy.xlsx"
GROUPS = {"Normal": "Normal brain", "IDH1 MUT GBM": "IDH-Mut HGG", "WT GBM": "WT HGG"}
ORDER = ["Normal brain", "IDH-Mut HGG", "WT HGG"]
COL = {"Normal brain": "#2e8b3d", "IDH-Mut HGG": "#e07b39", "WT HGG": "#7b5ea7"}


def repo_root(start=None):
    d = os.path.abspath(start or os.path.dirname(__file__))
    while not os.path.exists(os.path.join(d, "BCAT1_figure_palette.json")):
        p = os.path.dirname(d)
        if p == d:
            raise RuntimeError("repository root not found")
        d = p
    return d


ROOT = repo_root()
HERE = os.path.join(ROOT, "MDT patient data")


def set_style(size=9):
    plt.rcParams.update({"font.family": "sans-serif",
                         "font.sans-serif": ["DejaVu Sans"], "font.size": size,
                         "axes.linewidth": 0.9, "svg.fonttype": "none",
                         "pdf.fonttype": 42, "ps.fonttype": 42})


def save(fig, path_stem):
    for e, kw in (("png", dict(dpi=400)), ("svg", {}), ("pdf", {})):
        fig.savefig(f"{path_stem}.{e}", bbox_inches="tight", **kw)
    print("wrote", os.path.relpath(path_stem, ROOT) + ".{png,svg,pdf}")


def load_tissue(sheet="MassNormImpData (tissue)"):
    """Return (long dataframe: metabolite, super, sub, sample, group, value)."""
    import openpyxl, warnings
    warnings.filterwarnings("ignore")
    wb = openpyxl.load_workbook(os.path.join(HERE, XLSX), data_only=True)
    rows = list(wb[sheet].iter_rows(values_only=True))
    hdr, grp = rows[0], rows[1]
    cols = [(i, str(hdr[i]), GROUPS.get(str(grp[i])))
            for i in range(3, len(hdr)) if hdr[i] and GROUPS.get(str(grp[i]))]
    rec = []
    for r in rows[3:]:
        if not r or not r[0]:
            continue
        for i, samp, g in cols:
            v = r[i]
            try:
                v = float(v)
            except (TypeError, ValueError):
                continue
            rec.append((str(r[0]).strip(), r[1], r[2], samp, g, v))
    return pd.DataFrame(rec, columns=["metabolite", "super_pathway",
                                      "sub_pathway", "sample", "group", "value"])
