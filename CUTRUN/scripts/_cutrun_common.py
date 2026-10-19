"""Shared helpers for the CUT&RUN manuscript figure scripts."""
import os
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def repo_root(start=None):
    d = os.path.abspath(start or os.path.dirname(__file__))
    while not os.path.exists(os.path.join(d, "BCAT1_figure_palette.json")):
        p = os.path.dirname(d)
        if p == d:
            raise RuntimeError("repository root not found")
        d = p
    return d


ROOT = repo_root()
QC = os.path.join(ROOT, "CUTRUN", "analysis_QC_normalization")
MM10 = os.path.join(ROOT, "CUTRUN", "mm10_analysis")


def set_style(size=9):
    plt.rcParams.update({"font.family": "sans-serif",
                         "font.sans-serif": ["DejaVu Sans"], "font.size": size,
                         "axes.linewidth": 0.9, "svg.fonttype": "none",
                         "pdf.fonttype": 42, "ps.fonttype": 42})


def save(fig, stem):
    for e, kw in (("png", dict(dpi=400)), ("svg", {}), ("pdf", {})):
        fig.savefig(f"{stem}.{e}", bbox_inches="tight", **kw)
    print("wrote", os.path.relpath(stem, ROOT) + ".{png,svg,pdf}")


def bins_log2_vs_igg(npz_path, condition, pseudo=1.0):
    """CPM-normalize the 10 kb bin matrix and return {mark: log2((mark+1)/(IgG+1))}.

    The pseudocount of 1.0 CPM applied to every bin, with NO bin filtering, is what
    reproduces the correlation coefficients printed on Figures 6B and S6A.
    """
    z = np.load(npz_path, allow_pickle=True)
    M = z["M"].astype(float)
    labels = [str(x) for x in z["labels"]]
    cpm = M / M.sum(0) * 1e6
    marks = {}
    for lab, i in zip(labels, range(len(labels))):
        if not lab.startswith(condition + "-"):
            continue
        mark = lab[len(condition) + 1:].split("-")[0]
        marks.setdefault(mark, []).append(i)
    mean = {m: cpm[:, ix].mean(1) for m, ix in marks.items()}
    if "IgG" not in mean:
        raise KeyError("no IgG sample for condition " + condition)
    igg = mean["IgG"]
    return {m: np.log2((v + pseudo) / (igg + pseudo))
            for m, v in mean.items() if m != "IgG"}
