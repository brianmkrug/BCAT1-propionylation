#!/usr/bin/env python
"""
Multi-mark CUT&RUN relatedness: PCA + clustered Spearman correlation heatmap.

Sample set (18): H3K4me3 x4, H3K9pr x4, ENL x4, H3K27ac x4, IgG x2.
RNAPII-Ser2P (Nat + Xlink, 5 libraries) is EXCLUDED -- those IPs never cleared
IgG background (1.0-2.0x enrichment; H3K4me3/reports/phase6_findings_report.md).

Method follows analysis_QC_normalization/QC_and_normalization_report.md so the
new figures are directly comparable to the 2026-07-08 pass:
  * 10 kb bins, chr1-22 + chrX, one fragment per read pair
  * log2 CPM
  * PCA on the top 50,000 most variable bins
  * Spearman correlation on all retained bins, average-linkage clustering
Two departures, both deliberate:
  * ENCODE hg38 blacklist v2 bins are dropped BEFORE anything else (the 2026-07-08
    pass ran unmasked and flagged that as a caveat).
  * a depth-equalized PCA is produced alongside the log2-CPM one, because the
    H3K27ac replicates differ ~7-10x in depth and shot noise alone can generate
    an apparent axis of variation at 10 kb.

Usage:  python multimark_pca_corr.py [--outdir DIR]
"""
from __future__ import annotations
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root

import argparse
import gzip
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform
from scipy.stats import rankdata

ROOT = os.path.expanduser(REPO_ROOT + "/CUTRUN")
DEFAULT_OUT = os.path.join(ROOT, "analysis_QC_normalization", "multimark_2026-07")
GRID_NPZ = os.path.join(ROOT, "analysis_QC_normalization", "data", "counts_10kb_matrix.npz")
BLACKLIST = os.path.join(ROOT, "hg38-blacklist.v2.bed.gz")
BIN = 10_000
SEED = 42

# display order = mark blocks, shNT before shBCAT1 within a mark
ORDER = [
    "shNT-H3K4me3-r1", "shNT-H3K4me3-r2", "shBCAT1-H3K4me3-r1", "shBCAT1-H3K4me3-r2",
    "shNT-H3K27ac-r1", "shNT-H3K27ac-r2", "shBCAT1-H3K27ac-r1", "shBCAT1-H3K27ac-r2",
    "shNT-H3K9pr-r1", "shNT-H3K9pr-r2", "shBCAT1-H3K9pr-r1", "shBCAT1-H3K9pr-r2_60M",
    "shNT-ENL-r1", "shNT-ENL-r2", "shBCAT1-ENL-r1", "shBCAT1-ENL-r2",
    "shNT-IgG-r1", "shBCAT1-IgG-r1",
]
MARK_COLOR = {           # distinct hues; IgG deliberately grey
    "H3K4me3": "#2166ac",
    "H3K27ac": "#e8862a",
    "H3K9pr":  "#6c3483",
    "ENL":     "#1e8449",
    "IgG":     "#808080",
}
COND_MARKER = {"shNT": "o", "shBCAT1": "^"}


def mark_of(label: str) -> str:
    for m in MARK_COLOR:
        if f"-{m}-" in label:
            return m
    raise ValueError(label)


def cond_of(label: str) -> str:
    return "shBCAT1" if label.startswith("shBCAT1") else "shNT"


# --------------------------------------------------------------------- load --
def load_grid():
    d = np.load(GRID_NPZ, allow_pickle=True)
    return np.asarray(d["chroms"], dtype="U6"), np.asarray(d["starts"], dtype=np.int64)


def load_counts(outdir, chroms, starts, labels):
    idx = {(c, s): i for i, (c, s) in enumerate(zip(chroms, starts))}
    M = np.zeros((len(chroms), len(labels)), dtype=np.float32)
    off_grid = {}
    for j, lab in enumerate(labels):
        p = os.path.join(outdir, "counts", f"{lab}.10kb.tsv.gz")
        if not os.path.exists(p):
            sys.exit(f"[fatal] missing counts for {lab}: {p}")
        miss = 0
        with gzip.open(p, "rt") as fh:
            for line in fh:
                c, s, n = line.split()
                k = idx.get((c, int(s)))
                if k is None:
                    miss += int(n)
                else:
                    M[k, j] = float(n)
        off_grid[lab] = miss
        print(f"[load] {lab:24s} {M[:, j].sum():>12,.0f} fragments on grid"
              f"  ({miss:,} off-grid)")
    return M, off_grid


def blacklist_mask(chroms, starts):
    """True = keep. A bin is dropped if it overlaps any blacklist interval."""
    keep = np.ones(len(chroms), dtype=bool)
    by_chrom = {}
    for i, (c, s) in enumerate(zip(chroms, starts)):
        by_chrom.setdefault(c, {})[s] = i
    n_reg = 0
    with gzip.open(BLACKLIST, "rt") as fh:
        for line in fh:
            if line.startswith(("#", "track")):
                continue
            f = line.split()
            c, a, b = f[0], int(f[1]), int(f[2])
            if c not in by_chrom:
                continue
            n_reg += 1
            for s in range((a // BIN) * BIN, ((b - 1) // BIN) * BIN + BIN, BIN):
                i = by_chrom[c].get(s)
                if i is not None:
                    keep[i] = False
    print(f"[blacklist] {n_reg} regions on analysed chromosomes -> "
          f"{(~keep).sum():,}/{len(keep):,} bins dropped ({100*(~keep).mean():.2f}%)")
    return keep


# ----------------------------------------------------------------- analysis --
def log2cpm(C):
    tot = C.sum(axis=0, keepdims=True)
    return np.log2(C / tot * 1e6 + 1.0), tot.ravel()


def equalize_depth(C, rng):
    """Binomial thinning of every column to the smallest library size.

    Removes depth-driven shot noise without touching the BAMs. Applied to raw
    counts, so a bin's expected value is unchanged -- only its variance rises to
    the level of the shallowest library, which is what makes the libraries
    comparable rather than the deepest ones artificially smooth.
    """
    tot = C.sum(axis=0)
    target = tot.min()
    out = np.empty_like(C)
    for j in range(C.shape[1]):
        p = target / tot[j]
        out[:, j] = C[:, j] if p >= 0.999 else rng.binomial(C[:, j].astype(np.int64), p)
    print(f"[depth-eq] all libraries thinned to {target:,.0f} fragments")
    return out


def pca(X, ntop=50_000):
    """X = bins x samples of log2 CPM. Returns scores (samples x PC), % var."""
    v = X.var(axis=1)
    top = np.argsort(v)[::-1][:min(ntop, (v > 0).sum())]
    Y = X[top].T                      # samples x bins
    Y = Y - Y.mean(axis=0, keepdims=True)
    U, S, _ = np.linalg.svd(Y, full_matrices=False)
    scores = U * S
    var = S**2 / np.sum(S**2) * 100
    return scores, var, len(top)


def spearman_matrix(X):
    R = np.apply_along_axis(rankdata, 0, X)
    R = R - R.mean(axis=0, keepdims=True)
    R /= np.linalg.norm(R, axis=0, keepdims=True)
    return R.T @ R


# ------------------------------------------------------------------- plots --
def _offsets(x, y, xr, yr):
    """Greedy label placement: try 8 positions around each point, keep the first
    that does not collide with an already-placed label. Cheap stand-in for
    adjustText (not in this env) -- 18 points, so exhaustive is unnecessary."""
    cand = [(7, 4), (-7, 4), (7, -10), (-7, -10), (7, 14), (-7, 14),
            (7, -20), (-7, -20), (7, 24), (-7, 24), (7, -30), (-7, -30)]
    # collision box in data units, sized for a ~20-char label at 6.4 pt
    bw, bh = 0.24 * xr, 0.035 * yr
    placed, out = [], []
    for xi, yi in zip(x, y):
        for dx, dy in cand:
            cx = xi + dx / 100.0 * xr * (1 if dx > 0 else 1)
            cy = yi + dy / 100.0 * yr
            anchor = (cx + (bw / 2 if dx > 0 else -bw / 2), cy)
            if all(abs(anchor[0] - p[0]) > bw or abs(anchor[1] - p[1]) > bh
                   for p in placed):
                placed.append(anchor)
                out.append((dx, dy))
                break
        else:
            placed.append((xi, yi))
            out.append(cand[0])
    return out


def plot_pca(scores, var, labels, path, pcx=0, pcy=1, title=""):
    fig, ax = plt.subplots(figsize=(8.0, 6.4))
    x, y = scores[:, pcx], scores[:, pcy]
    xr, yr = np.ptp(x), np.ptp(y)
    offs = _offsets(x, y, xr, yr)
    for i, lab in enumerate(labels):
        m, c = mark_of(lab), cond_of(lab)
        ax.scatter(x[i], y[i], s=118, zorder=3,
                   color=MARK_COLOR[m], marker=COND_MARKER[c],
                   edgecolor="white", linewidth=1.1)
        dx, dy = offs[i]
        ax.annotate(lab.replace("_60M", ""), (x[i], y[i]),
                    textcoords="offset points", xytext=(dx, dy), fontsize=6.6,
                    ha="left" if dx > 0 else "right", color="#333333", zorder=4)
    ax.margins(0.16)
    ax.axhline(0, color="#dddddd", lw=0.8, zorder=0)
    ax.axvline(0, color="#dddddd", lw=0.8, zorder=0)
    ax.set_xlabel(f"PC{pcx+1} ({var[pcx]:.1f}%)")
    ax.set_ylabel(f"PC{pcy+1} ({var[pcy]:.1f}%)")
    ax.set_title(title, fontsize=10)
    handles = [Line2D([], [], marker="s", ls="", color=MARK_COLOR[m], label=m)
               for m in MARK_COLOR]
    handles += [Line2D([], [], marker=COND_MARKER[c], ls="", color="#444444", label=c)
                for c in COND_MARKER]
    ax.legend(handles=handles, fontsize=7.4, frameon=False, loc="best", ncol=2)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    fig.savefig(path.replace(".png", ".pdf"))
    plt.close(fig)
    print("[plot]", path)


def plot_corr(S, labels, path, title=""):
    D = squareform(np.clip(1 - S, 0, None), checks=False)
    Z = linkage(D, method="average")
    order = dendrogram(Z, no_plot=True)["leaves"]
    Sx = S[np.ix_(order, order)]
    lx = [labels[i].replace("_60M", "") for i in order]

    fig = plt.figure(figsize=(10.4, 9.4))
    # dendrogram directly above the heatmap; colorbar to the RIGHT so it cannot
    # clip the row labels (it did when placed in the left column)
    gs = fig.add_gridspec(2, 2, width_ratios=[40, 1], height_ratios=[1, 5.0],
                          wspace=0.03, hspace=0.02, left=0.20, right=0.93,
                          top=0.93, bottom=0.17)
    axd = fig.add_subplot(gs[0, 0])
    dendrogram(Z, ax=axd, color_threshold=0, above_threshold_color="#666666",
               no_labels=True)
    axd.set_axis_off()
    ax = fig.add_subplot(gs[1, 0])
    im = ax.imshow(Sx, cmap="RdYlBu_r", vmin=np.floor(S.min() * 20) / 20, vmax=1.0)
    ax.set_xticks(range(len(lx)))
    ax.set_xticklabels(lx, rotation=90, fontsize=7)
    ax.set_yticks(range(len(lx)))
    ax.set_yticklabels(lx, fontsize=7)
    for i in range(len(lx)):
        for j in range(len(lx)):
            ax.text(j, i, f"{Sx[i, j]:.2f}", ha="center", va="center", fontsize=4.9,
                    color="white" if Sx[i, j] > 0.80 or Sx[i, j] < 0.35 else "#222222")
    for t, lab in zip(ax.get_xticklabels(), lx):
        t.set_color(MARK_COLOR[mark_of(lab)])
    for t, lab in zip(ax.get_yticklabels(), lx):
        t.set_color(MARK_COLOR[mark_of(lab)])
    cax = fig.add_subplot(gs[1, 1])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label("Spearman r (10 kb bins)", fontsize=8)
    cb.ax.tick_params(labelsize=7)
    fig.suptitle(title, fontsize=10, y=0.975)
    fig.savefig(path, dpi=200)
    fig.savefig(path.replace(".png", ".pdf"))
    plt.close(fig)
    print("[plot]", path)
    return [labels[i] for i in order]


def plot_fingerprint(C, labels, path):
    fig, ax = plt.subplots(figsize=(7.0, 5.6))
    fr = {}
    for j, lab in enumerate(labels):
        x = np.sort(C[:, j])[::-1]
        cum = np.cumsum(x) / x.sum()
        frac = np.arange(1, len(x) + 1) / len(x)
        ax.plot(frac * 100, cum * 100, lw=1.5, color=MARK_COLOR[mark_of(lab)],
                ls="--" if cond_of(lab) == "shBCAT1" else "-", alpha=0.9)
        fr[lab] = cum[max(int(0.01 * len(x)) - 1, 0)] * 100
    ax.set_xscale("log")
    ax.set_xlabel("top fraction of 10 kb bins (%, log scale)")
    ax.set_ylabel("cumulative % of fragments")
    ax.axvline(1, color="#999999", lw=0.8, ls=":")
    ax.set_title("Enrichment fingerprint (blacklist-masked bins)", fontsize=10)
    handles = [Line2D([], [], color=MARK_COLOR[m], label=m) for m in MARK_COLOR]
    handles += [Line2D([], [], color="#444", ls="-", label="shNT"),
                Line2D([], [], color="#444", ls="--", label="shBCAT1")]
    ax.legend(handles=handles, fontsize=7.4, frameon=False, ncol=2)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    fig.savefig(path.replace(".png", ".pdf"))
    plt.close(fig)
    print("[plot]", path)
    return fr


# -------------------------------------------------------------------- main --
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=DEFAULT_OUT)
    a = ap.parse_args()
    plots = os.path.join(a.outdir, "plots")
    data = os.path.join(a.outdir, "data")
    os.makedirs(plots, exist_ok=True)
    os.makedirs(data, exist_ok=True)

    chroms, starts = load_grid()
    C, off_grid = load_counts(a.outdir, chroms, starts, ORDER)
    keep = blacklist_mask(chroms, starts)
    Ck = C[keep]
    print(f"[matrix] {Ck.shape[0]:,} bins x {Ck.shape[1]} samples after masking")

    X, tot = log2cpm(Ck)
    scores, var, ntop = pca(X)
    print("[pca] top %d variable bins; var%% = %s" % (ntop, np.round(var[:5], 1)))
    plot_pca(scores, var, ORDER, os.path.join(plots, "PCA_PC1_PC2.png"), 0, 1,
             "PCA — log2 CPM, 10 kb bins, top 50k variable (Pol2 excluded)")
    plot_pca(scores, var, ORDER, os.path.join(plots, "PCA_PC1_PC3.png"), 0, 2,
             "PCA — log2 CPM, 10 kb bins, top 50k variable (Pol2 excluded)")

    rng = np.random.default_rng(SEED)
    Xe, _ = log2cpm(equalize_depth(Ck, rng))
    se, ve, _ = pca(Xe)
    print("[pca-depth-eq] var%% =", np.round(ve[:5], 1))
    plot_pca(se, ve, ORDER, os.path.join(plots, "PCA_PC1_PC2_depthEqualized.png"), 0, 1,
             "PCA — depth-equalized counts (all libraries thinned to the shallowest)")

    S = spearman_matrix(X)
    order = plot_corr(S, ORDER, os.path.join(plots, "corr_heatmap.png"),
                      "Genome-wide Spearman correlation, 10 kb bins "
                      "(blacklist-masked, Pol2 excluded)")
    fr = plot_fingerprint(Ck, ORDER, os.path.join(plots, "fingerprint_10kb.png"))

    np.savez_compressed(
        os.path.join(data, "multimark_counts_10kb.npz"),
        M=C, labels=np.array(ORDER), chroms=chroms, starts=starts,
        blacklist_keep=keep)
    np.savetxt(os.path.join(data, "spearman_matrix.csv"),
               S, delimiter=",", header=",".join(ORDER), comments="")
    with open(os.path.join(data, "sample_metrics.csv"), "w") as fh:
        fh.write("sample,mark,condition,fragments_chr1_22_X,fragments_after_blacklist,"
                 "pct_signal_top1pct_bins,PC1,PC2,PC3\n")
        for i, lab in enumerate(ORDER):
            fh.write(f"{lab},{mark_of(lab)},{cond_of(lab)},{C[:, i].sum():.0f},"
                     f"{Ck[:, i].sum():.0f},{fr[lab]:.2f},"
                     f"{scores[i,0]:.2f},{scores[i,1]:.2f},{scores[i,2]:.2f}\n")
    with open(os.path.join(data, "corr_cluster_order.txt"), "w") as fh:
        fh.write("\n".join(order) + "\n")
    print("[ok] wrote", data)


if __name__ == "__main__":
    main()
