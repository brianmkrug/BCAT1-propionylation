#!/usr/bin/env python
"""
Replicate-pooled, spike-in-scaled profiles: correlation heatmap + occupancy view.

Pooling: replicates are SUMMED at the count level (not averaged), so the pooled
profile uses all available read depth, which is the point. Ten profiles:
H3K4me3 / H3K27ac / H3K9pr / ENL x {shNT, shBCAT1}, plus one IgG per condition.

Spike-in scaling: applied to the shBCAT1 side after CPM, so the tracks carry
genome-wide occupancy change rather than composition only. CPM forces every
library to sum to 1e6 and therefore *destroys* exactly the abundance information
the spike-in restores -- which is why the factor is applied after, not instead of,
CPM.

    ENL      shBCAT1 x 0.80
    H3K9pr   shBCAT1 x 0.62
    H3K27ac  shBCAT1 x 1.00   (invariant)
    H3K4me3  shBCAT1 x 1.00   ASSUMED -- no factor was supplied for H3K4me3
    IgG      x 1.00           control, not spike-normalized

IMPORTANT -- what the correlation heatmap can and cannot show:
Spearman is invariant to any monotone rescaling of a sample, and Pearson to any
affine one. Multiplying shBCAT1 by 0.80 changes NEITHER. The correlation matrix
below is therefore numerically identical with and without spike-in scaling, and
this script asserts that rather than leaving it implicit. To see the occupancy
change you need a scale-sensitive view, so the log2-ratio panel is produced
alongside: there, the spike-in shifts each mark's whole distribution to
log2(factor), and bins departing from that centre are the ones moving beyond the
genome-wide trend.
"""
from __future__ import annotations
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root

import gzip
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform
from scipy.stats import rankdata

ROOT = os.path.expanduser(REPO_ROOT + "/CUTRUN")
OUT = os.path.join(ROOT, "analysis_QC_normalization", "multimark_2026-07")
GRID = os.path.join(ROOT, "analysis_QC_normalization", "data", "counts_10kb_matrix.npz")
BLACKLIST = os.path.join(ROOT, "hg38-blacklist.v2.bed.gz")
BIN = 10_000

SPIKEIN = {"ENL": 0.80, "H3K9pr": 0.62, "H3K27ac": 1.00, "H3K4me3": 1.00, "IgG": 1.00}
ASSUMED = {"H3K4me3"}          # no factor supplied; flagged on the figure

POOLS = {
    "H3K4me3":  {"shNT": ["shNT-H3K4me3-r1", "shNT-H3K4me3-r2"],
                 "shBCAT1": ["shBCAT1-H3K4me3-r1", "shBCAT1-H3K4me3-r2"]},
    "H3K27ac":  {"shNT": ["shNT-H3K27ac-r1", "shNT-H3K27ac-r2"],
                 "shBCAT1": ["shBCAT1-H3K27ac-r1", "shBCAT1-H3K27ac-r2"]},
    "H3K9pr":   {"shNT": ["shNT-H3K9pr-r1", "shNT-H3K9pr-r2"],
                 "shBCAT1": ["shBCAT1-H3K9pr-r1", "shBCAT1-H3K9pr-r2_60M"]},
    "ENL":      {"shNT": ["shNT-ENL-r1", "shNT-ENL-r2"],
                 "shBCAT1": ["shBCAT1-ENL-r1", "shBCAT1-ENL-r2"]},
    "IgG":      {"shNT": ["shNT-IgG-r1"], "shBCAT1": ["shBCAT1-IgG-r1"]},
}
ORDER = [f"{c}-{m}" for m in ("H3K4me3", "H3K27ac", "H3K9pr", "ENL", "IgG")
         for c in ("shNT", "shBCAT1")]
MARK_COLOR = {"H3K4me3": "#2166ac", "H3K27ac": "#e8862a", "H3K9pr": "#6c3483",
              "ENL": "#1e8449", "IgG": "#808080"}


def mark_of(p):
    return p.split("-", 1)[1]


def cond_of(p):
    return p.split("-", 1)[0]


def sp_matrix(X):
    R = np.apply_along_axis(rankdata, 0, X)
    R = R - R.mean(axis=0, keepdims=True)
    R /= np.linalg.norm(R, axis=0, keepdims=True)
    return R.T @ R


def load_grid_and_counts():
    d = np.load(GRID, allow_pickle=True)
    chroms = np.asarray(d["chroms"], dtype="U6")
    starts = np.asarray(d["starts"], dtype=np.int64)
    idx = {(c, s): i for i, (c, s) in enumerate(zip(chroms, starts))}

    pooled = np.zeros((len(chroms), len(ORDER)), dtype=np.float64)
    depths = {}
    for j, prof in enumerate(ORDER):
        m, c = mark_of(prof), cond_of(prof)
        tot = 0
        for lab in POOLS[m][c]:
            p = os.path.join(OUT, "counts", f"{lab}.10kb.tsv.gz")
            with gzip.open(p, "rt") as fh:
                for line in fh:
                    ch, st, n = line.split()
                    k = idx.get((ch, int(st)))
                    if k is not None:
                        pooled[k, j] += float(n)
                        tot += int(n)
        depths[prof] = tot
        print(f"[pool] {prof:20s} {'+'.join(POOLS[m][c]):45s} = {tot:>12,} fragments")
    return chroms, starts, pooled, depths


def blacklist_mask(chroms, starts):
    keep = np.ones(len(chroms), dtype=bool)
    by = {}
    for i, (c, s) in enumerate(zip(chroms, starts)):
        by.setdefault(c, {})[s] = i
    with gzip.open(BLACKLIST, "rt") as fh:
        for line in fh:
            if line.startswith(("#", "track")):
                continue
            f = line.split()
            c, a, b = f[0], int(f[1]), int(f[2])
            if c not in by:
                continue
            for s in range((a // BIN) * BIN, ((b - 1) // BIN) * BIN + BIN, BIN):
                i = by[c].get(s)
                if i is not None:
                    keep[i] = False
    print(f"[blacklist] {(~keep).sum():,}/{len(keep):,} bins dropped "
          f"({100*(~keep).mean():.2f}%)")
    return keep


def plot_corr(S, labels, path, title, subtitle):
    D = squareform(np.clip(1 - S, 0, None), checks=False)
    Z = linkage(D, method="average")
    order = dendrogram(Z, no_plot=True)["leaves"]
    Sx = S[np.ix_(order, order)]
    lx = [labels[i] for i in order]

    fig = plt.figure(figsize=(8.8, 8.4))
    gs = fig.add_gridspec(2, 2, width_ratios=[40, 1], height_ratios=[1, 5.0],
                          wspace=0.03, hspace=0.02, left=0.24, right=0.92,
                          top=0.90, bottom=0.20)
    axd = fig.add_subplot(gs[0, 0])
    dendrogram(Z, ax=axd, color_threshold=0, above_threshold_color="#666666",
               no_labels=True)
    axd.set_axis_off()
    ax = fig.add_subplot(gs[1, 0])
    im = ax.imshow(Sx, cmap="RdYlBu_r", vmin=np.floor(S.min() * 20) / 20, vmax=1.0)
    ax.set_xticks(range(len(lx)))
    ax.set_xticklabels(lx, rotation=90, fontsize=8.4)
    ax.set_yticks(range(len(lx)))
    ax.set_yticklabels(lx, fontsize=8.4)
    for i in range(len(lx)):
        for j in range(len(lx)):
            ax.text(j, i, f"{Sx[i, j]:.2f}", ha="center", va="center", fontsize=7,
                    color="white" if Sx[i, j] > 0.80 or Sx[i, j] < 0.35 else "#222222")
    for t, lab in zip(ax.get_xticklabels(), lx):
        t.set_color(MARK_COLOR[mark_of(lab)])
    for t, lab in zip(ax.get_yticklabels(), lx):
        t.set_color(MARK_COLOR[mark_of(lab)])
    cax = fig.add_subplot(gs[1, 1])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label("Spearman r (10 kb bins)", fontsize=8)
    cb.ax.tick_params(labelsize=7)
    fig.suptitle(title, fontsize=11, y=0.975)
    fig.text(0.5, 0.935, subtitle, ha="center", fontsize=8, color="#555555")
    fig.savefig(path, dpi=200)
    fig.savefig(path.replace(".png", ".pdf"))
    plt.close(fig)
    print("[plot]", path)
    return lx


def plot_occupancy(cpm, keep_hi, path):
    """The scale-SENSITIVE view the correlation matrix cannot give."""
    marks = ["H3K4me3", "H3K27ac", "H3K9pr", "ENL"]
    fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.9),
                             gridspec_kw={"width_ratios": [1.25, 1]})
    ax = axes[0]
    for m in marks:
        nt = cpm[:, ORDER.index(f"shNT-{m}")]
        bc = cpm[:, ORDER.index(f"shBCAT1-{m}")]
        ok = keep_hi[m]
        lr = np.log2((bc[ok] + 1e-3) / (nt[ok] + 1e-3))
        ax.hist(lr, bins=160, range=(-3, 3), histtype="step", lw=1.7,
                color=MARK_COLOR[m], density=True,
                label=f"{m} (×{SPIKEIN[m]:.2f}, median {np.median(lr):+.2f})")
        ax.axvline(np.log2(SPIKEIN[m]), color=MARK_COLOR[m], ls=":", lw=1.2)
    ax.axvline(0, color="#999999", lw=0.9)
    ax.set_xlabel("log2( shBCAT1 / shNT ), spike-in scaled — 10 kb bins with signal")
    ax.set_ylabel("density")
    ax.set_title("Occupancy change after spike-in scaling\n"
                 "dotted line = log2(spike-in factor), the expected global shift",
                 fontsize=9.5)
    ax.legend(fontsize=7.4, frameon=False)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    ax = axes[1]
    x = np.arange(len(marks))
    med_unscaled, med_scaled = [], []
    for m in marks:
        nt = cpm[:, ORDER.index(f"shNT-{m}")]
        bc = cpm[:, ORDER.index(f"shBCAT1-{m}")]
        ok = keep_hi[m]
        med_unscaled.append(np.median(np.log2((bc[ok] / SPIKEIN[m] + 1e-3)
                                              / (nt[ok] + 1e-3))))
        med_scaled.append(np.median(np.log2((bc[ok] + 1e-3) / (nt[ok] + 1e-3))))
    ax.bar(x - 0.2, med_unscaled, 0.4, color="#bbbbbb", label="CPM only")
    ax.bar(x + 0.2, med_scaled, 0.4,
           color=[MARK_COLOR[m] for m in marks], label="+ spike-in")
    for i, m in enumerate(marks):
        ax.plot([i + 0.0, i + 0.4], [np.log2(SPIKEIN[m])] * 2, color="k", ls=":", lw=1.2)
    ax.axhline(0, color="#999999", lw=0.9)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{m}\n×{SPIKEIN[m]:.2f}" + ("*" if m in ASSUMED else "")
                        for m in marks], fontsize=8.5)
    ax.set_ylabel("median log2( shBCAT1 / shNT )")
    ax.set_title("Genome-wide shift restored by the spike-in\n"
                 "dotted = log2(factor); * = factor assumed, not supplied",
                 fontsize=9.5)
    ax.legend(fontsize=7.6, frameon=False)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    fig.savefig(path.replace(".png", ".pdf"))
    plt.close(fig)
    print("[plot]", path)
    return dict(zip(marks, zip(med_unscaled, med_scaled)))


def main():
    plots = os.path.join(OUT, "plots")
    data = os.path.join(OUT, "data")
    chroms, starts, C, depths = load_grid_and_counts()
    keep = blacklist_mask(chroms, starts)
    Ck = C[keep]

    cpm_raw = Ck / Ck.sum(axis=0, keepdims=True) * 1e6
    scale = np.array([SPIKEIN[mark_of(p)] if cond_of(p) == "shBCAT1" else 1.0
                      for p in ORDER])
    cpm = cpm_raw * scale
    print("[scale] applied:",
          {p: s for p, s in zip(ORDER, scale) if s != 1.0})

    X_raw = np.log2(cpm_raw + 1)
    X = np.log2(cpm + 1)
    S_raw = sp_matrix(X_raw)
    S = sp_matrix(X)

    # Assert the invariance rather than assume it: rank correlation cannot see a
    # positive rescaling of a sample, so these two matrices must agree.
    dmax = float(np.abs(S - S_raw).max())
    print(f"[check] max |Spearman(scaled) - Spearman(unscaled)| = {dmax:.2e} "
          f"-> scaling is invisible to the correlation matrix, as expected")

    plot_corr(S, ORDER, os.path.join(plots, "corr_heatmap_pooled_spikein.png"),
              "Replicate-pooled profiles — genome-wide Spearman, 10 kb bins",
              "spike-in scaled (ENL ×0.80, H3K9pr ×0.62); rank correlation is "
              "invariant to this scaling — identical to CPM-only")

    # high-signal bins per mark, so the ratio panel is not dominated by empty bins
    keep_hi = {}
    for m in ("H3K4me3", "H3K27ac", "H3K9pr", "ENL"):
        nt = cpm_raw[:, ORDER.index(f"shNT-{m}")]
        bc = cpm_raw[:, ORDER.index(f"shBCAT1-{m}")]
        keep_hi[m] = ((nt + bc) / 2) > np.percentile((nt + bc) / 2, 75)
        print(f"[bins] {m:8s} {keep_hi[m].sum():,} bins above the 75th pct of mean signal")
    meds = plot_occupancy(cpm, keep_hi, os.path.join(plots, "occupancy_shift_spikein.png"))

    np.savez_compressed(os.path.join(data, "pooled_spikein_10kb.npz"),
                        M=C, labels=np.array(ORDER), chroms=chroms, starts=starts,
                        blacklist_keep=keep, spikein=scale)
    np.savetxt(os.path.join(data, "spearman_matrix_pooled_spikein.csv"), S,
               delimiter=",", header=",".join(ORDER), comments="")
    with open(os.path.join(data, "pooled_profiles.csv"), "w") as fh:
        fh.write("profile,mark,condition,replicates_pooled,fragments,spikein_factor,"
                 "median_log2ratio_cpm_only,median_log2ratio_spikein\n")
        for p in ORDER:
            m, c = mark_of(p), cond_of(p)
            mu, ms = meds.get(m, ("", ""))
            fh.write(f"{p},{m},{c},{'+'.join(POOLS[m][c])},{depths[p]},"
                     f"{SPIKEIN[m] if c == 'shBCAT1' else 1.0},"
                     f"{mu if c == 'shBCAT1' else ''},"
                     f"{ms if c == 'shBCAT1' else ''}\n")
    print("[ok] wrote", data)


if __name__ == "__main__":
    main()
