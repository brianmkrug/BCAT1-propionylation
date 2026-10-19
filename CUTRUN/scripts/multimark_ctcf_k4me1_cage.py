#!/usr/bin/env python
"""
CORRECTION (2026-10-05): CTCF and H3K4me1 are CUT&RUN from shNT DIPG-XIII cells
(Diagenode C15410210; CST 5326). The xChIP/ChIPRx/DMSO/K27M-EV strings in the original
fastq names are legacy labels, not sample metadata; text below that calls them ChIP or a
different treatment arm is superseded.

Correlation of the CUT&RUN panel with CTCF, H3K4me1 and CAGE, on two grids, plus
a quantitative test of the proposed two-cluster structure.

TWO GRIDS, ON PURPOSE
  genome-wide 10 kb : the usual view, but CAGE covers 0.088% of the genome at
                      1 bp, so most bins are exactly zero and Spearman is
                      tie-dominated and attenuated FOR EVERY MARK EQUALLY. Ranks
                      between marks stay meaningful; absolute values do not.
  promoter windows  : core.bed, TSS +/- 1 kb, where CAGE actually lives. This is
                      the informative grid for "does H3K4me3 track initiation
                      more closely than H3K27ac does" -- and that particular
                      comparison is ASSAY-CONTROLLED, since H3K4me3 and H3K27ac
                      are both CUT&RUN from the same cells.

THE CONFOUND (see PLAN_ctcf_k4me1_cage.md 0)
[Superseded -- see CORRECTION above.] Originally read as ChIP from DMSO and K27M-EV
samples rather than shNT/shBCAT1. The proposed "other" group is exactly where both ChIPs sit, so a
split is confounded with assay. Two things are reported to address it:
  * H3K27ac is CUT&RUN and already separates from H3K4me3/ENL/H3K9pr, so the
    split exists within CUT&RUN alone;
  * the pre-specified discriminator -- if ASSAY drives clustering the two ChIPs
    pair with each other; if BIOLOGY drives it, H3K4me1 joins H3K27ac and CTCF
    sits apart from both. This is printed as a verdict.
"""
from __future__ import annotations
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root

import gzip
import itertools
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform
from scipy.stats import rankdata

ROOT = os.path.expanduser(REPO_ROOT + "/CUTRUN")
QC = os.path.join(ROOT, "analysis_QC_normalization", "multimark_2026-07")
GRID = os.path.join(ROOT, "analysis_QC_normalization", "data", "counts_10kb_matrix.npz")
BLACKLIST = os.path.join(ROOT, "hg38-blacklist.v2.bed.gz")
CORE = os.path.join(ROOT, "H3K4me3", "reference", "core.bed")
CAGE_S = "S.K27M"
BIN = 10_000
NBOOT = 100
SEED = 7

# pooled shNT CUT&RUN profiles + the three additions
POOLS = {
    "H3K4me3": ["shNT-H3K4me3-r1", "shNT-H3K4me3-r2"],
    "H3K27ac": ["shNT-H3K27ac-r1", "shNT-H3K27ac-r2"],
    "H3K9pr":  ["shNT-H3K9pr-r1", "shNT-H3K9pr-r2"],
    "ENL":     ["shNT-ENL-r1", "shNT-ENL-r2"],
    "IgG":     ["shNT-IgG-r1"],
    "CTCF":    ["CTCF"],
    "H3K4me1": ["H3K4me1"],
}
ORDER = ["H3K4me3", "ENL", "H3K9pr", "H3K27ac", "H3K4me1", "CTCF", "IgG", "CAGE"]
GROUP_A = ["H3K4me3", "ENL", "H3K9pr"]      # proposed promoter/elongation cluster
GROUP_B = ["H3K27ac", "H3K4me1", "CTCF"]    # proposed enhancer/insulator cluster
ASSAY = {"H3K4me3": "CUT&RUN", "ENL": "CUT&RUN", "H3K9pr": "CUT&RUN",
         "H3K27ac": "CUT&RUN", "IgG": "CUT&RUN",
         "H3K4me1": "CUT&RUN", "CTCF": "CUT&RUN", "CAGE": "CAGE"}
COLOR = {"H3K4me3": "#2166ac", "ENL": "#1e8449", "H3K9pr": "#6c3483",
         "H3K27ac": "#e8862a", "H3K4me1": "#d94a8c", "CTCF": "#7d6608",
         "IgG": "#808080", "CAGE": "#c0392b"}


def sp_matrix(X):
    R = np.apply_along_axis(rankdata, 0, X)
    R = R - R.mean(axis=0, keepdims=True)
    R /= np.linalg.norm(R, axis=0, keepdims=True)
    return R.T @ R


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
            if f[0] not in by:
                continue
            for s in range((int(f[1]) // BIN) * BIN,
                           ((int(f[2]) - 1) // BIN) * BIN + BIN, BIN):
                i = by[f[0]].get(s)
                if i is not None:
                    keep[i] = False
    return keep


def load_10kb():
    d = np.load(GRID, allow_pickle=True)
    chroms = np.asarray(d["chroms"], dtype="U6")
    starts = np.asarray(d["starts"], dtype=np.int64)
    idx = {(c, s): i for i, (c, s) in enumerate(zip(chroms, starts))}
    M = np.zeros((len(chroms), len(ORDER)))
    for j, prof in enumerate(ORDER):
        if prof == "CAGE":
            p = os.path.join(QC, "cage", f"{CAGE_S}.cage_10kb.tsv.gz")
            with gzip.open(p, "rt") as fh:
                for line in fh:
                    c, s, v = line.split()
                    k = idx.get((c, int(s)))
                    if k is not None:
                        M[k, j] += float(v)
            continue
        for lab in POOLS[prof]:
            with gzip.open(os.path.join(QC, "counts", f"{lab}.10kb.tsv.gz"), "rt") as fh:
                for line in fh:
                    c, s, n = line.split()
                    k = idx.get((c, int(s)))
                    if k is not None:
                        M[k, j] += float(n)
    return chroms, starts, M


def load_promoter():
    n = sum(1 for _ in open(CORE))
    M = np.zeros((n, len(ORDER)))
    for j, prof in enumerate(ORDER):
        if prof == "CAGE":
            p = os.path.join(QC, "cage", f"{CAGE_S}.cage_promoters.tsv.gz")
            t = pd.read_csv(p, sep="\t", header=None, names=["row", "v"])
            M[t.row.to_numpy() - 1, j] = t.v.to_numpy()
            continue
        for lab in POOLS[prof]:
            t = pd.read_csv(os.path.join(QC, "promcat", f"{lab}.promoter.tsv.gz"), sep="\t")
            M[t.promoter_row.to_numpy() - 1, j] += t.fragments.to_numpy()
    return M


def cpm_log(M):
    tot = M.sum(axis=0, keepdims=True)
    tot[tot == 0] = 1
    return np.log2(M / tot * 1e6 + 1)


# ------------------------------------------------------------- cluster test --
def cluster_stats(S, labels, rng, X=None):
    ia = [labels.index(m) for m in GROUP_A]
    ib = [labels.index(m) for m in GROUP_B]

    def stat(Sx):
        w = ([Sx[i, j] for i, j in itertools.combinations(ia, 2)] +
             [Sx[i, j] for i, j in itertools.combinations(ib, 2)])
        b = [Sx[i, j] for i in ia for j in ib]
        return float(np.mean(w) - np.mean(b))

    obs = stat(S)

    # exact enumeration of every 3/3 split of the six marks. Only 10 distinct
    # bipartitions exist, so the finest achievable p-value is 0.1 -- reported as
    # a rank rather than dressed up as significance.
    six = GROUP_A + GROUP_B
    parts, seen = [], set()
    for combo in itertools.combinations(six, 3):
        key = frozenset([frozenset(combo), frozenset(set(six) - set(combo))])
        if key in seen:
            continue
        seen.add(key)
        ga = [labels.index(m) for m in combo]
        gb = [labels.index(m) for m in set(six) - set(combo)]
        w = ([S[i, j] for i, j in itertools.combinations(ga, 2)] +
             [S[i, j] for i, j in itertools.combinations(gb, 2)])
        b = [S[i, j] for i in ga for j in gb]
        parts.append((float(np.mean(w) - np.mean(b)), tuple(sorted(combo))))
    parts.sort(reverse=True)
    rank = 1 + [p[1] for p in parts].index(tuple(sorted(GROUP_A)))

    # silhouette on correlation distance
    D = 1 - S
    sil = []
    for m in six:
        i = labels.index(m)
        own, oth = (ia, ib) if m in GROUP_A else (ib, ia)
        a = np.mean([D[i, j] for j in own if j != i])
        b = np.mean([D[i, j] for j in oth])
        sil.append((m, (b - a) / max(a, b)))

    # bootstrap over bins -> CI on the separation statistic
    boot = []
    if X is not None:
        n = X.shape[0]
        for _ in range(NBOOT):
            idx = rng.integers(0, n, n)
            boot.append(stat(sp_matrix(X[idx])))
        boot = np.array(boot)
    return obs, parts, rank, sil, boot


def discriminator(S, labels):
    """Pre-specified: assay-driven -> the two ChIPs pair; biology-driven ->
    H3K4me1 joins H3K27ac and CTCF sits apart from both."""
    g = lambda a, b: S[labels.index(a), labels.index(b)]
    chip_chip = g("CTCF", "H3K4me1")
    k4me1_k27ac = g("H3K4me1", "H3K27ac")
    ctcf_k27ac = g("CTCF", "H3K27ac")
    if k4me1_k27ac > chip_chip and k4me1_k27ac > ctcf_k27ac:
        verdict = ("BIOLOGY — H3K4me1 sits with H3K27ac (%.3f) above the "
                   "ChIP-ChIP pairing (%.3f); the two ChIPs are split by biology"
                   % (k4me1_k27ac, chip_chip))
    elif chip_chip > k4me1_k27ac:
        verdict = ("ASSAY — the two ChIP samples pair with each other (%.3f) "
                   "above H3K4me1-H3K27ac (%.3f); ChIP data cannot settle the "
                   "grouping" % (chip_chip, k4me1_k27ac))
    else:
        verdict = "AMBIGUOUS"
    return chip_chip, k4me1_k27ac, ctcf_k27ac, verdict


def plot_corr(S, labels, path, title, subtitle):
    D = squareform(np.clip(1 - S, 0, None), checks=False)
    Z = linkage(D, method="average")
    order = dendrogram(Z, no_plot=True)["leaves"]
    Sx = S[np.ix_(order, order)]
    lx = [labels[i] for i in order]
    fig = plt.figure(figsize=(8.6, 8.0))
    gs = fig.add_gridspec(2, 2, width_ratios=[40, 1], height_ratios=[1, 5.0],
                          wspace=0.03, hspace=0.02, left=0.22, right=0.91,
                          top=0.89, bottom=0.20)
    axd = fig.add_subplot(gs[0, 0])
    dendrogram(Z, ax=axd, color_threshold=0, above_threshold_color="#666666", no_labels=True)
    axd.set_axis_off()
    ax = fig.add_subplot(gs[1, 0])
    im = ax.imshow(Sx, cmap="RdYlBu_r", vmin=np.floor(S.min() * 20) / 20, vmax=1.0)
    ticks = [f"{l}\n({ASSAY[l]})" for l in lx]
    ax.set_xticks(range(len(lx))); ax.set_xticklabels(ticks, rotation=90, fontsize=7.6)
    ax.set_yticks(range(len(lx))); ax.set_yticklabels(ticks, fontsize=7.6)
    for i in range(len(lx)):
        for j in range(len(lx)):
            ax.text(j, i, f"{Sx[i, j]:.2f}", ha="center", va="center", fontsize=6.6,
                    color="white" if Sx[i, j] > 0.80 or Sx[i, j] < 0.35 else "#222")
    for t, l in zip(ax.get_xticklabels(), lx): t.set_color(COLOR[l])
    for t, l in zip(ax.get_yticklabels(), lx): t.set_color(COLOR[l])
    cax = fig.add_subplot(gs[1, 1])
    cb = fig.colorbar(im, cax=cax); cb.set_label("Spearman r", fontsize=8)
    cb.ax.tick_params(labelsize=7)
    fig.suptitle(title, fontsize=11, y=0.965)
    fig.text(0.5, 0.925, subtitle, ha="center", fontsize=8, color="#555")
    fig.savefig(path, dpi=200); fig.savefig(path.replace(".png", ".pdf"))
    plt.close(fig)
    print("[plot]", path)


def report(name, S, labels, X, rng, fh):
    obs, parts, rank, sil, boot = cluster_stats(S, labels, rng, X)
    cc, kk, ck, verdict = discriminator(S, labels)
    lo, hi = (np.percentile(boot, [2.5, 97.5]) if len(boot) else (np.nan, np.nan))
    fh.write(f"\n=== {name} ===\n")
    fh.write(f"separation (mean within-group r - mean between-group r) = {obs:+.4f}\n")
    fh.write(f"  bootstrap 95% CI over bins: [{lo:+.4f}, {hi:+.4f}]  "
             f"({'excludes 0' if lo > 0 else 'includes 0'})\n")
    fh.write(f"  rank of the proposed split among all 10 possible 3/3 splits: {rank}/10\n")
    fh.write("  best splits by separation:\n")
    for v, c in parts[:3]:
        fh.write(f"    {v:+.4f}  {{{', '.join(c)}}}\n")
    fh.write("  silhouette width per mark:\n")
    for m, s in sil:
        fh.write(f"    {m:9s} {s:+.3f}\n")
    fh.write(f"  discriminator: CTCF-H3K4me1 {cc:.3f} | H3K4me1-H3K27ac {kk:.3f} | "
             f"CTCF-H3K27ac {ck:.3f}\n    -> {verdict}\n")
    if "CAGE" in labels:
        j = labels.index("CAGE")
        fh.write("  CAGE correlation with each mark (assay-controlled pairs in bold):\n")
        for m in ORDER:
            if m == "CAGE":
                continue
            star = " *" if m in ("H3K4me3", "H3K27ac") else ""
            fh.write(f"    CAGE-{m:9s} {S[j, labels.index(m)]:+.3f}{star}\n")
        fh.write("    * H3K4me3 and H3K27ac are both CUT&RUN from the same cells,\n"
                 "      so their difference here is not an assay artefact\n")
    return obs, lo, hi, rank


def main():
    rng = np.random.default_rng(SEED)
    plots, data = os.path.join(QC, "plots"), os.path.join(QC, "data")
    chroms, starts, M10 = load_10kb()
    keep = blacklist_mask(chroms, starts)
    M10 = M10[keep]
    print(f"[10kb] {M10.shape[0]:,} bins after blacklist")
    nzc = (M10[:, ORDER.index("CAGE")] > 0).sum()
    print(f"[10kb] CAGE non-zero in {nzc:,} bins ({100*nzc/M10.shape[0]:.1f}%) "
          f"-> genome-wide Spearman is tie-attenuated")

    Mp = load_promoter()
    nzp = (Mp[:, ORDER.index("CAGE")] > 0).sum()
    print(f"[prom] {Mp.shape[0]:,} promoters; CAGE non-zero in {nzp:,} "
          f"({100*nzp/Mp.shape[0]:.1f}%)")

    X10, Xp = cpm_log(M10), cpm_log(Mp)
    S10, Sp = sp_matrix(X10), sp_matrix(Xp)

    plot_corr(S10, ORDER, os.path.join(plots, "corr_heatmap_withCTCF_K4me1_CAGE_10kb.png"),
              "Genome-wide Spearman, 10 kb bins — CUT&RUN panel + CTCF, H3K4me1, CAGE",
              f"CAGE non-zero in only {100*nzc/M10.shape[0]:.1f}% of bins — its values here are tie-attenuated")
    plot_corr(Sp, ORDER, os.path.join(plots, "corr_heatmap_withCTCF_K4me1_CAGE_promoters.png"),
              "Spearman over promoter windows (TSS ±1 kb) — the informative grid for CAGE",
              "H3K4me3 vs H3K27ac against CAGE is assay-controlled: both CUT&RUN, same cells")

    with open(os.path.join(data, "cluster_separation_test.txt"), "w") as fh:
        fh.write("Test of the proposed grouping\n")
        fh.write(f"  A = {GROUP_A}\n  B = {GROUP_B}\n")
        fh.write("NOTE: CTCF and H3K4me1 are single shNT CUT&RUN libraries (no replicate);\n"
                 "they have a longer fragment-length profile than the rest of the panel.\n")
        report("genome-wide 10 kb", S10, ORDER, X10, rng, fh)
        report("promoter windows", Sp, ORDER, Xp, rng, fh)
    print(open(os.path.join(data, "cluster_separation_test.txt")).read())

    for nm, S in (("10kb", S10), ("promoters", Sp)):
        pd.DataFrame(S, index=ORDER, columns=ORDER).to_csv(
            os.path.join(data, f"spearman_withCTCF_K4me1_CAGE_{nm}.csv"))
    print("[ok] wrote", data)


if __name__ == "__main__":
    main()
