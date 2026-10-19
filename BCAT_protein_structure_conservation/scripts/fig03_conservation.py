#!/usr/bin/env python3
"""Figure 3 - vertebrate conservation (UCSC hg38 100-way phyloP / phastCons)
across the complete coding sequence of BCAT1 and BCAT2."""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
from scipy.stats import mannwhitneyu
from bcat_common import *

aln, sf, cons = load_all()
NT = {"BCAT1": 25, "BCAT2": 36}
COL = {"BCAT1": C1, "BCAT2": C2}


def codon_tracks(g):
    c = cons[g]
    n = len(c["protein"])
    pp = np.array(c["phyloP"], dtype=float)
    pc = np.array(c["phastCons"], dtype=float)
    ppc = pp[:n * 3].reshape(n, 3)
    pcc = pc[:n * 3].reshape(n, 3)
    # exon boundary positions in CDS coordinates (transcript order)
    sizes = [b - a for a, b in c["blocks"]]
    if c["strand"] == "-":
        sizes = sizes[::-1]
    bnds = np.cumsum(sizes)[:-1] / 3.0
    return n, ppc, pcc, bnds, sizes


def smooth(v, w=9):
    return np.convolve(v, np.ones(w) / w, mode="same")


fig = plt.figure(figsize=(7.2, 8.4))
gs = fig.add_gridspec(6, 2, height_ratios=[1.0, 0.8, 1.0, 0.8, 1.3, 1.25],
                      hspace=0.78, wspace=0.30,
                      left=0.085, right=0.975, top=0.945, bottom=0.065)

summary = {}
for r, g in enumerate(["BCAT1", "BCAT2"]):
    n, ppc, pcc, bnds, sizes = codon_tracks(g)
    x = np.arange(1, n + 1)
    mp = ppc.mean(axis=1)
    mc = pcc.mean(axis=1)
    col = COL[g]
    nt = NT[g]

    # ---- phyloP
    ax = fig.add_subplot(gs[r * 2, :])
    ax.set_title("ac"[r], loc="left", fontweight="bold")
    ax.axhline(0, color="#999", lw=0.5)
    ax.axvspan(1, nt, color=col, alpha=0.13, lw=0, zorder=0)
    pos = np.clip(mp, 0, None); neg = np.clip(mp, None, 0)
    ax.bar(x, pos, width=1.0, color="#6f6f7a", lw=0, zorder=2)
    ax.bar(x, neg, width=1.0, color="#c0392b", lw=0, zorder=2)
    ax.plot(x, smooth(mp), color=col, lw=1.1, zorder=4)
    for b in bnds:
        ax.axvline(b, color="#b9b9c2", lw=0.4, ls=(0, (2, 2)), zorder=1)
    edges = np.concatenate([[0], np.cumsum(sizes) / 3.0])
    for i in range(len(sizes)):
        ax.add_patch(Rectangle((edges[i], 8.5), edges[i + 1] - edges[i], 0.75,
                               fc=col if i % 2 == 0 else "#c9ccd1", ec="w", lw=0.4,
                               clip_on=False, zorder=5))
    ax.text(-0.012, 8.9, "CDS\nexons", transform=ax.get_yaxis_transform(),
            ha="right", va="center", fontsize=6, color=GREY)
    ax.set_ylim(-2.2, 8.2); ax.set_xlim(0, 400)
    ax.set_ylabel("phyloP\n(100-way)")
    ax.spines["bottom"].set_bounds(0, n)
    ax.text(0.998, 1.18, f"{g}  ({cons[g]['transcript']}, {cons[g]['chrom']}, "
            f"{cons[g]['strand']} strand, {n*3} coding bases)",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=6.5, color=col)
    ax.text(nt + 3, 7.7, f"N-terminal 1–{nt}", fontsize=6, color=col, va="top")
    ax.tick_params(labelbottom=False)
    if r == 0:
        hh = [Patch(fc="#6f6f7a"), Patch(fc="#c0392b")]
        ax.legend(hh, ["conserved (phyloP > 0)", "accelerated (phyloP < 0)"],
                  ncol=2, loc="lower right", bbox_to_anchor=(1.0, -0.04),
                  handlelength=1.0, handleheight=0.8, fontsize=6)

    # ---- phastCons + exon map
    ax2 = fig.add_subplot(gs[r * 2 + 1, :])
    ax2.fill_between(x, 0, mc, color="#3c3c44", lw=0, alpha=0.85, step="mid")
    ax2.axvspan(1, nt, color=col, alpha=0.13, lw=0, zorder=0)
    ax2.set_ylim(0, 1.05); ax2.set_xlim(0, 400)
    ax2.set_yticks([0, 0.5, 1])
    ax2.set_ylabel("phastCons")
    ax2.set_xlabel(f"{g} residue  (codon position along the complete CDS)")
    ax2.spines["bottom"].set_bounds(0, n)
    ax2.set_title("bd"[r], loc="left", fontweight="bold")

    # ---- stats
    core = mp[nt:]
    ntv = mp[:nt]
    u, p = mannwhitneyu(ntv, core, alternative="two-sided")
    summary[g] = dict(n=n, nt=nt, mean_nt=ntv.mean(), mean_core=core.mean(),
                      med_nt=np.median(ntv), med_core=np.median(core), p=p,
                      pc_nt=mc[:nt].mean(), pc_core=mc[nt:].mean(),
                      by_codon=ppc.mean(axis=0), mp=mp, mc=mc)
    print(f"{g}: phyloP N-term 1-{nt} mean {ntv.mean():.2f} vs core {core.mean():.2f} "
          f"(Mann-Whitney p = {p:.2e}); phastCons {mc[:nt].mean():.2f} vs {mc[nt:].mean():.2f}")

# ------------------------------------------------ e: N-term vs core distribution
ax = fig.add_subplot(gs[4, 0]); ax.set_title("e", loc="left", fontweight="bold")
data, labels, colours = [], [], []
for g in ["BCAT1", "BCAT2"]:
    s = summary[g]
    data += [s["mp"][:s["nt"]], s["mp"][s["nt"]:]]
    labels += [f"{g}\nN-term", f"{g}\ncore"]
    colours += [COL[g], "#c9ccd1"]
bp = ax.boxplot(data, widths=0.6, showfliers=False, patch_artist=True,
                medianprops=dict(color="#16161a", lw=1.0),
                whiskerprops=dict(lw=0.7), capprops=dict(lw=0.7),
                boxprops=dict(lw=0.7))
for patch, c in zip(bp["boxes"], colours):
    patch.set_facecolor(c); patch.set_alpha(0.85)
rng = np.random.default_rng(0)
for i, d in enumerate(data):
    ax.plot(i + 1 + rng.normal(0, 0.07, len(d)), d, ".", ms=1.1, color="#3c3c44",
            alpha=0.35, zorder=0)
ax.set_xticklabels(labels, fontsize=6)
ax.axhline(0, color="#999", lw=0.5)
ax.set_ylabel("per-codon mean phyloP")
ax.set_ylim(-2.4, 8.5)
for i, g in enumerate(["BCAT1", "BCAT2"]):
    p = summary[g]["p"]
    y = 7.6
    ax.plot([1 + 2 * i, 2 + 2 * i], [y, y], color="#16161a", lw=0.6)
    ax.text(1.5 + 2 * i, y + 0.15, f"p = {p:.1e}", ha="center", va="bottom", fontsize=5.8)

# ------------------------------------------------ f: conservation by codon position
ax = fig.add_subplot(gs[4, 1]); ax.set_title("f", loc="left", fontweight="bold")
w = 0.36
for i, g in enumerate(["BCAT1", "BCAT2"]):
    v = summary[g]["by_codon"]
    ax.bar(np.arange(3) + (i - 0.5) * w, v, width=w, color=COL[g], label=g, lw=0)
ax.set_xticks(range(3)); ax.set_xticklabels(["1st", "2nd", "3rd\n(wobble)"])
ax.set_xlabel("codon position")
ax.set_ylabel("mean phyloP")
ax.legend(loc="upper right", handlelength=1.0, handleheight=0.8)

# ------------------------------------------------ g: aligned comparison
ax = fig.add_subplot(gs[5, :]); ax.set_title("g", loc="left", fontweight="bold")
a1, a2 = aln["a1"], aln["a2"]
m1, m2 = col_maps(a1, a2)
L = len(a1)
for g, m in [("BCAT1", m1), ("BCAT2", m2)]:
    v = np.full(L, np.nan)
    mp = summary[g]["mp"]
    for c in range(L):
        if m[c]:
            v[c] = mp[m[c] - 1]
    ok = ~np.isnan(v)
    sm = np.full(L, np.nan)
    idx = np.where(ok)[0]
    sm[idx] = np.convolve(v[idx], np.ones(11) / 11, mode="same")
    ax.plot(np.arange(L), sm, color=COL[g], lw=1.0, label=g)
ntcols = [c for c in range(L) if (m1[c] and m1[c] <= NT["BCAT1"]) or
          (m2[c] and m2[c] <= NT["BCAT2"])]
ax.axvspan(0, max(ntcols) + 1, color="#f0c9a8", alpha=0.45, lw=0, zorder=0)
ax.text(max(ntcols) / 2, 7.0, "divergent N-terminal region", ha="center", fontsize=6,
        color="#a85a12")
ax.axhline(0, color="#999", lw=0.5)
ax.set_xlim(-4, L + 4); ax.set_ylim(-1.0, 7.6)
ax.set_xlabel("alignment column (BCAT1 reference)")
ax.set_ylabel("phyloP\n(11-codon window)")
ax.legend(loc="upper right", ncol=2, handlelength=1.2)

fig.suptitle("Vertebrate conservation across the BCAT1 and BCAT2 coding sequences "
             "(UCSC hg38, 100-way phyloP / phastCons)", y=0.995, fontsize=8.5)
save(fig, "Fig3_BCAT_conservation_UCSC100way")

# text summary for the methods note
import json
json.dump({g: {k: (float(v) if np.isscalar(v) else np.asarray(v).tolist())
               for k, v in s.items() if k not in ("mp", "mc")}
           for g, s in summary.items()},
          open("/root/bc/out/conservation_summary.json", "w"), indent=1)
