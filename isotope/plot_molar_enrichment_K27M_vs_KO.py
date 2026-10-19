#!/usr/bin/env python3
"""Figure 4C -- 13C molar enrichment of acyl-CoA species, K27M vs isogenic KO.

Reads the tidy `K27M-KO molar enrichments` sheet (per-replicate molar enrichment %,
one row per replicate, one column per acyl-CoA) and draws a 3 x 2 grid of horizontal
bar panels: rows are the three tracers (isoleucine, valine, leucine), columns are the
two cell lines.

Statistics (recovered and verified against the published q values, see
PROVENANCE_AUDIT.md): Welch unequal-variance two-tailed t-test per acyl-CoA, then
Benjamini-Hochberg FDR across the metabolites quantified within each cell line x tracer
block. q is annotated where q < 0.05.

The upstream derivation (natural-abundance correction and the molar-enrichment formula)
is NOT reproduced here -- it is not recorded anywhere in the workbook.

Run from anywhere:  python isotope/plot_molar_enrichment_K27M_vs_KO.py
"""
import os
import numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
from _isotope_common import load, set_style, save, PALETTE, HERE

set_style()
RNG = np.random.default_rng(3)

C = {"K27M": PALETTE["conditions"]["K27M"]["hex"],
     "KO":   PALETTE["conditions"]["KO"]["hex"]}

TRACERS = [("Isoleucine-C13", "Isoleucine"),
           ("Valine-C13",     "Valine"),
           ("Leucine-C13",    "Leucine")]
# plotted bottom-to-top, so list in reverse of the figure's top-to-bottom order
METS = ["Succinyl-CoA", "Propionyl-CoA", "Malonyl-CoA", "Glutaryl-CoA",
        "Butyryl-CoA", "Acetyl-CoA", "3HMG-CoA", "3HB-CoA"]

df = load("K27M-KO molar enrichments")
for m in METS:
    df[m] = pd.to_numeric(df[m], errors="coerce")
LINES = [l for l in ["BT245", "DIPGXIII", "DIPG13"] if l in set(df["Cell line"])][:2]

fig, axes = plt.subplots(len(TRACERS), len(LINES), figsize=(6.4, 7.2),
                         sharex=True, gridspec_kw=dict(wspace=0.30, hspace=0.16))
rows = []

for r, (tracer, tlabel) in enumerate(TRACERS):
    for c, line in enumerate(LINES):
        ax = axes[r, c]
        sub = df[(df["Cell line"] == line) & (df["Isotope"] == tracer)]

        pvals, keep = [], []
        for m in METS:
            k = sub[sub.Genotype == "K27M"][m].dropna().to_numpy(float)
            o = sub[sub.Genotype == "KO"][m].dropna().to_numpy(float)
            if len(k) > 1 and len(o) > 1 and (k.std() + o.std()) > 0:
                pvals.append(stats.ttest_ind(k, o, equal_var=False).pvalue)
                keep.append(m)
        q = dict(zip(keep, multipletests(pvals, method="fdr_bh")[1])) if keep else {}

        for i, m in enumerate(METS):
            for off, g in ((0.19, "K27M"), (-0.19, "KO")):
                v = sub[sub.Genotype == g][m].dropna().to_numpy(float)
                if not len(v):
                    continue
                ax.barh(i + off, v.mean(), height=0.34, color=C[g],
                        edgecolor="none", zorder=2,
                        xerr=v.std(ddof=1) / np.sqrt(len(v)) if len(v) > 1 else None,
                        error_kw=dict(elinewidth=0.9, capsize=2, capthick=0.9))
                ax.scatter(v, i + off + RNG.uniform(-0.07, 0.07, len(v)), s=11,
                           facecolor="white", edgecolor="black", linewidth=0.5,
                           zorder=4)
            if m in q:
                kk = sub[sub.Genotype == "K27M"][m].dropna().to_numpy(float)
                oo = sub[sub.Genotype == "KO"][m].dropna().to_numpy(float)
                rows.append(dict(figure="Figure 4C", cell_line=line, tracer=tracer,
                                 metabolite=m, n_K27M=len(kk), n_KO=len(oo),
                                 mean_K27M=kk.mean(), mean_KO=oo.mean(),
                                 p=pvals[keep.index(m)], q=q[m]))
                if q[m] < 0.05:
                    xm = max(kk.mean(), oo.mean())
                    ax.text(xm * 1.04 + 0.4, i, f"q={q[m]:.2g}", va="center",
                            fontsize=6.4)

        ax.set_yticks(range(len(METS)))
        ax.set_yticklabels(METS if c == 0 else [], fontsize=7.5)
        ax.set_ylim(-0.6, len(METS) - 0.4)
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(width=0.9, length=3, labelsize=7.5)
        if r == 0:
            ax.set_title(line, fontsize=9.5, fontweight="bold", pad=4)
        if c == len(LINES) - 1:
            ax.text(1.02, 0.5, tlabel, transform=ax.transAxes, rotation=270,
                    va="center", ha="left", fontsize=9)
        if r == len(TRACERS) - 1:
            ax.set_xlabel("Molar enrichment (%)", fontsize=8.5)

h = [plt.Rectangle((0, 0), 1, 1, fc=C[g], ec="none") for g in ("K27M", "KO")]
fig.legend(h, ["K27M", "KO"], loc="upper center", bbox_to_anchor=(0.5, 0.95),
           ncol=2, frameon=False, fontsize=8.5, handlelength=1.1,
           columnspacing=1.4)

save(fig, "K27M_vs_KO_molar_enrichment")
st = pd.DataFrame(rows)
st["test"] = ("Welch unequal-variance two-tailed t-test per acyl-CoA; "
              "Benjamini-Hochberg FDR across metabolites within each "
              "cell line x tracer block")
st.to_csv(os.path.join(HERE, "K27M_vs_KO_molar_enrichment_stats.csv"), index=False)
print()
print(st[st.q < 0.05].to_string(index=False))
