#!/usr/bin/env python3
"""Figure 3G -- GSEA running enrichment for the Tabula Muris OPC ageing signature at P6.

Input:  BCAT1_KO_analysis/out/deseq2_p6_noWTr3_KOvsWT.csv   ranked gene list
        BCAT1_KO_analysis/gmt/MSigDB_m8_mouse.gmt           gene set
        BCAT1_KO_analysis/gsea/MSigDB_m8_mouse_p6_noWTr3.csv  stored NES / FDR
Output: figures/gsea_enrichment_OPC_P6.{png,svg,pdf}

Gene set: TABULA_MURIS_SENIS_BRAIN_NON_MYELOID_OLIGODENDROCYTE_PRECURSOR_CELL_AGEING.

The running enrichment curve is recomputed here deterministically from the ranked list
(genes ordered by the DESeq2 Wald statistic, descending = up in KO first; classic
weighted Kolmogorov-Smirnov walk with weight p=1). This reproduces the stored
ES = 0.5358 exactly, which confirms the ranking metric.

NES and FDR are READ FROM the stored gseapy output rather than recomputed, because they
depend on a random permutation null that would not reproduce bit-for-bit. The script
asserts that its recomputed ES matches the stored ES before plotting.

Run from anywhere:  python BCAT1_KO_analysis/plot_gsea_enrichment_OPC_P6.py
"""
import os
import numpy as np, pandas as pd, matplotlib
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


ROOT = repo_root(); HERE = os.path.join(ROOT, "BCAT1_KO_analysis")
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans"],
                     "font.size": 9, "axes.linewidth": 0.9, "svg.fonttype": "none",
                     "pdf.fonttype": 42, "ps.fonttype": 42})

TERM = "TABULA_MURIS_SENIS_BRAIN_NON_MYELOID_OLIGODENDROCYTE_PRECURSOR_CELL_AGEING"
COLOR = "#1a9e77"

# ---- ranked list (DESeq2 Wald stat, descending) ----
d = (pd.read_csv(os.path.join(HERE, "out/deseq2_p6_noWTr3_KOvsWT.csv"))
     .dropna(subset=["stat"]).sort_values("stat", ascending=False))

# ---- gene set ----
gs = set()
for line in open(os.path.join(HERE, "gmt/MSigDB_m8_mouse.gmt")):
    f = line.rstrip("\n").split("\t")
    if f[0] == TERM:
        gs = set(f[2:]); break
if not gs:
    raise SystemExit("gene set not found in gmt")

# ---- weighted KS walk ----
hit = d.gene.isin(gs).to_numpy()
w = np.abs(d["stat"].to_numpy())
Phit = np.cumsum(np.where(hit, w, 0.0)) / w[hit].sum()
Pmiss = np.cumsum(~hit) / (~hit).sum()
run = Phit - Pmiss
peak = int(np.argmax(np.abs(run)))
ES = run[peak]

stored = pd.read_csv(os.path.join(HERE, "gsea/MSigDB_m8_mouse_p6_noWTr3.csv"))
row = stored[stored.Term == TERM].iloc[0]
assert abs(ES - row["ES"]) < 1e-3, \
    f"recomputed ES {ES:.4f} != stored {row['ES']:.4f}; ranking metric has changed"
NES, FDR = row["NES"], row["FDR q-val"]
print(f"recomputed ES = {ES:+.4f}  (stored {row['ES']:+.4f})  "
      f"NES = {NES:.2f}  FDR = {FDR:.3g}  hits = {int(hit.sum())}/{len(gs)}")

fig, (a1, a2) = plt.subplots(2, 1, figsize=(3.9, 3.2), sharex=True,
                             gridspec_kw=dict(height_ratios=[4, 1], hspace=0.06))
x = np.arange(len(run))
a1.plot(x, run, color=COLOR, lw=1.4, zorder=3)
a1.fill_between(x[:peak + 1], run[:peak + 1], 0, color=COLOR, alpha=0.18,
                linewidth=0, zorder=2)
a1.axhline(0, color="black", lw=0.9)
a1.axvline(peak, color="0.35", ls=":", lw=0.9)
a1.text(peak + len(x) * 0.02, ES * 0.55,
        f"NES {NES:.1f}, FDR{'<0.001' if FDR < 1e-3 else f'={FDR:.3f}'}",
        fontsize=8, va="center")
a1.set_ylabel("ES")
a1.set_title("Tabula Muris OPC ageing signature (P6)", fontsize=9, pad=4)
a1.spines[["top", "right"]].set_visible(False)
a1.tick_params(width=0.9, length=3)

a2.vlines(x[hit], 0, 1, color="0.25", lw=0.35)
a2.set_ylim(0, 1); a2.set_yticks([])
a2.set_xlim(0, len(x))
a2.set_xlabel("Gene rank (up→down in KO)")
a2.spines[["top", "right", "left"]].set_visible(False)
a2.tick_params(width=0.9, length=3)

for e, kw in (("png", dict(dpi=400)), ("svg", {}), ("pdf", {})):
    fig.savefig(os.path.join(HERE, "figures", f"gsea_enrichment_OPC_P6.{e}"),
                bbox_inches="tight", **kw)
print("wrote figures/gsea_enrichment_OPC_P6.{png,svg,pdf}")
