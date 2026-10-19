#!/usr/bin/env python3
"""Figure 3F -- differential expression in Bcat1-KO vs WT mouse brain at P6 and P90.

Input:  BCAT1_KO_analysis/out/deseq2_p6_noWTr3_KOvsWT.csv   (P6, WT replicate 3 excluded)
        BCAT1_KO_analysis/out/deseq2_p90_KOvsWT.csv         (P90)
        BCAT1_KO_analysis/out/volcano_custom_genesets.json  (RGC / OPC_P / OXPHOS)
Output: figures/volcano_P6_P90_customsets.{png,svg,pdf}
        figures/volcano_P6_P90_customsets_zoom12.{png,svg,pdf}   (y clipped at 12)

Back-to-back volcanoes. Points are coloured by membership of three curated sets:
OXPHOS (red), proliferative OPC (purple), pooled RGC (teal); other genes with
padj < 0.05 are dark grey and non-significant genes light grey. Bcat1 itself sits far
above the plotted y range and is annotated with a leader line.

The P6 panel uses the outlier-corrected DESeq2 result (WT replicate 3 dropped); this is
the version whose statistics match the published figure. No P90 outlier-corrected file
exists, so P90 uses the standard result.

Run from anywhere:  python BCAT1_KO_analysis/plot_volcano_P6_P90.py
"""
import os, json
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

SETCOL = {"OXPHOS": "#d1362f", "OPC_P": "#8e44ad", "RGC": "#1abc9c"}
SETLAB = {"OXPHOS": "OXPHOS", "OPC_P": "OPC-P (prolif.)", "RGC": "RGC (pooled)"}
LABEL = {"P6": ["Smc4", "Zbtb20", "Trim59", "Top2a", "Nfib", "Nusap1",
                "Atp6v0a1", "Atp6v0c", "Uqcrc1"], "P90": []}
PADJ = 0.05

sets = json.load(open(os.path.join(HERE, "out/volcano_custom_genesets.json")))
member = {}
for k, genes in sets.items():
    for g in genes:
        member.setdefault(g, k)          # first set wins if a gene is in two


def load(fn):
    d = pd.read_csv(os.path.join(HERE, "out", fn)).dropna(subset=["padj"])
    d["nlp"] = -np.log10(d["padj"].clip(lower=1e-300))
    d["set"] = d.gene.map(member)
    return d


def draw(axes, ymax):
    for ax, (tag, d) in zip(axes, DATA):
        ns = d[d.padj >= PADJ]
        ax.scatter(ns.log2FoldChange, ns.nlp, s=4, color="#dcdcdc",
                   edgecolor="none", zorder=1)
        other = d[(d.padj < PADJ) & d.set.isna()]
        ax.scatter(other.log2FoldChange, other.nlp, s=5, color="#8a8a8a",
                   edgecolor="none", zorder=2)
        for k in ("RGC", "OPC_P", "OXPHOS"):
            s = d[d.set == k]
            ax.scatter(s.log2FoldChange, s.nlp, s=13, color=SETCOL[k],
                       edgecolor="none", zorder=3)
        for g in LABEL[tag]:
            r = d[d.gene == g]
            if len(r):
                r = r.iloc[0]
                ax.annotate(g, (r.log2FoldChange, r.nlp), textcoords="offset points",
                            xytext=(5, 4), fontsize=7, style="italic",
                            color=SETCOL.get(r["set"], "black"), fontweight="bold")
        b = d[d.gene == "Bcat1"]
        if len(b):
            b = b.iloc[0]
            ax.annotate("Bcat1", (b.log2FoldChange, ymax * 0.985),
                        xytext=(b.log2FoldChange + 0.55, ymax * 0.80),
                        arrowprops=dict(arrowstyle="-", lw=0.8, color="black"),
                        fontsize=7.5, style="italic", ha="left")
        ax.axhline(-np.log10(PADJ), ls="--", color="0.55", lw=0.8)
        ax.axvline(0, color="0.55", lw=0.8)
        ax.set_xlim(-2.05, 2.05); ax.set_ylim(0, ymax)
        ax.set_title(f"{tag}  (WT vs KO)", fontsize=9.5)
        ax.set_xlabel("log2 fold change (KO / WT)")
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(width=0.9, length=3)
    axes[0].set_ylabel("-log$_{10}$ adjusted p-value")


DATA = [("P6", load("deseq2_p6_noWTr3_KOvsWT.csv")),
        ("P90", load("deseq2_p90_KOvsWT.csv"))]

for stem, ymax in (("volcano_P6_P90_customsets", 12.5),
                   ("volcano_P6_P90_customsets_zoom12", 12.0)):
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.4), sharey=True,
                             gridspec_kw=dict(wspace=0.10))
    draw(axes, ymax)
    h = ([plt.Line2D([], [], ls="", marker="o", ms=4.5, color=SETCOL[k])
          for k in ("OXPHOS", "OPC_P", "RGC")] +
         [plt.Line2D([], [], ls="", marker="o", ms=4, color="#8a8a8a"),
          plt.Line2D([], [], ls="", marker="o", ms=3.5, color="#dcdcdc")])
    lab = [SETLAB[k] for k in ("OXPHOS", "OPC_P", "RGC")] + \
          [f"other, padj<{PADJ}", "ns"]
    axes[1].legend(h, lab, frameon=False, fontsize=7.5, loc="upper right",
                   handletextpad=0.2, labelspacing=0.3)
    for e, kw in (("png", dict(dpi=400)), ("svg", {}), ("pdf", {})):
        fig.savefig(os.path.join(HERE, "figures", f"{stem}.{e}"),
                    bbox_inches="tight", **kw)
    plt.close(fig); print("wrote figures/%s.{png,svg,pdf}" % stem)

for tag, d in DATA:
    print(f"\n{tag}: {len(d)} genes, {int((d.padj<PADJ).sum())} at padj<{PADJ}; "
          f"set members plotted: " +
          ", ".join(f"{k}={int((d.set==k).sum())}" for k in SETCOL))
