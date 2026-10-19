#!/usr/bin/env python3
"""
Final compact version of the adjacent xenograft BCAA volcano plots.
- ~7 cm x 3 cm, 3 pt text.
- red/blue up/down colouring REMOVED: all non-BCAA points neutral grey.
- BCAA-pathway metabolites (Leucine/Isoleucine/Valine sub-pathway): green
  (filled = significant, hollow = n.s.). Top BCAA hits labelled.
- shared symmetric x-range; exports SVG + PNG.
"""
import os, re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42

from matplotlib.lines import Line2D
from scipy import stats
from statsmodels.stats.multitest import multipletests
from adjustText import adjust_text

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(HERE)
XLSX = os.path.join(BASE, "norm_cell_data copy.xlsx")
FIG = os.path.join(HERE, "figures")
META = ["Cell line", "Genotype", "Source"]
LINES = ["HSJ019", "DIPGXIII"]
DROP_LABEL = ("HSJ019", "K27M", 2)
FDR_THR, LFC_THR = 0.05, 1.0
BCAA_SUB = "Leucine, Isoleucine and Valine Metabolism"
N_LABEL = 4
NS_C, BCAA_C = "#bcbcbc", "#2ca02c"
FS = 3  # text size (pt)

CM = 1 / 2.54
plt.rcParams.update({
    "font.size": FS, "axes.titlesize": FS, "axes.labelsize": FS,
    "xtick.labelsize": FS, "ytick.labelsize": FS, "legend.fontsize": FS,
    "axes.linewidth": 0.3, "xtick.major.width": 0.3, "ytick.major.width": 0.3,
    "xtick.major.size": 1.3, "ytick.major.size": 1.3, "lines.linewidth": 0.4,
    "svg.fonttype": "none"})


def norm(s): return re.sub(r"\s+", " ", str(s).strip().lower())


def load():
    df = pd.read_excel(XLSX, sheet_name="xenograft")
    cnt, keep = {}, []
    for _, r in df.iterrows():
        k = (r["Cell line"], r["Genotype"]); cnt[k] = cnt.get(k, 0) + 1
        keep.append(not (r["Cell line"] == DROP_LABEL[0] and
                         r["Genotype"] == DROP_LABEL[1] and cnt[k] == DROP_LABEL[2]))
    df = df.loc[keep].reset_index(drop=True)
    ann = pd.read_excel(os.path.join(BASE, "xenograft_annotation.xlsx"))
    ann["key"] = ann["Metabolite"].map(norm)
    bcaa = set(ann.loc[ann["Sub pathway"] == BCAA_SUB, "key"])
    return df, bcaa


def da_line(df, line, bcaa):
    m = df["Cell line"] == line
    sub = df.loc[m.values, [c for c in df.columns if c not in META]]
    g = df.loc[m.values, "Genotype"].values
    Xz = sub.replace(0, np.nan)
    Xz = Xz.loc[:, Xz.notna().mean(axis=0) >= 0.5]
    half = Xz.min(axis=0, skipna=True) / 2.0
    Xl = np.log10(Xz.fillna(half))
    kl, ol = Xl.loc[g == "K27M"], Xl.loc[g == "KO"]
    rows = []
    for met in Xl.columns:
        a, b = kl[met].values, ol[met].values
        t, p = stats.ttest_ind(a, b, equal_var=False)
        rows.append((met, (a.mean() - b.mean()) / np.log10(2), p))
    r = pd.DataFrame(rows, columns=["metabolite", "log2FC", "p"]).dropna(subset=["p"])
    r["FDR"] = multipletests(r["p"], method="fdr_bh")[1]
    r["sig"] = (r["FDR"] < FDR_THR) & (r["log2FC"].abs() > LFC_THR)
    r["is_bcaa"] = r["metabolite"].map(norm).isin(bcaa)
    return r


def volcano(ax, r, line, xlim):
    x, y = r["log2FC"].values, -np.log10(r["p"].values)
    other = (~r["is_bcaa"]).values
    bc_sig = (r["is_bcaa"] & r["sig"]).values
    bc_ns = (r["is_bcaa"] & ~r["sig"]).values
    ax.scatter(x[other], y[other], s=1.3, c=NS_C, alpha=.6, linewidths=0)
    ax.scatter(x[bc_ns], y[bc_ns], s=3, facecolors="none", edgecolors=BCAA_C,
               linewidths=0.3, alpha=.8)
    ax.scatter(x[bc_sig], y[bc_sig], s=5, c=BCAA_C, edgecolors="k",
               linewidths=0.2, zorder=6)
    ax.axvline(LFC_THR, color="gray", ls="--", lw=.25)
    ax.axvline(-LFC_THR, color="gray", ls="--", lw=.25)
    ax.axhline(-np.log10(0.05), color="gray", ls=":", lw=.3)
    ax.set_xlim(-xlim, xlim)
    # annotate the log2FC cutoff on the dashed threshold lines
    ax.text(LFC_THR, 0.99, f"log2FC = ±{LFC_THR:.0f}", transform=ax.get_xaxis_transform(),
            rotation=90, ha="right", va="top", fontsize=FS - 0.3, color="gray")
    lab = r[r["sig"] & r["is_bcaa"] & (r["log2FC"] > 0)].nlargest(N_LABEL, "log2FC")
    texts = [ax.text(rr["log2FC"], -np.log10(rr["p"]), str(rr["metabolite"]),
                     fontsize=FS, color=BCAA_C) for _, rr in lab.iterrows()]
    adjust_text(texts, ax=ax, x=list(lab["log2FC"]), y=list(-np.log10(lab["p"])),
                expand_points=(1.4, 1.6), force_points=(0.3, 0.5),
                arrowprops=dict(arrowstyle="-", color="#888888", lw=.2))
    ax.set_xlabel("log2 FC (K27M / KO)")
    ax.set_ylabel("-log10 p")
    ax.set_title(f"{line}  (BCAA sig = {int((r['is_bcaa'] & r['sig']).sum())})")
    ax.set_xticks([-8, -4, 0, 4, 8])
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def main():
    df, bcaa = load()
    res = {ln: da_line(df, ln, bcaa) for ln in LINES}
    xlim = float(np.ceil(min(r["log2FC"].abs().max() for r in res.values())))
    fig, axes = plt.subplots(1, 2, figsize=(7 * CM, 3 * CM))
    for ax, ln in zip(axes, LINES):
        volcano(ax, res[ln], ln, xlim)
    handles = [Line2D([0], [0], marker="o", color="w", markerfacecolor=BCAA_C,
                      markeredgecolor="k", markersize=2.5, label="BCAA (sig)"),
               Line2D([0], [0], marker="o", color="w", markerfacecolor="none",
                      markeredgecolor=BCAA_C, markersize=2.5, label="BCAA (n.s.)")]
    fig.legend(handles=handles, loc="upper center", ncol=2, frameon=False,
               handletextpad=0.3, columnspacing=0.8, bbox_to_anchor=(0.5, 1.06))
    fig.tight_layout(w_pad=0.6)
    for ext in ("svg", "png", "pdf"):
        out = os.path.join(FIG, f"xenograft_volcano_BCAA_adjacent_final.{ext}")
        fig.savefig(out, dpi=600, bbox_inches="tight")
        print("wrote", out)
    plt.close(fig)


if __name__ == "__main__":
    main()
