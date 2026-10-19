#!/usr/bin/env python3
"""
Pentose & glucuronate interconversions (KEGG hsa00040) — the top positively
selected pathway in the DIPGXIII K27M arm.

Panel A  bar chart of K27M beta for all 18 pathway genes present in the screen,
         shaded by whether the gene is actually expressed in DIPGXIII RNA-seq.
Panel B  K27M vs KO beta for the same genes.
Panel C  schematic of the two arms of the pathway.

Outputs fig_pentose_glucuronate.{png,pdf,svg} + pentose_glucuronate_table.csv
"""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Patch

SCR = REPO_ROOT + "/CRISPR screen"
RNA = REPO_ROOT + "/xDIPGXIII_K27M_vs_KO_analysis/gene_counts_raw_clean.csv"
OUT = f"{SCR}/depmap_classification"
CPM_MIN = 1.0

d = pd.read_excel(f"{SCR}/beta_correlation.xlsx")
gmt = pd.read_pickle(f"{SCR}/gsea_out/gmt.pkl")
path = gmt[[k for k in gmt if "entose" in k][0]]

cnt = pd.read_csv(RNA, index_col=0)
cpm = (cnt / cnt.sum() * 1e6).mean(axis=1)
d["CPM"] = d.Gene.map(cpm)
d["rank_pos"] = d["K27M.beta"].rank(ascending=False).astype(int)

sub = d[d.Gene.isin(path)].sort_values("K27M.beta").reset_index(drop=True)
sub["expressed"] = sub.CPM > CPM_MIN
sub.sort_values("K27M.beta", ascending=False).to_csv(
    f"{OUT}/pentose_glucuronate_table.csv", index=False)

C_EXPR, C_NOT = "#b2182b", "#dcdcdc"

fig = plt.figure(figsize=(13.5, 7.4))
gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.85], width_ratios=[1.0, 1.0],
                      hspace=0.42, wspace=0.26)

# ---------------- A: beta bars ----------------
axa = fig.add_subplot(gs[:, 0])
y = np.arange(len(sub))
axa.barh(y, sub["K27M.beta"], color=np.where(sub.expressed, C_EXPR, C_NOT),
         height=0.74, linewidth=0)
axa.axvline(0, color="black", lw=0.8)
axa.set_yticks(y)
axa.set_yticklabels(sub.Gene, fontsize=8.6)
for t, e in zip(axa.get_yticklabels(), sub.expressed):
    if not e:
        t.set_color("#8a8a8a"); t.set_style("italic")
for i, r in sub.iterrows():
    x = r["K27M.beta"]
    axa.text(x + (0.03 if x >= 0 else -0.03), i,
             f"{r.CPM:.1f}" if r.CPM >= 1 else f"{r.CPM:.2f}",
             va="center", ha="left" if x >= 0 else "right",
             fontsize=6.6, color="#555555")
axa.set_xlim(-0.75, 1.22)
axa.set_xlabel("K27M MAGeCK MLE $\\beta$   (grey numerals = mean CPM in DIPGXIII)", fontsize=9)
axa.set_title("Pentose & glucuronate interconversions\nranked by K27M $\\beta$",
              fontsize=11, fontweight="bold")
axa.spines[["top", "right", "left"]].set_visible(False)
axa.tick_params(axis="x", labelsize=8.5)
axa.legend(handles=[Patch(facecolor=C_EXPR, label=f"expressed (CPM > {CPM_MIN:.0f})"),
                    Patch(facecolor=C_NOT, label="not expressed — β uninterpretable")],
           fontsize=8, loc="lower right", frameon=False)

# ---------------- B: K27M vs KO ----------------
axb = fig.add_subplot(gs[0, 1])
axb.axhline(0, color="#bbbbbb", lw=0.7)
axb.axvline(0, color="#bbbbbb", lw=0.7)
axb.scatter(sub["KO.beta"], sub["K27M.beta"],
            c=np.where(sub.expressed, C_EXPR, C_NOT), s=42,
            edgecolor="#666666", linewidth=0.4, zorder=3)
for _, r in sub.iterrows():
    if r.expressed and (r["K27M.beta"] > 0.3 or abs(r.difference) > 0.8):
        axb.annotate(r.Gene, (r["KO.beta"], r["K27M.beta"]), xytext=(4, 3),
                     textcoords="offset points", fontsize=7.6)
lim = [-2.1, 1.5]
axb.plot(lim, lim, ls="--", color="#999999", lw=0.7)
axb.set_xlim(lim); axb.set_ylim(-0.75, 1.1)
axb.set_xlabel("KO (H3-WT isogenic) $\\beta$", fontsize=9)
axb.set_ylabel("K27M $\\beta$", fontsize=9)
axb.set_title("K27M vs isogenic KO", fontsize=10, fontweight="bold")
axb.spines[["top", "right"]].set_visible(False)
axb.tick_params(labelsize=8)

# ---------------- C: schematic ----------------
axc = fig.add_subplot(gs[1, 1])
axc.set_xlim(0.1, 10.55); axc.set_ylim(0.2, 5.4); axc.axis("off")
axc.set_title("Two arms of the pathway", fontsize=10, fontweight="bold", y=0.99)

beta = dict(zip(sub.Gene, sub["K27M.beta"]))


def node(x, y, label, w=1.6, h=0.52, fc="#f2f2f2"):
    axc.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h,
                                 boxstyle="round,pad=0.05", fc=fc, ec="#888888", lw=0.7))
    axc.text(x, y, label, ha="center", va="center", fontsize=7.4)


def arrow(x1, y1, x2, y2, gene, dy=0.20):
    axc.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                                  mutation_scale=9, lw=1.0, color="#555555"))
    b = beta.get(gene, np.nan)
    col = C_EXPR if (gene in beta and sub.set_index("Gene").loc[gene, "expressed"]) else "#8a8a8a"
    axc.text((x1 + x2) / 2, (y1 + y2) / 2 + dy, f"{gene}\nβ={b:+.2f}",
             ha="center", va="bottom", fontsize=6.8, color=col, fontweight="bold")


# upper arm: UDP-sugar / GAG
node(1.2, 4.2, "Glc-1-P"); node(4.0, 4.2, "UDP-glucose"); node(7.0, 4.2, "UDP-glucuronate")
node(9.45, 4.2, "GAGs / PGs\nglucuronidation", w=1.9, h=0.72, fc="#fdece8")
arrow(2.0, 4.2, 3.2, 4.2, "UGP2")
arrow(4.8, 4.2, 6.2, 4.2, "UGDH")
axc.add_patch(FancyArrowPatch((7.8, 4.2), (8.45, 4.2), arrowstyle="-|>",
                              mutation_scale=9, lw=1.0, color="#555555"))

# lower arm: uronate / glucuronate cycle back into PPP
node(1.2, 1.5, "D-glucuronate"); node(4.0, 1.5, "L-xylulose")
node(6.6, 1.5, "xylitol"); node(9.45, 1.5, "D-xylulose-5-P\n→ PPP", w=1.9, h=0.72, fc="#fdece8")
arrow(2.0, 1.5, 3.2, 1.5, "CRYL1")
arrow(4.8, 1.5, 5.8, 1.5, "DCXR")
arrow(7.4, 1.5, 8.45, 1.5, "SORD")
axc.text(6.6, 0.72, "XYLB  β=%+.2f" % beta.get("XYLB", np.nan),
         ha="center", fontsize=6.8, color=C_EXPR, fontweight="bold")
axc.text(4.9, 2.55, "NADPH consumed", fontsize=6.6, color="#1a5276", style="italic")
axc.text(7.3, 2.15, "NAD$^+$ consumed", fontsize=6.6, color="#1a5276", style="italic")
axc.text(5.4, 4.86, "2 NAD$^+$ consumed", fontsize=6.6, color="#1a5276", style="italic")
axc.annotate("", xy=(7.0, 3.9), xytext=(1.4, 1.85),
             arrowprops=dict(arrowstyle="-|>", lw=0.8, color="#aaaaaa", ls=":"))
axc.text(2.4, 2.9, "GUSB β=%+.2f\n(lysosomal release)" % beta.get("GUSB", np.nan),
         fontsize=6.6, color="#777777")

fig.text(0.005, 0.005,
         "GSEA prerank (K27M β): NES 1.73, nominal p = 0.0014, FDR q = 0.16 — top positively selected KEGG pathway of 76 tested. "
         "Mann–Whitney of pathway vs all screened genes: p = 0.0065 (all 18 genes), p = 0.043 restricted to expressed genes. "
         "Expression = mean CPM across 5 DIPGXIII RNA-seq libraries.",
         fontsize=6.6, color="#555555")

fig.subplots_adjust(left=0.075, right=0.985, top=0.92, bottom=0.075)
for ext in ("png", "pdf", "svg"):
    fig.savefig(f"{OUT}/fig_pentose_glucuronate.{ext}",
                dpi=300 if ext == "png" else None, bbox_inches="tight")
print("saved")
