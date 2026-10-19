"""Compare the single-cell pseudobulk DE against the hg38 bulk DE, on the
shared gene universe."""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import numpy as np, pandas as pd
from scipy.stats import spearmanr, pearsonr, fisher_exact
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

B = REPO_ROOT + ""
OUT = f"{B}/scXenograft"; FIG = f"{OUT}/figures"

scr = pd.read_csv(f"{OUT}/sc_pseudobulk_deseq2_K27M_vs_KO.csv", index_col=0)
bk = pd.read_csv(f"{B}/xHSJ019_K27M_vs_KO_analysis/deseq2_K27M_vs_KO.csv", index_col=0)
sc_pb = pd.read_csv(f"{OUT}/sc_pseudobulk_counts_by_library.csv", index_col=0)
bulk_c = pd.read_csv(f"{B}/xHSJ019_K27M_vs_KO_analysis/gene_counts_raw.csv", index_col=0)

sg, bg = set(sc_pb.index), set(bulk_c.index)
print("=== Gene universe (hg19 sc vs hg38 bulk symbols) ===")
print(f"  sc (hg19, after min_cells filter): {len(sg)}")
print(f"  bulk (hg38, all quantified):       {len(bg)}")
print(f"  shared:                            {len(sg & bg)}")
print(f"  sc-only: {len(sg-bg)}    bulk-only: {len(bg-sg)}")

common = scr.index.intersection(bk.index)
print(f"\n  tested in BOTH DE analyses: {len(common)}")
a, b = scr.loc[common], bk.loc[common]

print("\n=== Global concordance of log2FC ===")
ok = np.isfinite(a.log2FoldChange) & np.isfinite(b.log2FoldChange)
print(f"  Pearson  r = {pearsonr(a.log2FoldChange[ok], b.log2FoldChange[ok])[0]:+.3f}")
print(f"  Spearman r = {spearmanr(a.log2FoldChange[ok], b.log2FoldChange[ok]).statistic:+.3f}")
print(f"  Spearman on Wald stat = {spearmanr(a.stat[ok], b.stat[ok]).statistic:+.3f}")
same = np.sign(a.log2FoldChange[ok]) == np.sign(b.log2FoldChange[ok])
print(f"  same direction: {100*same.mean():.1f}% of {ok.sum()} genes")

for th in [0.05, 0.01]:
    sb = (b.padj < th) & (b.log2FoldChange.abs() > 1)
    sub = ok & sb
    same2 = np.sign(a.log2FoldChange[sub]) == np.sign(b.log2FoldChange[sub])
    print(f"  restricted to bulk padj<{th} & |LFC|>1 (n={int(sub.sum())}): "
          f"same direction {100*same2.mean():.1f}%, "
          f"Spearman {spearmanr(a.log2FoldChange[sub], b.log2FoldChange[sub]).statistic:+.3f}")

print("\n=== Overlap of significant calls ===")
A_ = (a.padj < 0.05) & (a.log2FoldChange.abs() > 1)
B_ = (b.padj < 0.05) & (b.log2FoldChange.abs() > 1)
tab = [[int((A_ & B_).sum()), int((A_ & ~B_).sum())],
       [int((~A_ & B_).sum()), int((~A_ & ~B_).sum())]]
orr, p = fisher_exact(tab)
print(f"  sc sig {int(A_.sum())}, bulk sig {int(B_.sum())}, both {tab[0][0]}")
print(f"  Fisher OR {orr:.1f}, p {p:.2e}")
conc = ((A_ & B_) & (np.sign(a.log2FoldChange) == np.sign(b.log2FoldChange))).sum()
print(f"  of the {tab[0][0]} shared calls, {int(conc)} agree in direction "
      f"({100*conc/max(tab[0][0],1):.1f}%)")

print("\n=== Genes of interest ===")
GOI = ["BCAT1", "BCAT2", "SQLE", "FDFT1", "CYP51A1", "MSMO1", "DHCR24", "NSDHL",
       "EBP", "LSS", "SC5D", "HMGCS1", "MVD", "FDPS", "IDI1", "PMVK", "ACAT2",
       "HMGCR", "SREBF1", "SREBF2", "ANXA1", "IRX4", "PEG3", "FOXD1", "PDGFRA",
       "OLIG1", "OLIG2", "MBP", "PLP1", "APOE", "GFAP"]
rows = []
for g in GOI:
    if g in common:
        rows.append(dict(gene=g, sc_LFC=a.log2FoldChange[g], sc_padj=a.padj[g],
                         bulk_LFC=b.log2FoldChange[g], bulk_padj=b.padj[g]))
goi = pd.DataFrame(rows).set_index("gene")
goi["agree"] = np.sign(goi.sc_LFC) == np.sign(goi.bulk_LFC)
print(goi.round(4).to_string())
goi.round(5).to_csv(f"{OUT}/sc_vs_bulk_genes_of_interest.csv")

# mevalonate/sterol axis as a block
MEV = ["SQLE","FDFT1","CYP51A1","MSMO1","DHCR24","NSDHL","EBP","LSS","SC5D",
       "HMGCS1","MVD","FDPS","IDI1","PMVK","ACAT2","HMGCR"]
m = [g for g in MEV if g in common]
print(f"\nMevalonate/sterol axis ({len(m)} genes present in both):")
print(f"  bulk  mean log2FC {b.log2FoldChange[m].mean():+.3f}  "
      f"({int((b.log2FoldChange[m]>0).sum())}/{len(m)} up in K27M)")
print(f"  sc    mean log2FC {a.log2FoldChange[m].mean():+.3f}  "
      f"({int((a.log2FoldChange[m]>0).sum())}/{len(m)} up in K27M)")
from scipy.stats import wilcoxon
print(f"  sc vs 0: Wilcoxon p {wilcoxon(a.log2FoldChange[m]).pvalue:.4f}")

merged = a.join(b, lsuffix="_sc", rsuffix="_bulk")
merged.round(5).to_csv(f"{OUT}/sc_vs_bulk_deseq2_merged.csv")

# ---------- figure ----------
fig, axes = plt.subplots(1, 3, figsize=(16, 5.2))
ax = axes[0]
ax.scatter(b.log2FoldChange[ok], a.log2FoldChange[ok], s=4, alpha=.2,
           c="#888888", linewidths=0, rasterized=True)
sel = ok & B_
ax.scatter(b.log2FoldChange[sel], a.log2FoldChange[sel], s=6, alpha=.45,
           c="#C1272D", linewidths=0, rasterized=True, label="bulk padj<0.05 & |LFC|>1")
lim = 10
ax.plot([-lim, lim], [-lim, lim], lw=.8, ls="--", c="k")
ax.axhline(0, lw=.5, c="k"); ax.axvline(0, lw=.5, c="k")
ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
ax.set_xlabel("bulk log2FC (hg38, n=3v3)")
ax.set_ylabel("sc pseudobulk log2FC (hg19, n=2v2)")
ax.set_title(f"log2FC concordance\nPearson r = "
             f"{pearsonr(a.log2FoldChange[ok], b.log2FoldChange[ok])[0]:+.2f}", fontsize=10)
ax.legend(fontsize=7, frameon=False, markerscale=2)

ax = axes[1]
for g in m:
    ax.annotate(g, (b.log2FoldChange[g], a.log2FoldChange[g]), fontsize=6.5)
ax.scatter(b.log2FoldChange[m], a.log2FoldChange[m], s=40, c="#2B6CB0")
if "BCAT1" in common:
    ax.scatter([b.log2FoldChange["BCAT1"]], [a.log2FoldChange["BCAT1"]],
               s=90, c="#C1272D", marker="D", zorder=5)
    ax.annotate("BCAT1", (b.log2FoldChange["BCAT1"], a.log2FoldChange["BCAT1"]),
                fontsize=9, color="#C1272D", fontweight="bold",
                xytext=(5, -10), textcoords="offset points")
lo, hi = -1.5, 2.6
ax.plot([lo, hi], [lo, hi], lw=.8, ls="--", c="k")
ax.axhline(0, lw=.5, c="k"); ax.axvline(0, lw=.5, c="k")
ax.set_xlabel("bulk log2FC"); ax.set_ylabel("sc pseudobulk log2FC")
ax.set_title("mevalonate / sterol axis + BCAT1", fontsize=10)

ax = axes[2]
v = a.log2FoldChange
p = -np.log10(a.padj.clip(lower=1e-300))
sigm = (a.padj < 0.05) & (v.abs() > 1)
ax.scatter(v[~sigm], p[~sigm], s=4, c="#cccccc", linewidths=0, rasterized=True)
ax.scatter(v[sigm], p[sigm], s=6, c="#C1272D", linewidths=0, rasterized=True)
ax.set_xlabel("log2FC (K27M vs KO)"); ax.set_ylabel("-log10 padj")
ax.set_title(f"sc pseudobulk volcano\n{int(sigm.sum())} genes padj<0.05 & |LFC|>1",
             fontsize=10)
for a_ in axes:
    for sp in ["top", "right"]: a_.spines[sp].set_visible(False)
fig.suptitle("HSJ019 single-cell pseudobulk vs bulk tissue DE, shared gene universe",
             fontsize=12)
fig.tight_layout(); fig.savefig(f"{FIG}/hsj019_scpseudobulk_vs_bulk.png", dpi=170)
print("\nfigure written")
