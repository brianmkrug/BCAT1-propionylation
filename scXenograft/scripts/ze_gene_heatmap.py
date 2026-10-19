"""Gene-level paired K27M | KO heatmaps for cholesterol biosynthesis and
handling, genes ordered by biological stage. Per-gene metric = POOLED z-score
of (genotype,bin) mean log-norm expression (both genotypes share the center/scale
per gene). Right strip = sc-pseudobulk log2FC (K27M vs KO). Genes with <12%
detection dropped (z would be dropout noise)."""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import numpy as np, pandas as pd, scanpy as sc
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

B = REPO_ROOT + ""
W = f"{B}/scXenograft/work"; OUT = f"{B}/scXenograft"; FIG = f"{OUT}/figures"
NB = 20; DET_MIN = 12

# curated order by biological stage; (group, [genes])
GROUPS = [
    ("Mevalonate arm", ["ACAT2","HMGCS1","HMGCR","MVK","PMVK","MVD"]),
    ("Isoprenoid / FPP", ["IDI1","FDPS","GGPS1"]),
    ("Squalene", ["FDFT1","SQLE"]),
    ("Post-squalene sterol ring", ["LSS","TM7SF2","MSMO1","NSDHL","HSD17B7",
                                    "SC5D","DHCR7","DHCR24","EBP","LBR","ARV1"]),
    ("Import / receptors", ["LDLR","SCARB1","LRP1","NPC1","NPC2","SORT1"]),
    ("Esterification / hydrolysis", ["SOAT1","LIPA"]),
    ("Efflux", ["ABCA1"]),
    ("Apolipoprotein / lipoprotein", ["APOE","APOC1","PLTP","LPL"]),
    ("Intracellular transport / reg.", ["MYLIP","OSBPL5","STARD3","VAPA","VAPB","VDAC1"]),
]
BLOCK = [("BIOSYNTHESIS (de novo)", ["Mevalonate arm","Isoprenoid / FPP","Squalene",
                                     "Post-squalene sterol ring"]),
         ("HANDLING (uptake / efflux / transport)",
          ["Import / receptors","Esterification / hydrolysis","Efflux",
           "Apolipoprotein / lipoprotein","Intracellular transport / reg."])]

emb = pd.read_csv(f"{W}/pooled_dpt.csv.gz", index_col=0)
emb["bin"] = pd.qcut(emb.dpt, NB, labels=False, duplicates="drop")
A = sc.read_h5ad(f"{W}/hsj019_norm.h5ad")[emb.index]
de = pd.read_csv(f"{OUT}/sc_pseudobulk_deseq2_K27M_vs_KO.csv", index_col=0)
isK = (emb.genotype == "K27M").values

ordered, gsub, glab = [], [], []
for gname, genes in GROUPS:
    for g in genes:
        if g in A.var_names:
            v = np.asarray(A[:, g].X.todense()).ravel()
            if 100*(v > 0).mean() >= DET_MIN:
                ordered.append(g); gsub.append(gname); glab.append(g)
E = pd.DataFrame(np.asarray(A[:, ordered].X.todense()), index=emb.index, columns=ordered)

def zmat(gt):
    m = np.vstack([[E[g].values[(emb.genotype==gt).values & (emb.bin==b).values].mean()
                    for b in range(NB)] for g in ordered])
    return m
Kr, Or = zmat("K27M"), zmat("K27M-KO")
mu = np.array([np.r_[Kr[i], Or[i]].mean() for i in range(len(ordered))])[:, None]
sd = np.array([np.r_[Kr[i], Or[i]].std() for i in range(len(ordered))])[:, None]
K = (Kr - mu)/sd; O = (Or - mu)/sd
lfc = np.array([de.log2FoldChange.get(g, np.nan) for g in ordered])

# group boundaries within the ordered list
bounds, off = [], 0
gseen = []
for gname, genes in GROUPS:
    present = [g for g in genes if g in ordered]
    if present:
        gseen.append((gname, off, off+len(present)))
        off += len(present)

VZ = 2.5
fig = plt.figure(figsize=(15.5, 13))
gs = GridSpec(1, 4, width_ratios=[1, 1, 0.06, 0.045], wspace=0.05)
axK = fig.add_subplot(gs[0, 0]); axO = fig.add_subplot(gs[0, 1])
for ax, Mx, ttl in [(axK, K, "K27M"), (axO, O, "K27M-KO (r1+r3)")]:
    im = ax.imshow(Mx, aspect="auto", cmap="RdBu_r", vmin=-VZ, vmax=VZ)
    ax.set_xticks([0, NB//2, NB-1]); ax.set_xticklabels(["early","mid","late"])
    ax.set_xlabel("pooled pseudotime bin"); ax.set_title(ttl, fontsize=11)
    for _, a, bnd in gseen:
        if a > 0: ax.axhline(a-0.5, color="k", lw=1.0)
axK.set_yticks(range(len(ordered))); axK.set_yticklabels(glab, fontsize=7.5)
axO.set_yticks([])
# log2FC strip
axL = fig.add_subplot(gs[0, 2])
lm = axL.imshow(lfc[:, None], aspect="auto", cmap="PuOr_r",
                vmin=-2, vmax=2)
axL.set_xticks([]); axL.set_yticks([]); axL.set_title("log2FC\nK27M/KO", fontsize=7)
for _, a, bnd in gseen:
    if a > 0: axL.axhline(a-0.5, color="k", lw=1.0)
fig.colorbar(im, cax=fig.add_subplot(gs[0, 3]), label="pooled z (per gene)")
# stage-group labels (left) and biosynthesis/handling brackets
for gname, a, bnd in gseen:
    axK.text(-0.055*NB, (a+bnd-1)/2, gname, transform=axK.transData if False else axK.get_yaxis_transform(),
             ha="right", va="center", fontsize=7.2, style="italic", color="#333")
fig.suptitle("Cholesterol biosynthesis and handling — gene level, K27M vs KO across pooled pseudotime\n"
             "pooled z per gene (shared scale); right strip = pseudobulk log2FC; genes >=12% detected, "
             "ordered by pathway stage", fontsize=11)
fig.subplots_adjust(left=0.22, top=0.93)
fig.savefig(f"{FIG}/cholesterol_gene_heatmap.png", dpi=170, bbox_inches="tight")
plt.close(fig)

# export
outdf = pd.DataFrame(K, index=ordered, columns=[f"K27M_bin{b}" for b in range(NB)]).join(
        pd.DataFrame(O, index=ordered, columns=[f"KO_bin{b}" for b in range(NB)]))
outdf["stage"] = gsub; outdf["pb_log2FC"] = lfc
outdf.round(3).to_csv(f"{OUT}/cholesterol_gene_heatmap_values.csv")
print(f"{len(ordered)} genes shown (>= {DET_MIN}% detected)")
print("dropped (low detection):",
      [g for _,gg in GROUPS for g in gg if g in A.var_names and g not in ordered])
print("\nby-stage genotype lean (mean pooled-z K27M minus KO across bins):")
for gname, a, bnd in gseen:
    d = float(np.nanmean(K[a:bnd] - O[a:bnd]))
    print(f"  {gname:<32} {d:+.2f}  ({'K27M' if d>0 else 'KO'} higher)")
