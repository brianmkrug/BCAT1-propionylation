"""BCAA metabolism gene-level paired K27M | KO heatmap along pooled pseudotime
(4-library pooled, matching the TCA/OxPhos figure). Genes ordered by catabolic
step from BCAA-degradation + propanoate-metabolism KEGG sets. The propionyl-CoA
degradative node (PCCA, PCCB, MMUT, MCEE) is asterisked.
Pooled z per gene; right strip = pseudobulk log2FC (K27M/KO)."""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import numpy as np, pandas as pd, scanpy as sc
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

B = REPO_ROOT + ""; W = f"{B}/scXenograft/work"
OUT = f"{B}/scXenograft"; FIG = f"{OUT}/figures"
NB = 20; DET = 12
ALIAS = {"MMUT": "MUT"}   # hg19 symbol
STAR = {"PCCA", "PCCB", "MMUT", "MCEE"}   # propionyl-CoA -> succinyl-CoA node

# curated catabolic order (union of BCAA degradation + propanoate metabolism)
GROUPS = [
 ("Transamination (BCAA→BCKA)", ["BCAT1","BCAT2"]),
 ("BCKDH complex (BCKA→acyl-CoA)", ["BCKDHA","BCKDHB","DBT","DLD"]),
 ("Acyl-CoA dehydrogenases", ["IVD","ACAD8","ACADSB","ACADM","ACADS"]),
 ("Leucine → ketogenic (HMG-CoA)", ["MCCC1","MCCC2","AUH","HMGCL","HMGCLL1","HMGCS1",
                                    "OXCT1","OXCT2","ACAT1","ACAT2"]),
 ("Val/Ile → CoA (hydratase/aldehyde dh)", ["ECHS1","HIBCH","HIBADH","ALDH6A1","ABAT",
        "HADH","HADHA","HADHB","HSD17B10","ACAA1","ACAA2","EHHADH","AACS","ALDH2",
        "ALDH3A2","ALDH7A1","ALDH9A1","ALDH1B1","AOX1","IL4I1","ACSF3","AGXT2"]),
 ("Propionyl-CoA → succinyl-CoA (anaplerotic) *", ["PCCA","PCCB","MCEE","MMUT"]),
 ("Propanoate / other", ["SUCLG1","SUCLG2","SUCLA2","ACSS1","ACSS2","ACSS3","ACACA",
        "ACACB","MLYCD","ACOX1","ACOX3","ECHDC1","LDHA","LDHB","LDHC"]),
]

emb = pd.read_csv(f"{W}/pooled_dpt.csv.gz", index_col=0)
emb["bin"] = pd.qcut(emb.dpt, NB, labels=False, duplicates="drop")
A = sc.read_h5ad(f"{W}/hsj019_norm.h5ad")[emb.index]
de = pd.read_csv(f"{OUT}/sc_pseudobulk_deseq2_K27M_vs_KO.csv", index_col=0)
vset = set(A.var_names)
def ref(g): return ALIAS.get(g, g)
def det(g):
    r = ref(g)
    return 100*(np.asarray(A[:, r].X.todense()).ravel() > 0).mean() if r in vset else 0

ordered, labels, gseen, off = [], [], [], 0
seen = set()
for gname, genes in GROUPS:
    keep = []
    for g in genes:
        r = ref(g)
        if r in vset and det(g) >= DET and r not in seen:
            keep.append(r); seen.add(r)
            lab = ("MMUT/MUT" if g == "MMUT" else g) + (" ✱" if g in STAR else "")
            labels.append(lab)
    if keep:
        gseen.append((gname, off, off+len(keep))); off += len(keep)
        ordered += keep
star_group = [n for n, _ in GROUPS].index("Propionyl-CoA → succinyl-CoA (anaplerotic) *")

E = pd.DataFrame(np.asarray(A[:, ordered].X.todense()), index=emb.index, columns=ordered)
def zmat(gt):
    return np.vstack([[E[g].values[(emb.genotype==gt).values & (emb.bin==b).values].mean()
                       for b in range(NB)] for g in ordered])
Kr, Or = zmat("K27M"), zmat("K27M-KO")
mu = np.array([np.r_[Kr[i],Or[i]].mean() for i in range(len(ordered))])[:,None]
sd = np.array([np.r_[Kr[i],Or[i]].std() for i in range(len(ordered))])[:,None]
K, O = (Kr-mu)/sd, (Or-mu)/sd
lfc = np.array([de.log2FoldChange.get(g, np.nan) for g in ordered])

VZ = 2.5
fig = plt.figure(figsize=(14.5, max(11, 0.26*len(ordered))))
gs = GridSpec(1, 4, width_ratios=[1,1,0.06,0.045], wspace=0.05)
axK, axO = fig.add_subplot(gs[0,0]), fig.add_subplot(gs[0,1])
for ax, Mx, ttl in [(axK,K,"K27M"),(axO,O,"K27M-KO (r1+r3)")]:
    im = ax.imshow(Mx, aspect="auto", cmap="RdBu_r", vmin=-VZ, vmax=VZ)
    ax.set_xticks([0,NB//2,NB-1]); ax.set_xticklabels(["early","mid","late"]); ax.set_xlabel("pooled pseudotime bin")
    ax.set_title(ttl, fontsize=11)
    for _,a,b in gseen:
        if a>0: ax.axhline(a-0.5, color="k", lw=0.8)
    # box the asterisked node
    a0, a1 = gseen[star_group][1], gseen[star_group][2]
    ax.add_patch(plt.Rectangle((-0.5, a0-0.5), NB, a1-a0, fill=False, edgecolor="#111", lw=2.2))
axK.set_yticks(range(len(ordered))); axK.set_yticklabels(labels, fontsize=7); axO.set_yticks([])
for t, lab in zip(axK.get_yticklabels(), labels):
    if "✱" in lab: t.set_fontweight("bold"); t.set_color("#7a0a0a")
axL = fig.add_subplot(gs[0,2]); axL.imshow(lfc[:,None], aspect="auto", cmap="PuOr_r", vmin=-1.5, vmax=1.5)
axL.set_xticks([]); axL.set_yticks([]); axL.set_title("log2FC\nK27M/KO", fontsize=6.5)
for _,a,b in gseen:
    if a>0: axL.axhline(a-0.5, color="k", lw=0.8)
fig.colorbar(im, cax=fig.add_subplot(gs[0,3]), label="pooled z (per gene)")
for gname,a,b in gseen:
    axK.text(-0.085*NB, (a+b-1)/2, gname, transform=axK.get_yaxis_transform(),
             ha="right", va="center", fontsize=6.6, style="italic",
             color="#7a0a0a" if "*" in gname else "#333",
             fontweight="bold" if "*" in gname else "normal")
fig.suptitle("BCAA metabolism — gene level, K27M vs KO across pooled pseudotime\n"
             "pooled z per gene (shared scale); right strip = pseudobulk log2FC; "
             "✱ / boxed = propionyl-CoA → succinyl-CoA node (PCCA, PCCB, MMUT, MCEE)",
             fontsize=11)
fig.subplots_adjust(left=0.26, top=0.95)
fig.savefig(f"{FIG}/bcaa_gene_heatmap.png", dpi=170, bbox_inches="tight"); plt.close(fig)

outdf = pd.DataFrame(K, index=ordered, columns=[f"K27M_bin{b}" for b in range(NB)]).join(
        pd.DataFrame(O, index=ordered, columns=[f"KO_bin{b}" for b in range(NB)]))
outdf["group"] = [g for g,a,b in gseen for _ in range(b-a)]; outdf["pb_log2FC"] = lfc
outdf["propionyl_node"] = [lab.replace(" ✱","") in {"PCCA","PCCB","MCEE","MMUT/MUT"} and "✱" in lab for lab in labels]
outdf.round(3).to_csv(f"{OUT}/bcaa_gene_heatmap_values.csv")
print(f"{len(ordered)} genes shown (>= {DET}% detected)")
print("dropped (low detection):", [ref(g) for _,gg in GROUPS for g in gg
      if ref(g) in vset and det(g) < DET])
print("\npropionyl-CoA node (pooled-z K27M vs KO across bins, + log2FC):")
for g in ["PCCA","PCCB","MUT","MCEE"]:
    if g in ordered:
        i = ordered.index(g)
        print(f"  {g:<6} K27M {K[i].mean():+.2f}  KO {O[i].mean():+.2f}  "
              f"Δ {K[i].mean()-O[i].mean():+.2f}  log2FC {lfc[i]:+.2f}  det {det(g):.0f}%")
print("\ngenotype lean by group (mean pooled-z K27M-KO):")
for name,a,b in gseen:
    print(f"  {name:<44} n={b-a:<2} {np.nanmean(K[a:b]-O[a:b]):+.2f}")
