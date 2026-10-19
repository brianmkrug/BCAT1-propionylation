"""Gene-level paired K27M | KO heatmap for TCA cycle + OxPhos, ordered by
cycle step / ETC complex. Pooled z per gene; log2FC strip. Large exploratory list."""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import re
import numpy as np, pandas as pd, scanpy as sc
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

B = REPO_ROOT + ""
W = f"{B}/scXenograft/work"; OUT = f"{B}/scXenograft"; FIG = f"{OUT}/figures"
NB = 20; DET = 15; CAP = 14   # cap the biggest complexes for legibility

def read_set(name):
    for ln in open(f"{B}/xHSJ019_K27M_vs_KO_analysis/KEGG_metabolism_expanded.gmt"):
        f = ln.rstrip("\n").split("\t")
        if f[0] == name: return [x for x in f[2:] if x]
    return []
TCA = read_set("Citrate cycle (TCA cycle)"); OXP = read_set("Oxidative phosphorylation")

emb = pd.read_csv(f"{W}/pooled_dpt.csv.gz", index_col=0)
emb["bin"] = pd.qcut(emb.dpt, NB, labels=False, duplicates="drop")
A = sc.read_h5ad(f"{W}/hsj019_norm.h5ad")[emb.index]

# hg38/2016 -> hg19 symbol aliases for F-type ATP synthase (renamed 2016) and misc
ALIAS = {"ATP5F1A":"ATP5A1","ATP5F1B":"ATP5B","ATP5F1C":"ATP5C1","ATP5F1D":"ATP5D",
         "ATP5F1E":"ATP5E","ATP5MC1":"ATP5G1","ATP5MC2":"ATP5G2","ATP5MC3":"ATP5G3",
         "ATP5ME":"ATP5I","ATP5MF":"ATP5J2","ATP5MG":"ATP5L","ATP5PB":"ATP5F1",
         "ATP5PD":"ATP5H","ATP5PF":"ATP5J","ATP5PO":"ATP5O","ATP5MPL":"USMG5"}
vset = set(A.var_names)
def remap(gl):
    out = []
    for g in gl:
        if g in vset: out.append(g)
        elif g in ALIAS and ALIAS[g] in vset: out.append(ALIAS[g])
    return out
OXP = remap(OXP)
de = pd.read_csv(f"{OUT}/sc_pseudobulk_deseq2_K27M_vs_KO.csv", index_col=0)
isK = (emb.genotype == "K27M").values
def det(g): return 100*(np.asarray(A[:, g].X.todense()).ravel() > 0).mean() if g in A.var_names else 0

# ---- TCA ordered by step ----
TCA_STEP = [
 ("PDH feed-in", ["PDHA1","PDHB","DLAT","DLD","PDHX","PC"]),
 ("Citrate synthase", ["CS"]),
 ("Aconitase", ["ACO2","ACO1"]),
 ("Isocitrate dh (NAD, CIII-feed)", ["IDH3A","IDH3B","IDH3G","IDH2","IDH1"]),
 ("a-KG dh", ["OGDH","DLST"]),
 ("Succinyl-CoA synth", ["SUCLG1","SUCLA2","SUCLG2"]),
 ("Succinate dh (CII)", ["SDHA","SDHB","SDHC","SDHD"]),
 ("Fumarase", ["FH"]),
 ("Malate dh", ["MDH2","MDH1"]),
 ("Anaplerosis/export", ["PC","PCK2","ACLY"]),
]
# ---- OxPhos ordered by complex ----
def oxc(g):
    if g.startswith("NDUF"): return "CI (NADH dh)"
    if g.startswith("SDH"): return "CII (SDH)"
    if g.startswith("UQCR") or g in ("CYC1",): return "CIII (bc1)"
    if g.startswith("COX") or g.startswith("COA") or g=="NDUFA4": return "CIV (COX)"
    if g.startswith("ATP5"): return "CV (ATP synthase)"
    if g.startswith("ATP6V") or g.startswith("ATP6AP") or g=="TCIRG1": return "V-ATPase (lysosomal)"
    return "other"
# mtDNA-encoded ETC subunits (present in hg19 ref, absent from the KEGG set) -
# added directly for completeness, shown as their own block per complex
MT_BY_COMPLEX = {
 "CI (NADH dh)": ["MT-ND1","MT-ND2","MT-ND3","MT-ND4","MT-ND4L","MT-ND5","MT-ND6"],
 "CIII (bc1)": ["MT-CYB"],
 "CIV (COX)": ["MT-CO1","MT-CO2","MT-CO3"],
 "CV (ATP synthase)": ["MT-ATP6","MT-ATP8"],
}
MT_ALL = [g for v in MT_BY_COMPLEX.values() for g in v]
# CII/SDH shown under TCA (Succinate dh step), so excluded from the OxPhos block
OX_ORDER = ["CI (NADH dh)","CIII (bc1)","CIV (COX)","CV (ATP synthase)",
            "V-ATPase (lysosomal)"]
ox_by = {c: [] for c in OX_ORDER}
for g in OXP:
    c = oxc(g)
    if c in ox_by and det(g) >= DET: ox_by[c].append(g)
for c in ox_by:                       # order by detection, cap the big ones
    ox_by[c] = sorted(ox_by[c], key=det, reverse=True)
    if c in ("CI (NADH dh)","CIV (COX)"): ox_by[c] = ox_by[c][:CAP]
# append mtDNA-encoded subunits to their complex (they are highly expressed)
for c, mts in MT_BY_COMPLEX.items():
    ox_by[c] = ox_by[c] + [g for g in mts if g in vset]
print("CV (ATP synthase) genes:", ox_by["CV (ATP synthase)"])

GROUPS = []
seen = set()
for name, genes in TCA_STEP:
    keep = [g for g in genes if det(g) >= 12 and g not in seen]
    for g in keep: seen.add(g)
    if keep: GROUPS.append((name, keep))
TCA_END = len(GROUPS)
for c in OX_ORDER:
    if ox_by[c]: GROUPS.append((c, ox_by[c]))

# defensive dedup preserving order
ordered, _s = [], set()
for _, gg in GROUPS:
    for g in gg:
        if g not in _s: ordered.append(g); _s.add(g)
E = pd.DataFrame(np.asarray(A[:, ordered].X.todense()), index=emb.index, columns=ordered)
def zmat(gt):
    return np.vstack([[E[g].values[(emb.genotype==gt).values & (emb.bin==b).values].mean()
                       for b in range(NB)] for g in ordered])
Kr, Or = zmat("K27M"), zmat("K27M-KO")
mu = np.array([np.r_[Kr[i],Or[i]].mean() for i in range(len(ordered))])[:,None]
sd = np.array([np.r_[Kr[i],Or[i]].std() for i in range(len(ordered))])[:,None]
K, O = (Kr-mu)/sd, (Or-mu)/sd
lfc = np.array([de.log2FoldChange.get(g,np.nan) for g in ordered])

bounds, off, gseen = [], 0, []
for name, genes in GROUPS:
    gseen.append((name, off, off+len(genes))); off += len(genes)
tca_rows = sum(len(gg) for _,gg in GROUPS[:TCA_END])

VZ = 2.5
fig = plt.figure(figsize=(15.5, max(12, 0.23*len(ordered))))
gs = GridSpec(1, 4, width_ratios=[1,1,0.05,0.04], wspace=0.05)
axK, axO = fig.add_subplot(gs[0,0]), fig.add_subplot(gs[0,1])
for ax, Mx, ttl in [(axK,K,"K27M"),(axO,O,"K27M-KO (r1+r3)")]:
    im = ax.imshow(Mx, aspect="auto", cmap="RdBu_r", vmin=-VZ, vmax=VZ)
    ax.set_xticks([0,NB//2,NB-1]); ax.set_xticklabels(["early","mid","late"])
    ax.set_xlabel("pooled pseudotime bin"); ax.set_title(ttl, fontsize=11)
    for _,a,b in gseen:
        if a>0: ax.axhline(a-0.5, color="k", lw=0.7)
    ax.axhline(tca_rows-0.5, color="k", lw=2.2)      # TCA | OxPhos divider
axK.set_yticks(range(len(ordered))); axK.set_yticklabels(ordered, fontsize=6.5); axO.set_yticks([])
for t, g in zip(axK.get_yticklabels(), ordered):   # bold + red the mtDNA-encoded rows
    if g.startswith("MT-"): t.set_fontweight("bold"); t.set_color("#7a0a0a")
axL = fig.add_subplot(gs[0,2]); axL.imshow(lfc[:,None], aspect="auto", cmap="PuOr_r", vmin=-1.5, vmax=1.5)
axL.set_xticks([]); axL.set_yticks([]); axL.set_title("log2FC\nK27M/KO", fontsize=6.5)
for _,a,b in gseen:
    if a>0: axL.axhline(a-0.5, color="k", lw=0.7)
axL.axhline(tca_rows-0.5, color="k", lw=2.2)
fig.colorbar(im, cax=fig.add_subplot(gs[0,3]), label="pooled z (per gene)")
for name,a,b in gseen:
    axK.text(-0.075*NB, (a+b-1)/2, name, transform=axK.get_yaxis_transform(),
             ha="right", va="center", fontsize=6.8, style="italic", color="#333")
axK.text(-0.16*NB, (tca_rows-1)/2, "TCA CYCLE", transform=axK.get_yaxis_transform(),
         rotation=90, ha="right", va="center", fontsize=10, fontweight="bold", color="#C1272D")
axK.text(-0.16*NB, (tca_rows+len(ordered)-1)/2, "OXPHOS / ETC", transform=axK.get_yaxis_transform(),
         rotation=90, ha="right", va="center", fontsize=10, fontweight="bold", color="#2B6CB0")
fig.suptitle("TCA cycle and oxidative phosphorylation — gene level, K27M vs KO across pooled pseudotime\n"
             "pooled z per gene (shared scale); right strip = pseudobulk log2FC; ordered by cycle step / ETC complex; "
             "bold red = mtDNA-encoded subunits",
             fontsize=11)
fig.subplots_adjust(left=0.24, top=0.95)
fig.savefig(f"{FIG}/tca_oxphos_gene_heatmap.png", dpi=170, bbox_inches="tight"); plt.close(fig)

outdf = pd.DataFrame(K, index=ordered, columns=[f"K27M_bin{b}" for b in range(NB)]).join(
        pd.DataFrame(O, index=ordered, columns=[f"KO_bin{b}" for b in range(NB)]))
outdf["group"] = [g for g,a,b in gseen for _ in range(b-a)]; outdf["pb_log2FC"] = lfc
outdf.round(3).to_csv(f"{OUT}/tca_oxphos_gene_heatmap_values.csv")
print(f"\n{len(ordered)} genes ({tca_rows} TCA + {len(ordered)-tca_rows} OxPhos)")
print("\ngenotype lean by group (mean pooled-z K27M-KO):")
for name,a,b in gseen:
    d = float(np.nanmean(K[a:b]-O[a:b]))
    print(f"  {name:<34} n={b-a:<2} {d:+.2f} ({'K27M' if d>0 else 'KO'})")
