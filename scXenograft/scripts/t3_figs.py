"""Figures for the within-genotype Harmony integrations of HSJ019.
Meta-scores are re-z-scored within each genotype subset so they describe
structure inside that genotype, not the (confounded) K27M-vs-KO contrast.
"""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import re, os
import numpy as np, pandas as pd, scanpy as sc
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

W = REPO_ROOT + "/scXenograft/work"
FIG = REPO_ROOT + "/scXenograft/figures"
os.makedirs(FIG, exist_ok=True)
CLS = ["RGC", "GLIP", "OPC", "OL", "ASTRO", "EPEN", "NEURON", "NON_NEURAL"]

def cls(nm):
    s = re.sub(r"^[FP]-[a-z0-9]+_", "", nm)
    if re.match(r"^(V|D|G)?RGC", s): return "RGC"
    if s.startswith("OPC"): return "OPC"
    if s.startswith(("NFOL", "MOL")): return "OL"
    if s.startswith(("ASTR", "ASTP", "OPAS")): return "ASTRO"
    if s.startswith("GLIP"): return "GLIP"
    if s.startswith("EPEN"): return "EPEN"
    if s.startswith(("CEXN","PEXN","NEUR","CINHN","PINH","MGINH","GABAN","EXIP",
                     "INIP","UBCN","SERN","ACHN","THLN","MFN","CJRZN","SMSN")): return "NEURON"
    if s.startswith(("PERI","ENDO","VSMC","VLM","MNG","MGL","MAC","CPLX","SCHW",
                     "OTIC","CHEM","ASEP","FLRPT")): return "NON_NEURAL"
    return "OTHER"

sig_all = pd.read_csv(f"{W}/hsj019_sigscores.csv.gz", index_col=0)
MK = ["PDGFRA","OLIG1","OLIG2","SOX10","MBP","PLP1","BCAS1","GFAP","APOE","FABP7",
      "PTPRZ1","MKI67","TOP2A","BCAT1","DLL3","SOX2"]
Afull = sc.read_h5ad(f"{W}/hsj019_norm.h5ad")

D = {}
for which, gt in [("k27m", "K27M"), ("ko", "K27M-KO")]:
    e = pd.read_csv(f"{W}/hsj019_{which}_embed.csv.gz", index_col=0)
    s = sig_all.loc[e.index]
    z = (s - s.mean()) / s.std().replace(0, np.nan)
    for c in CLS:
        cols = [n for n in s.columns if cls(n) == c]
        e[c] = z[cols].mean(1)
    e["PROLIF_sig"] = z[[n for n in s.columns if n.endswith("-P")]].mean(1)
    sub = Afull[e.index]
    for g in MK:
        if g in sub.var_names:
            e[g] = np.asarray(sub[:, g].X.todense()).ravel()
    e["label"] = e[CLS[:7]].idxmax(1)
    e.to_csv(f"{W}/hsj019_{which}_annotated.csv.gz")
    D[which] = e

def scat(ax, x, y, c, s=3, **kw):
    ax.scatter(x, y, c=c, s=s, linewidths=0, rasterized=True, **kw)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_visible(False)

# ---- Fig A: uncorrected vs within-genotype Harmony, coloured by library ----
fig, axes = plt.subplots(2, 3, figsize=(14, 9))
for row, (which, ttl) in enumerate([("k27m", "HSJ019 K27M (2 libraries)"),
                                    ("ko", "HSJ019 K27M-KO (3 libraries)")]):
    e = D[which]
    libs = sorted(e["sample"].unique())
    cmap = dict(zip(libs, plt.cm.Dark2(np.linspace(0, .8, len(libs)))))
    for col, (key, kt) in enumerate([("uncorr", "uncorrected"), ("harmony", "Harmony on library")]):
        ax = axes[row, col]
        for l in libs:
            sub = e[e["sample"] == l]
            scat(ax, sub[f"UMAP1_{key}"], sub[f"UMAP2_{key}"], [cmap[l]],
                 label=l.split("_", 1)[1])
        ax.set_title(f"{ttl}\n{kt}", fontsize=10)
        if col == 0: ax.legend(markerscale=5, frameon=False, fontsize=7)
    ax = axes[row, 2]
    lb = e.label
    lc = dict(zip(sorted(lb.unique()), plt.cm.Set2(np.linspace(0, 1, lb.nunique()))))
    for l, sub in e.groupby(lb):
        scat(ax, sub.UMAP1_harmony, sub.UMAP2_harmony, [lc[l]], label=l)
    ax.set_title("Harmony UMAP — dominant lineage meta-score", fontsize=10)
    ax.legend(markerscale=5, frameon=False, fontsize=7)
fig.suptitle("HSJ019 — Harmony applied separately within each genotype", fontsize=13)
fig.tight_layout(); fig.savefig(f"{FIG}/hsj019_bygenotype_harmony_umaps.png", dpi=170)
plt.close(fig)

# ---- Fig B/C: meta-scores + markers on the within-genotype Harmony UMAP ----
for which, ttl in [("k27m", "HSJ019 K27M"), ("ko", "HSJ019 K27M-KO")]:
    e = D[which]
    panels = ["RGC", "GLIP", "OPC", "OL", "ASTRO", "PROLIF_sig",
              "PDGFRA", "OLIG1", "BCAS1", "MBP", "APOE", "BCAT1"]
    fig, axes = plt.subplots(3, 4, figsize=(16, 11))
    for ax, p in zip(axes.ravel(), panels):
        v = e[p].values
        lo, hi = np.percentile(v, [2, 98])
        o = np.argsort(v)
        sm = ax.scatter(e.UMAP1_harmony.values[o], e.UMAP2_harmony.values[o],
                        c=np.clip(v[o], lo, hi), s=3, cmap="RdYlBu_r",
                        linewidths=0, rasterized=True)
        ax.set_title(p + ("" if p in CLS or p == "PROLIF_sig" else "  (lognorm)"), fontsize=10)
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values(): sp.set_visible(False)
        fig.colorbar(sm, ax=ax, fraction=0.045)
    fig.suptitle(f"{ttl} — lineage meta-scores (top 6) and markers, Harmony UMAP", fontsize=13)
    fig.tight_layout(); fig.savefig(f"{FIG}/hsj019_{which}_metascores_harmony.png", dpi=170)
    plt.close(fig)

# ---- library mixing quantification ----
print("Library mixing after within-genotype Harmony")
for which in ["k27m", "ko"]:
    e = D[which]
    exp = e["sample"].value_counts(normalize=True)
    for key in ["uncorr", "harmony"]:
        ct = pd.crosstab(e[f"leiden_{key}"], e["sample"], normalize="index")
        # mean abs deviation from expected library proportion, per cluster
        dev = (ct - exp.reindex(ct.columns)).abs().sum(1).div(2)
        print(f"  {which:5s} {key:8s}: median cluster imbalance "
              f"{dev.median():.2f} (0 = perfectly mixed, 1 = library-pure), "
              f"clusters >0.75 pure: {(dev > 0.75).sum()}/{len(dev)}")
print("\nfigures written to", FIG)
