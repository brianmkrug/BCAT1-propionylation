"""Figures for HSJ019 KO with r2 excluded (r1 + r3 only), plus a side-by-side
with the K27M integration on matched panels."""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import re, os
import numpy as np, pandas as pd, scanpy as sc
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

W = REPO_ROOT + "/scXenograft/work"
FIG = REPO_ROOT + "/scXenograft/figures"
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
      "PTPRZ1","MKI67","BCAT1","DLL3","SOX2"]
Afull = sc.read_h5ad(f"{W}/hsj019_norm.h5ad")

e = pd.read_csv(f"{W}/hsj019_ko13_embed.csv.gz", index_col=0)
s = sig_all.loc[e.index]
z = (s - s.mean()) / s.std().replace(0, np.nan)
for c in CLS:
    e[c] = z[[n for n in s.columns if cls(n) == c]].mean(1)
e["PROLIF_sig"] = z[[n for n in s.columns if n.endswith("-P")]].mean(1)
sub = Afull[e.index]
for g in MK:
    if g in sub.var_names:
        e[g] = np.asarray(sub[:, g].X.todense()).ravel()
e["label"] = e[CLS[:7]].idxmax(1)
e.to_csv(f"{W}/hsj019_ko13_annotated.csv.gz")

def scat(ax, x, y, c, s=4, **kw):
    ax.scatter(x, y, c=c, s=s, linewidths=0, rasterized=True, **kw)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_visible(False)

# ---- Fig 1: uncorrected vs Harmony by library, + label ----
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
libs = sorted(e["sample"].unique())
cmap = dict(zip(libs, [plt.cm.Dark2(0.0), plt.cm.Dark2(0.6)]))
for ax, (key, kt) in zip(axes[:2], [("uncorr", "uncorrected"),
                                    ("harmony", "Harmony on library")]):
    for l in libs:
        d = e[e["sample"] == l]
        scat(ax, d[f"UMAP1_{key}"], d[f"UMAP2_{key}"], [cmap[l]],
             label=f"{l.split('_',1)[1]} (n={len(d)})")
    ax.set_title(f"HSJ019 K27M-KO, r2 excluded\n{kt}", fontsize=10)
    ax.legend(markerscale=4, frameon=False, fontsize=8)
lc = dict(zip(sorted(e.label.unique()), plt.cm.Set2(np.linspace(0, 1, e.label.nunique()))))
for l, d in e.groupby("label"):
    scat(axes[2], d.UMAP1_harmony, d.UMAP2_harmony, [lc[l]], label=l)
axes[2].set_title("Harmony UMAP — dominant lineage meta-score", fontsize=10)
axes[2].legend(markerscale=4, frameon=False, fontsize=8)
fig.tight_layout(); fig.savefig(f"{FIG}/hsj019_ko_r1r3_harmony_umaps.png", dpi=170)
plt.close(fig)

# ---- Fig 2: meta-scores + markers ----
panels = ["RGC","GLIP","OPC","OL","ASTRO","PROLIF_sig",
          "PDGFRA","OLIG1","BCAS1","MBP","APOE","BCAT1"]
fig, axes = plt.subplots(3, 4, figsize=(16, 11))
for ax, p in zip(axes.ravel(), panels):
    v = e[p].values; lo, hi = np.percentile(v, [2, 98]); o = np.argsort(v)
    sm = ax.scatter(e.UMAP1_harmony.values[o], e.UMAP2_harmony.values[o],
                    c=np.clip(v[o], lo, hi), s=4, cmap="RdYlBu_r",
                    linewidths=0, rasterized=True)
    ax.set_title(p if p in CLS + ["PROLIF_sig"] else f"{p}  (lognorm)", fontsize=10)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_visible(False)
    fig.colorbar(sm, ax=ax, fraction=0.045)
fig.suptitle("HSJ019 K27M-KO (r1 + r3) — lineage meta-scores and markers, Harmony UMAP",
             fontsize=13)
fig.tight_layout(); fig.savefig(f"{FIG}/hsj019_ko_r1r3_metascores_harmony.png", dpi=170)
plt.close(fig)

# ---- mixing quantification vs the 3-library version ----
print("Library mixing (0 = perfectly mixed, 1 = library-pure)")
for tag, f in [("KO r1+r2+r3", "hsj019_ko_embed.csv.gz"),
               ("KO r1+r3   ", "hsj019_ko13_embed.csv.gz"),
               ("K27M r1+r2 ", "hsj019_k27m_embed.csv.gz")]:
    d = pd.read_csv(f"{W}/{f}", index_col=0)
    exp = d["sample"].value_counts(normalize=True)
    for key in ["uncorr", "harmony"]:
        ct = pd.crosstab(d[f"leiden_{key}"], d["sample"], normalize="index")
        dev = (ct - exp.reindex(ct.columns)).abs().sum(1).div(2)
        print(f"  {tag} {key:8s}: median {dev.median():.2f}, "
              f"clusters >0.75 pure {(dev > 0.75).sum()}/{len(dev)}")

print("\nKO r1+r3: label composition per library (%)")
print((pd.crosstab(e.label, e["sample"], normalize="columns")*100).round(1).to_string())
print("\nKO r1+r3: mean meta-score per library")
print(e.groupby("sample")[CLS[:5] + ["PROLIF_sig", "BCAT1", "APOE", "PDGFRA",
                                     "OLIG1"]].mean().round(2).to_string())
