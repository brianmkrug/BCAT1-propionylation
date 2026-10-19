"""Metabolic ssGSEA along pseudotime, and K27M vs KO at the library level.

ssGSEA ranks genes within each cell against the same gene universe (23,065
genes, both subsets), so raw ES are comparable across subsets and libraries.
"""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import numpy as np, pandas as pd
from scipy.stats import spearmanr
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

W = REPO_ROOT + "/scXenograft/work"
FIG = REPO_ROOT + "/scXenograft/figures"
SETS = [("k27m", "HSJ019 K27M"), ("ko13", "HSJ019 K27M-KO (r1+r3)")]
KEY = ["CHOL_Hallmark_cholesterol_homeostasis", "CHOL_Reactome_cholesterol_biosynthesis",
       "CHOL_WP_mevalonate_pathway", "CHOL_KEGGc_terpenoid_backbone",
       "FA_Hallmark_fatty_acid_metabolism", "FA_Reactome_mito_beta_oxidation",
       "FA_KEGGc_fatty_acid_degradation", "FA_KEGGc_unsaturated_FA_biosynthesis",
       "CTX_KEGGc_BCAA_degradation", "CTX_KEGGc_propanoate_metabolism"]

DAT = {}
for w, _ in SETS:
    d = pd.read_csv(f"{W}/hsj019_{w}_dpt.csv.gz", index_col=0)
    r = pd.read_csv(f"{W}/hsj019_{w}_ssgsea_raw.csv.gz", index_col=0).loc[d.index]
    z = pd.read_csv(f"{W}/hsj019_{w}_ssgsea_z.csv.gz", index_col=0).loc[d.index]
    # composite scores: mean z over member sets (robust to per-set dropout)
    z["CHOL_composite"] = z[[c for c in z.columns if c.startswith("CHOL_")]].mean(1)
    z["FA_composite"] = z[[c for c in z.columns if c.startswith("FA_")]].mean(1)
    r["CHOL_composite"] = r[[c for c in r.columns if c.startswith("CHOL_")]].mean(1)
    r["FA_composite"] = r[[c for c in r.columns if c.startswith("FA_")]].mean(1)
    DAT[w] = (d, r, z)

# ---------- correlation with pseudotime and lineage ----------
print("Spearman rho with DPT pseudotime (and with RGC / OL / BCAT1):\n")
rows = []
for w, ttl in SETS:
    d, r, z = DAT[w]
    print(f"--- {ttl} ---")
    for k in KEY + ["CHOL_composite", "FA_composite"]:
        v = z[k].values
        ok = np.isfinite(v)
        rp = spearmanr(d.dpt_pseudotime.values[ok], v[ok]).statistic
        rr = spearmanr(d.RGC.values[ok], v[ok]).statistic
        ro = spearmanr(d.OL.values[ok], v[ok]).statistic
        rb = spearmanr(d.BCAT1.values[ok], v[ok]).statistic
        print(f"  {k[:42]:<42} DPT {rp:+.3f}   RGC {rr:+.3f}   OL {ro:+.3f}   BCAT1 {rb:+.3f}")
        rows.append(dict(subset=w, signature=k, rho_dpt=rp, rho_RGC=rr,
                         rho_OL=ro, rho_BCAT1=rb))
    print()
pd.DataFrame(rows).round(4).to_csv(f"{W}/ssgsea_pseudotime_correlations.csv", index=False)

# ---------- K27M vs KO at the library level (raw ES, comparable universe) ----------
print("K27M vs KO, per-library mean raw ssGSEA ES "
      "(separation = every KO library differs from every K27M library):\n")
lib = []
for w, _ in SETS:
    d, r, _ = DAT[w]
    g = r.groupby(d["sample"]).mean()
    g["genotype"] = "K27M" if w == "k27m" else "K27M-KO"
    lib.append(g)
lib = pd.concat(lib)
lib.round(1).to_csv(f"{W}/ssgsea_by_library.csv")
k = lib[lib.genotype == "K27M"]
o = lib[lib.genotype == "K27M-KO"]
for c in KEY + ["CHOL_composite", "FA_composite"]:
    kv, ov = k[c].values, o[c].values
    sep = ("K27M > KO" if kv.min() > ov.max() else
           "KO > K27M" if ov.min() > kv.max() else "overlaps")
    flag = "  <<<" if sep != "overlaps" else ""
    print(f"  {c[:42]:<42} K27M {np.round(kv,0)}  KO {np.round(ov,0)}  -> {sep}{flag}")

# ---------- Fig: composites + key sets along pseudotime ----------
fig, axes = plt.subplots(2, 3, figsize=(16, 9))
bins = np.linspace(0, 1, 26); ctr = (bins[:-1] + bins[1:]) / 2
for row, (w, ttl) in enumerate(SETS):
    d, r, z = DAT[w]
    idx = np.digitize(d.dpt_pseudotime, bins) - 1
    ax = axes[row, 0]
    for c, col in [("CHOL_composite", "#C1272D"), ("FA_composite", "#2B6CB0")]:
        m = [np.nanmean(z[c][idx == i]) if (idx == i).sum() > 15 else np.nan for i in range(25)]
        ax.plot(ctr, m, lw=2.4, color=col, label=c.replace("_composite", ""))
    for c, col, ls in [("CTX_KEGGc_BCAA_degradation", "#666666", "--"),
                       ("CTX_KEGGc_propanoate_metabolism", "#999999", ":")]:
        m = [np.nanmean(z[c][idx == i]) if (idx == i).sum() > 15 else np.nan for i in range(25)]
        ax.plot(ctr, m, lw=1.6, color=col, ls=ls, label=c.replace("CTX_KEGGc_", ""))
    ax.axhline(0, lw=.6, color="k"); ax.legend(fontsize=7, frameon=False)
    ax.set_xlabel("DPT pseudotime"); ax.set_ylabel("mean z-scored ssGSEA")
    ax.set_title(f"{ttl}\ncomposite metabolic scores", fontsize=10)

    for col_i, (comp, cm) in enumerate([("CHOL_composite", "RdBu_r"),
                                        ("FA_composite", "RdBu_r")]):
        ax = axes[row, col_i + 1]
        v = z[comp].values
        bad = ~np.isfinite(v)
        o_ = np.argsort(np.where(bad, -np.inf, v))[bad.sum():]
        sm = ax.scatter(d.UMAP1_harmony.values[o_], d.UMAP2_harmony.values[o_],
                        c=np.clip(v[o_], -2, 2), s=3, cmap=cm, vmin=-2, vmax=2,
                        linewidths=0, rasterized=True)
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values(): sp.set_visible(False)
        ax.set_title(f"{comp.replace('_composite','')} composite (z)", fontsize=10)
        fig.colorbar(sm, ax=ax, fraction=.045)
    for sp in ["top", "right"]: axes[row, 0].spines[sp].set_visible(False)
fig.suptitle("Cholesterol/mevalonate and fatty acid ssGSEA — UMAP and pseudotime",
             fontsize=13)
fig.tight_layout(); fig.savefig(f"{FIG}/hsj019_ssgsea_composites.png", dpi=170)
plt.close(fig)
print("\nfigures written")
