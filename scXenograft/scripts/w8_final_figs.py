"""Final metabolic figures using the depth-robust control-bin scores."""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

W = REPO_ROOT + "/scXenograft/work"
FIG = REPO_ROOT + "/scXenograft/figures"
SETS = [("k27m", "HSJ019 K27M"), ("ko13", "HSJ019 K27M-KO (r1+r3)")]
KEY = ["CHOL_composite", "CHOL_Hallmark_cholesterol_homeostasis",
       "CHOL_Reactome_cholesterol_biosynthesis", "CHOL_WP_mevalonate_pathway",
       "FA_composite", "FA_Hallmark_fatty_acid_metabolism",
       "FA_KEGGc_unsaturated_FA_biosynthesis", "CTX_KEGGc_BCAA_degradation"]

for w, ttl in SETS:
    d = pd.read_csv(f"{W}/hsj019_{w}_dpt.csv.gz", index_col=0)
    m = pd.read_csv(f"{W}/hsj019_{w}_metab_ctrlscore.csv.gz", index_col=0).loc[d.index]
    fig, axes = plt.subplots(2, 4, figsize=(17, 8.5))
    for ax, k in zip(axes.ravel(), KEY):
        v = m[k].values
        lo, hi = np.percentile(v, [2, 98])
        o = np.argsort(v)
        sm = ax.scatter(d.UMAP1_harmony.values[o], d.UMAP2_harmony.values[o],
                        c=np.clip(v[o], lo, hi), s=3, cmap="RdBu_r",
                        vmin=-max(abs(lo), abs(hi)), vmax=max(abs(lo), abs(hi)),
                        linewidths=0, rasterized=True)
        ax.set_title(k.replace("_", " ")[:38], fontsize=9)
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values(): sp.set_visible(False)
        fig.colorbar(sm, ax=ax, fraction=.045)
    fig.suptitle(f"{ttl} — metabolic panel, control-bin scores (depth-robust), "
                 f"colour = score centred on 0", fontsize=12)
    fig.tight_layout(); fig.savefig(f"{FIG}/hsj019_{w}_metabolic_ctrlscore_umaps.png", dpi=170)
    plt.close(fig)

# trends along pseudotime, per library
fig, axes = plt.subplots(1, 2, figsize=(13, 4.6))
bins = np.linspace(0, 1, 26); ctr = (bins[:-1]+bins[1:])/2
for ax, (w, ttl) in zip(axes, SETS):
    d = pd.read_csv(f"{W}/hsj019_{w}_dpt.csv.gz", index_col=0)
    m = pd.read_csv(f"{W}/hsj019_{w}_metab_ctrlscore.csv.gz", index_col=0).loc[d.index]
    libs = sorted(d["sample"].unique())
    for c, col in [("CHOL_composite", "#C1272D"), ("FA_composite", "#2B6CB0")]:
        for l, ls in zip(libs, ["-", "--"]):
            s = d["sample"] == l
            ii = np.digitize(d.dpt_pseudotime[s], bins) - 1
            vals = m[c][s].values
            y = [np.nanmean(vals[ii == i]) if (ii == i).sum() > 15 else np.nan
                 for i in range(25)]
            ax.plot(ctr, y, ls=ls, lw=2, color=col,
                    label=f"{c.replace('_composite','')} — {l.split('_',1)[1]}")
    ax.axhline(0, lw=.6, color="k")
    ax.set_xlabel("DPT pseudotime"); ax.set_ylabel("control-bin score")
    ax.set_title(ttl, fontsize=10); ax.legend(fontsize=7, frameon=False)
    for sp in ["top", "right"]: ax.spines[sp].set_visible(False)
fig.suptitle("Cholesterol and fatty acid programs along pseudotime, per library "
             "(depth-robust scoring)", fontsize=12)
fig.tight_layout(); fig.savefig(f"{FIG}/hsj019_metabolic_pseudotime_ctrlscore.png", dpi=170)
plt.close(fig)

# library-level genotype comparison on control-bin scores
lib = []
for w, _ in SETS:
    d = pd.read_csv(f"{W}/hsj019_{w}_dpt.csv.gz", index_col=0)
    m = pd.read_csv(f"{W}/hsj019_{w}_metab_ctrlscore.csv.gz", index_col=0).loc[d.index]
    g = m.groupby(d["sample"]).mean()
    g["genotype"] = "K27M" if w == "k27m" else "K27M-KO"
    lib.append(g)
lib = pd.concat(lib)
lib.round(4).to_csv(f"{W}/metab_ctrlscore_by_library.csv")
k = lib[lib.genotype == "K27M"]; o = lib[lib.genotype == "K27M-KO"]
print("K27M vs KO, per-library mean control-bin score:\n")
for c in KEY:
    kv, ov = k[c].values, o[c].values
    sep = ("K27M > KO" if kv.min() > ov.max() else
           "KO > K27M" if ov.min() > kv.max() else "overlaps")
    print(f"  {c[:42]:<42} K27M {np.round(kv,3)}  KO {np.round(ov,3)}  -> {sep}"
          f"{'  <<<' if sep != 'overlaps' else ''}")
print("\nfigures written")
