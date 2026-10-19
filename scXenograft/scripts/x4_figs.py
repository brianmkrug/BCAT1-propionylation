"""Figures for the pathway-variability analysis."""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

W = REPO_ROOT + "/scXenograft/work"
FIG = REPO_ROOT + "/scXenograft/figures"
NBIN = 20
SETS = [("k27m", "HSJ019 K27M"), ("ko13", "HSJ019 K27M-KO (r1+r3)")]

def short(s):
    return (s.replace("KEGG_", "").replace("HM_", "Hallmark ")
             .replace("_", " ")[:40])

CC = {w: pd.read_csv(f"{W}/hsj019_{w}_pathway_cellcycle_control.csv", index_col=0)
      for w, _ in SETS}

# ---- Fig 1: eta2 before vs after removing cycling cells ----
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
for ax, (w, ttl) in zip(axes, SETS):
    r = CC[w]
    keep = r.z_null > 3
    ax.scatter(r.eta2[~keep], r.eta2_G1only[~keep], s=18, c="#cccccc",
               label="not above random null")
    ax.scatter(r.eta2[keep], r.eta2_G1only[keep], s=34, c="#C1272D",
               label="above null (z>3)")
    lim = max(r.eta2.max(), r.eta2_G1only.max())*1.08
    ax.plot([0, lim], [0, lim], lw=.8, ls="--", c="k")
    for p in r[keep].index:
        ax.annotate(short(p), (r.eta2[p], r.eta2_G1only[p]), fontsize=6.5,
                    xytext=(3, 2), textcoords="offset points")
    for p in r.sort_values("eta2_G1only", ascending=False).head(8).index:
        if p not in r[keep].index:
            ax.annotate(short(p), (r.eta2[p], r.eta2_G1only[p]), fontsize=6.5,
                        color="#2B6CB0", xytext=(3, 2), textcoords="offset points")
    ax.set_xlabel("eta2 vs pseudotime (all cells)")
    ax.set_ylabel("eta2 vs pseudotime (G1 cells only)")
    ax.set_title(f"{ttl}\nbelow the diagonal = driven by the cell cycle", fontsize=10)
    ax.legend(fontsize=7, frameon=False)
    for sp in ["top", "right"]: ax.spines[sp].set_visible(False)
fig.tight_layout(); fig.savefig(f"{FIG}/hsj019_pathway_eta2_cellcycle.png", dpi=170)
plt.close(fig)

# ---- Fig 2: trajectory heatmaps of top pathways (G1-controlled ranking) ----
for w, ttl in SETS:
    d = pd.read_csv(f"{W}/hsj019_{w}_dpt.csv.gz", index_col=0)
    S = pd.read_csv(f"{W}/hsj019_{w}_allmetab_resid.csv.gz", index_col=0).loc[d.index]
    r = CC[w]
    top = r[(r.lib_consistency > 0.5)].sort_values("eta2_G1only", ascending=False).head(25).index
    pt_bin = pd.qcut(d.dpt_pseudotime, NBIN, labels=False, duplicates="drop")
    P = np.vstack([[np.nanmean(S[p][pt_bin == b]) for b in range(NBIN)] for p in top])
    P = (P - P.mean(1, keepdims=True)) / P.std(1, keepdims=True)
    ordr = np.argsort(np.argmax(P, 1))
    fig, ax = plt.subplots(figsize=(9, 8))
    im = ax.imshow(P[ordr], aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels([short(top[i]) for i in ordr], fontsize=7.5)
    ax.set_xticks([0, NBIN//2, NBIN-1]); ax.set_xticklabels(["early", "mid", "late"])
    ax.set_xlabel("DPT pseudotime bin")
    ax.set_title(f"{ttl}\ntop 25 metabolic pathways by G1-only trajectory variation\n"
                 f"(depth-residualised, row-scaled)", fontsize=10)
    fig.colorbar(im, ax=ax, fraction=.03, label="z (within pathway)")
    fig.tight_layout(); fig.savefig(f"{FIG}/hsj019_{w}_pathway_trajectory_heatmap.png", dpi=170)
    plt.close(fig)

# ---- Fig 3: K27M vs KO ----
k, o = CC["k27m"], CC["ko13"]
common = k.index.intersection(o.index)
fig, axes = plt.subplots(1, 2, figsize=(14, 6.2))
ax = axes[0]
ax.scatter(k.eta2_G1only[common], o.eta2_G1only[common], s=26, c="#444444")
lim = max(k.eta2_G1only[common].max(), o.eta2_G1only[common].max())*1.1
ax.plot([0, lim], [0, lim], lw=.8, ls="--", c="k")
lab = set(k.eta2_G1only[common].nlargest(8).index) | set(o.eta2_G1only[common].nlargest(8).index)
for p in lab:
    ax.annotate(short(p), (k.eta2_G1only[p], o.eta2_G1only[p]), fontsize=6.5,
                xytext=(3, 2), textcoords="offset points")
ax.set_xlabel("K27M — eta2 (G1 only)"); ax.set_ylabel("KO — eta2 (G1 only)")
ax.set_title("Trajectory variation, cell-cycle controlled", fontsize=10)

ax = axes[1]
delta = (o.eta2_G1only[common] - k.eta2_G1only[common]).sort_values()
sel = pd.concat([delta.head(12), delta.tail(12)])
ax.barh(range(len(sel)), sel.values,
        color=["#C1272D" if v < 0 else "#2B6CB0" for v in sel.values])
ax.set_yticks(range(len(sel)))
ax.set_yticklabels([short(p) for p in sel.index], fontsize=7)
ax.axvline(0, lw=.8, c="k")
ax.set_xlabel("eta2(KO) - eta2(K27M),  G1 cells only")
ax.set_title("red = varies more along K27M trajectory\nblue = varies more along KO",
             fontsize=10)
for a in axes:
    for sp in ["top", "right"]: a.spines[sp].set_visible(False)
fig.tight_layout(); fig.savefig(f"{FIG}/hsj019_pathway_variability_genotype.png", dpi=170)
plt.close(fig)
print("figures written")

# ---- export the ranked table ----
out = pd.DataFrame({
    "eta2_K27M": k.eta2.reindex(common), "eta2_G1_K27M": k.eta2_G1only.reindex(common),
    "rho_K27M": k.rho.reindex(common), "rho_partial_cc_K27M": k.rho_partial_cc.reindex(common),
    "libcons_K27M": k.lib_consistency.reindex(common), "z_null_K27M": k.z_null.reindex(common),
    "eta2_KO": o.eta2.reindex(common), "eta2_G1_KO": o.eta2_G1only.reindex(common),
    "rho_KO": o.rho.reindex(common), "rho_partial_cc_KO": o.rho_partial_cc.reindex(common),
    "libcons_KO": o.lib_consistency.reindex(common), "z_null_KO": o.z_null.reindex(common),
})
out["max_eta2_G1"] = out[["eta2_G1_K27M", "eta2_G1_KO"]].max(1)
out.sort_values("max_eta2_G1", ascending=False).round(4).to_csv(
    REPO_ROOT + "/scXenograft/pathway_trajectory_variability.csv")
print("\nTop 15 by cell-cycle-controlled trajectory variation (either genotype):")
print(out.sort_values("max_eta2_G1", ascending=False).head(15)[
    ["eta2_G1_K27M", "rho_K27M", "eta2_G1_KO", "rho_KO"]].round(3).to_string())
