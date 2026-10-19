"""Nucleotide / one-carbon-folate pathways along pseudotime.

These were the top per-genotype trajectory hits earlier and are shared by both
genotypes -- but were shown to be largely cell-cycle-driven. Plotted here:
  (1) pooled shared-axis genotype profiles (analogous to the metabolic profiles fig)
  (2) all-cells vs G1-only, per genotype, to show the cell-cycle dependence
  (3) score vs proliferation, to make the driver explicit
"""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

W = REPO_ROOT + "/scXenograft/work"
OUT = REPO_ROOT + "/scXenograft"; FIG = f"{OUT}/figures"
GT = {"K27M": "#C1272D", "K27M-KO": "#2B6CB0"}
NB = 20
PATH = ["KEGG_Purine metabolism", "KEGG_Pyrimidine metabolism",
        "KEGG_One carbon pool by folate", "KEGG_Folate biosynthesis"]
def short(s): return s.replace("KEGG_", "")

# ---------- pooled shared axis ----------
emb = pd.read_csv(f"{W}/pooled_dpt.csv.gz", index_col=0)
emb["bin"] = pd.qcut(emb.dpt, NB, labels=False, duplicates="drop")
D = pd.read_csv(f"{W}/pooled_metab_resid.csv.gz", index_col=0).loc[emb.index]
rep = pd.read_csv(f"{W}/pooled_bin_representation.csv", index_col=0)
BA = pd.read_csv(f"{OUT}/pooled_pathway_genotype_binconditional.csv", index_col=0)
mids = rep.dpt_mid.reindex(range(NB)).values

fig, axes = plt.subplots(1, 4, figsize=(18, 4.4))
for ax, c in zip(axes, PATH):
    v = D[c].values
    for g in ["K27M", "K27M-KO"]:
        prof = [np.nanmean(v[(emb.genotype == g).values & (emb.bin == b).values]) for b in range(NB)]
        ax.plot(mids, prof, lw=2.4, color=GT[g], label=g)
    # proliferation reference on a twin axis
    ax2 = ax.twinx()
    pr = [np.nanmean(emb.PROLIF_sig.values[(emb.bin == b).values]) for b in range(NB)]
    ax2.plot(mids, pr, lw=1.4, color="#888888", ls=":", label="proliferation (both gt)")
    ax2.set_yticks([]); ax2.set_ylabel("prolif score", color="#888888", fontsize=7)
    ax.set_title(f"{short(c)}\nbinadj Δ(K27M-KO)={BA.binadj_delta[c]:+.3f} "
                 f"({'ns/common' if abs(BA.binadj_delta[c])<0.012 else 'differs'})", fontsize=9)
    ax.set_xlabel("shared pseudotime"); ax.axhline(0, lw=.5, c="k")
    for sp in ["top"]: ax.spines[sp].set_visible(False)
axes[0].legend(frameon=False, fontsize=8, loc="upper right")
fig.suptitle("Nucleotide / one-carbon-folate metabolism along the SHARED pseudotime axis "
             "(early-enriched, common to both genotypes; grey dotted = proliferation)", fontsize=12)
fig.tight_layout(); fig.savefig(f"{FIG}/pooled_nucleotide_1C_profiles.png", dpi=160); plt.close(fig)

# ---------- per-genotype separate axes: all cells vs G1 only ----------
fig, axes = plt.subplots(2, 4, figsize=(18, 8.4))
for row, (w, gt, ttl) in enumerate([("k27m", "K27M", "HSJ019 K27M"),
                                    ("ko13", "K27M-KO", "HSJ019 KO (r1+r3)")]):
    d = pd.read_csv(f"{W}/hsj019_{w}_dpt.csv.gz", index_col=0)
    S = pd.read_csv(f"{W}/hsj019_{w}_allmetab_resid.csv.gz", index_col=0).loc[d.index]
    b_all = pd.qcut(d.dpt_pseudotime, NB, labels=False, duplicates="drop")
    g1 = (d.phase == "G1").values
    b_g1 = pd.qcut(d.dpt_pseudotime[g1], NB, labels=False, duplicates="drop")
    mid_all = [d.dpt_pseudotime[b_all == b].mean() for b in range(NB)]
    mid_g1 = [d.dpt_pseudotime[g1][b_g1 == b].mean() for b in range(NB)]
    for ax, c in zip(axes[row], PATH):
        v = S[c].values
        pa = [np.nanmean(v[b_all == b]) for b in range(NB)]
        pg = [np.nanmean(v[g1][b_g1 == b]) for b in range(NB)]
        ax.plot(mid_all, pa, lw=2.2, color=GT[gt], label="all cells")
        ax.plot(mid_g1, pg, lw=2.0, color=GT[gt], ls="--", label="G1 only")
        ax.axhline(0, lw=.5, c="k")
        ax.set_title(f"{ttl}\n{short(c)}", fontsize=9)
        ax.set_xlabel("pseudotime (own axis)")
        for sp in ["top", "right"]: ax.spines[sp].set_visible(False)
    axes[row][0].legend(frameon=False, fontsize=8)
fig.suptitle("Same pathways on each genotype's own pseudotime axis: solid = all cells, "
             "dashed = G1 only\n(the early peak flattens once cycling cells are removed "
             "-> largely cell-cycle-driven)", fontsize=12)
fig.tight_layout(); fig.savefig(f"{FIG}/nucleotide_1C_allcells_vs_G1.png", dpi=160); plt.close(fig)

# ---------- numbers ----------
print("Pooled shared axis: genotype difference at matched pseudotime (bin-adjusted)")
for c in PATH:
    print(f"  {short(c):<28} Δ(K27M-KO) {BA.binadj_delta[c]:+.4f}  "
          f"naive {BA.naive_delta[c]:+.4f}  shape_r {BA.shape_r[c]:+.3f}  "
          f"lib_consist {BA.lib_consist[c]:.2f}")
print("\nFraction of trajectory variance lost when cycling cells removed (from earlier):")
for w, gt in [("k27m", "K27M"), ("ko13", "KO")]:
    cc = pd.read_csv(f"{W}/hsj019_{w}_pathway_cellcycle_control.csv", index_col=0)
    print(f"  {gt}:")
    for c in PATH:
        if c in cc.index:
            r = cc.eta2_G1only[c]/cc.eta2[c] if cc.eta2[c] else np.nan
            print(f"    {short(c):<28} eta2 all {cc.eta2[c]:.3f} -> G1 {cc.eta2_G1only[c]:.3f} "
                  f"(retains {100*r:.0f}%)")
