#!/usr/bin/env python
"""
(1) Genomic-category composition of each mark -- where the signal actually sits,
    and specifically whether H3K27ac is more promoter-focused than the others.
(2) Pooled, spike-in-scaled promoter deltas, related to gene expression.

Pooling and scaling follow the settled convention: replicates SUMMED (all
available depth), CPM against the GENOME-WIDE pooled depth (not the promoter
total -- otherwise a mark that redistributes into promoters would look unchanged),
then the spike-in factor on the shBCAT1 side.

    ENL x0.80   H3K9pr x0.62   H3K27ac x1.00   H3K4me3 x1.00 (assumed)  IgG x1.00

Category enrichment = (% of fragments in category) / (% of genome in category),
so 1.0 means "as often as chance". Categories are disjoint, so these are shares
of one pie and cannot double-count a promoter that is also exonic.

Note on the spike-in and correlation: multiplying shBCAT1 by a constant shifts
every promoter's delta by log2(factor). That changes the SIGN of individual
promoters (whether they gain or lose) but not Spearman against expression, which
is shift-invariant. Both are reported.
"""
from __future__ import annotations
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
import pandas as pd
from scipy.stats import rankdata

ROOT = os.path.expanduser(REPO_ROOT + "/CUTRUN")
QC = os.path.join(ROOT, "analysis_QC_normalization", "multimark_2026-07")
PC = os.path.join(QC, "promcat")
ANN = os.path.join(QC, "annotation")
CORE = os.path.join(ROOT, "H3K4me3", "reference", "core.bed")
EXPR = os.path.join(ROOT, "H3K4me3", "expression", "expression_harmonised.csv")
PLOTS, DATA = os.path.join(QC, "plots"), os.path.join(QC, "data")

SPIKEIN = {"ENL": 0.80, "H3K9pr": 0.62, "H3K27ac": 1.00, "H3K4me3": 1.00, "IgG": 1.00}
ASSUMED = {"H3K4me3"}
POOLS = {
    "H3K4me3": {"shNT": ["shNT-H3K4me3-r1", "shNT-H3K4me3-r2"],
                "shBCAT1": ["shBCAT1-H3K4me3-r1", "shBCAT1-H3K4me3-r2"]},
    "H3K27ac": {"shNT": ["shNT-H3K27ac-r1", "shNT-H3K27ac-r2"],
                "shBCAT1": ["shBCAT1-H3K27ac-r1", "shBCAT1-H3K27ac-r2"]},
    "H3K9pr": {"shNT": ["shNT-H3K9pr-r1", "shNT-H3K9pr-r2"],
               "shBCAT1": ["shBCAT1-H3K9pr-r1", "shBCAT1-H3K9pr-r2_60M"]},
    "ENL": {"shNT": ["shNT-ENL-r1", "shNT-ENL-r2"],
            "shBCAT1": ["shBCAT1-ENL-r1", "shBCAT1-ENL-r2"]},
    "IgG": {"shNT": ["shNT-IgG-r1"], "shBCAT1": ["shBCAT1-IgG-r1"]},
}
MARKS = ["H3K4me3", "H3K27ac", "H3K9pr", "ENL", "IgG"]
CATS = ["promoter", "exonic", "intronic", "intergenic"]
CAT_COLOR = {"promoter": "#c0392b", "exonic": "#e8862a",
             "intronic": "#2166ac", "intergenic": "#95a5a6"}
MARK_COLOR = {"H3K4me3": "#2166ac", "H3K27ac": "#e8862a", "H3K9pr": "#6c3483",
              "ENL": "#1e8449", "IgG": "#808080"}


def sp(a, b):
    ra, rb = rankdata(a), rankdata(b)
    ra = ra - ra.mean(); rb = rb - rb.mean()
    return float(ra @ rb / (np.linalg.norm(ra) * np.linalg.norm(rb)))


def genome_fractions():
    p = pd.read_csv(os.path.join(ANN, "partition.bed"), sep="\t", header=None,
                    names=["c", "s", "e", "cat"])
    bp = p.assign(w=p.e - p.s).groupby("cat").w.sum()
    tot = bp.sum()
    return {c: bp[i + 1] / tot for i, c in enumerate(CATS)}, tot


def load_categories():
    rows = {}
    for m in MARKS:
        for cond in ("shNT", "shBCAT1"):
            agg = {c: 0 for c in CATS}
            tot = 0
            for lab in POOLS[m][cond]:
                d = pd.read_csv(os.path.join(PC, f"{lab}.category.tsv"), sep="\t",
                                index_col=0).fragments
                for c in CATS:
                    agg[c] += int(d[c])
                tot += int(d["total"])
            rows[(m, cond)] = (agg, tot)
    return rows


def load_promoter_matrix(n_prom):
    M = {}
    for m in MARKS:
        for cond in ("shNT", "shBCAT1"):
            v = np.zeros(n_prom, dtype=np.int64)
            for lab in POOLS[m][cond]:
                d = pd.read_csv(os.path.join(PC, f"{lab}.promoter.tsv.gz"), sep="\t")
                v[d.promoter_row.to_numpy() - 1] += d.fragments.to_numpy()
            M[(m, cond)] = v
    return M


# ---------------------------------------------------------------- figure 1 ----
def per_library_enrichment(gfrac):
    """Promoter/intronic enrichment for every library separately.

    The pooled numbers are depth-weighted, and for H3K27ac the r1 batch supplies
    ~90% of the depth -- so a pooled category profile is essentially the r1 answer.
    Splitting by library is the only way to see that.
    """
    rows = []
    for m in MARKS:
        for cond in ("shNT", "shBCAT1"):
            for lab in POOLS[m][cond]:
                d = pd.read_csv(os.path.join(PC, f"{lab}.category.tsv"), sep="\t",
                                index_col=0).fragments
                tot = int(d["total"])
                rep = "r2" if "-r2" in lab else "r1"
                rows.append({"label": lab, "mark": m, "cond": cond, "rep": rep,
                             "total": tot,
                             **{f"{c}_x": int(d[c]) / tot / gfrac[c] for c in CATS}})
    return pd.DataFrame(rows)


def fig_categories(cats, gfrac, per_lib, path):
    fig, axes = plt.subplots(1, 4, figsize=(20.2, 4.9),
                             gridspec_kw={"width_ratios": [1.15, 1.15, 1.0, 1.25]})

    # (a) composition, shNT
    ax = axes[0]
    bottoms = np.zeros(len(MARKS))
    for c in CATS:
        vals = np.array([cats[(m, "shNT")][0][c] / cats[(m, "shNT")][1] * 100
                         for m in MARKS])
        ax.bar(MARKS, vals, bottom=bottoms, color=CAT_COLOR[c], label=c,
               edgecolor="white", linewidth=0.6)
        bottoms += vals
    for i, c in enumerate(CATS):
        ax.axhline(sum(gfrac[CATS[j]] for j in range(i + 1)) * 100,
                   color="k", ls=":", lw=1.0, zorder=5)
    ax.set_ylabel("% of fragments")
    ax.set_title("Signal composition (shNT)\ndotted = genome composition "
                 "(what chance would give)", fontsize=9.5)
    ax.legend(fontsize=7.4, frameon=False, ncol=4, loc="upper center",
              bbox_to_anchor=(0.5, -0.10))
    ax.set_ylim(0, 100)

    # (b) enrichment over genome fraction
    ax = axes[1]
    w = 0.2
    x = np.arange(len(MARKS))
    for k, c in enumerate(CATS):
        vals = [cats[(m, "shNT")][0][c] / cats[(m, "shNT")][1] / gfrac[c] for m in MARKS]
        ax.bar(x + (k - 1.5) * w, vals, w, color=CAT_COLOR[c], label=c)
    ax.axhline(1, color="k", lw=1.0, ls="--")
    ax.set_xticks(x)
    ax.set_xticklabels(MARKS, fontsize=9)
    ax.set_ylabel("enrichment vs genome (x)")
    ax.set_title("Category enrichment (shNT)\n1.0 = as often as chance", fontsize=9.5)
    ax.legend(fontsize=7.4, frameon=False, ncol=2)

    # (c) promoter enrichment, the headline comparison
    ax = axes[2]
    for i, m in enumerate(MARKS):
        for j, cond in enumerate(("shNT", "shBCAT1")):
            e = cats[(m, cond)][0]["promoter"] / cats[(m, cond)][1] / gfrac["promoter"]
            ax.bar(i + (j - 0.5) * 0.38, e, 0.38, color=MARK_COLOR[m],
                   alpha=1.0 if cond == "shNT" else 0.55,
                   edgecolor="k", linewidth=0.5)
    ax.axhline(1, color="k", lw=1.0, ls="--")
    ax.set_xticks(range(len(MARKS)))
    ax.set_xticklabels(MARKS, fontsize=9)
    ax.set_ylabel("promoter enrichment (x)")
    ax.set_title("Promoter focus per mark, POOLED\nsolid = shNT, pale = shBCAT1",
                 fontsize=9.5)

    # (d) the same thing per library — the pooled panel hides a 2x batch split
    ax = axes[3]
    for i, m in enumerate(MARKS):
        sub = per_lib[per_lib["mark"] == m]
        for _, r in sub.iterrows():
            xj = i + (-0.18 if r["rep"] == "r1" else 0.18)
            ax.scatter(xj, r["promoter_x"], s=64, zorder=3,
                       color=MARK_COLOR[m],
                       marker="o" if r["cond"] == "shNT" else "^",
                       edgecolor="k", linewidth=0.6)
        for rep, xj in (("r1", -0.18), ("r2", 0.18)):
            v = sub[sub["rep"] == rep]["promoter_x"]
            if len(v):
                ax.plot([i + xj - 0.11, i + xj + 0.11], [v.mean()] * 2,
                        color="k", lw=1.4, zorder=4)
    ax.axhline(1, color="k", lw=1.0, ls="--")
    ax.set_xticks(range(len(MARKS)))
    ax.set_xticklabels(MARKS, fontsize=9)
    ax.set_ylabel("promoter enrichment (x)")
    ax.set_title("Per library — left of tick = r1, right = r2\n"
                 "H3K27ac splits 3.3x (r1) vs 6.0x (r2)", fontsize=9.5)
    ax.annotate("", xy=(1.18, 6.0), xytext=(1.18, 3.4),
                arrowprops=dict(arrowstyle="<->", color="#c0392b", lw=1.3))
    ax.text(1.30, 4.7, "batch", fontsize=8, color="#c0392b", rotation=90,
            va="center")
    handles = [Line2D([], [], marker="o", ls="", color="#444", label="shNT"),
               Line2D([], [], marker="^", ls="", color="#444", label="shBCAT1")]
    ax.legend(handles=handles, fontsize=7.6, frameon=False, loc="upper right")

    for a in axes:
        for s in ("top", "right"):
            a.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    fig.savefig(path.replace(".png", ".pdf"))
    plt.close(fig)
    print("[plot]", path)


# ---------------------------------------------------------------- figure 2 ----
def fig_expression(tab, path):
    marks = [m for m in MARKS if m != "IgG"]
    fig, axes = plt.subplots(1, 2, figsize=(12.4, 4.9),
                             gridspec_kw={"width_ratios": [1.0, 1.3]})

    ax = axes[0]
    x = np.arange(len(marks))
    for k, (plat, lab) in enumerate((("wholemRNA", "whole mRNA"),
                                     ("plasmidsaurus", "plasmidsaurus"))):
        vals = [tab[m][f"rho_{plat}"] for m in marks]
        ax.bar(x + (k - 0.5) * 0.38, vals, 0.38,
               color=[MARK_COLOR[m] for m in marks], alpha=1.0 if k == 0 else 0.55,
               edgecolor="k", linewidth=0.5)
    ax.axhline(0, color="k", lw=1.0)
    ax.set_xticks(x)
    ax.set_xticklabels(marks, fontsize=9)
    ax.set_ylabel("Spearman ρ (Δmark vs Δexpression)")
    ax.set_title("Does the promoter change predict expression?\n"
                 "solid = whole mRNA, pale = plasmidsaurus", fontsize=9.5)

    ax = axes[1]
    groups = ["down", "unchanged", "up"]
    w = 0.22
    x = np.arange(len(groups))
    for k, m in enumerate(marks):
        vals = [tab[m][f"mean_delta_{g}"] for g in groups]
        err = [tab[m][f"sem_delta_{g}"] for g in groups]
        ax.bar(x + (k - 1.5) * w, vals, w, yerr=err, capsize=2,
               color=MARK_COLOR[m], label=f"{m} (×{SPIKEIN[m]:.2f}"
                                          f"{'*' if m in ASSUMED else ''})")
    ax.axhline(0, color="k", lw=1.0)
    ax.set_xticks(x)
    ax.set_xticklabels(["down-regulated", "unchanged", "up-regulated"], fontsize=9)
    ax.set_ylabel("mean Δ promoter signal, spike-in scaled\nlog2(shBCAT1/shNT)")
    ax.set_title("Promoter mark change by expression class\n"
                 "(whole-mRNA padj<0.05)", fontsize=9.5)
    ax.legend(fontsize=7.4, frameon=False)
    for a in axes:
        for s in ("top", "right"):
            a.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    fig.savefig(path.replace(".png", ".pdf"))
    plt.close(fig)
    print("[plot]", path)


def main():
    gfrac, gtot = genome_fractions()
    print("[genome] " + "  ".join(f"{c} {gfrac[c]*100:.2f}%" for c in CATS))
    cats = load_categories()

    lines = ["mark,condition,fragments," + ",".join(f"pct_{c}" for c in CATS) + "," +
             ",".join(f"enrich_{c}" for c in CATS)]
    print("\n=== category composition (% of fragments | enrichment vs genome) ===")
    for m in MARKS:
        for cond in ("shNT", "shBCAT1"):
            agg, tot = cats[(m, cond)]
            pct = {c: agg[c] / tot * 100 for c in CATS}
            enr = {c: pct[c] / (gfrac[c] * 100) for c in CATS}
            print(f"  {m:8s} {cond:8s} " + "  ".join(
                f"{c[:4]} {pct[c]:5.2f}% ({enr[c]:4.2f}x)" for c in CATS))
            lines.append(f"{m},{cond},{tot}," +
                         ",".join(f"{pct[c]:.4f}" for c in CATS) + "," +
                         ",".join(f"{enr[c]:.4f}" for c in CATS))
    with open(os.path.join(DATA, "category_composition.csv"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    per_lib = per_library_enrichment(gfrac)
    per_lib.to_csv(os.path.join(DATA, "category_enrichment_per_library.csv"), index=False)
    fig_categories(cats, gfrac, per_lib,
                   os.path.join(PLOTS, "genomic_category_composition.png"))

    # ---------------- promoter deltas ----------------
    core = pd.read_csv(CORE, sep="\t", header=None,
                       names=["chrom", "start", "end", "name", "score", "strand"])
    core["ensembl"] = core.name.str.split("|").str[0].str.split(".").str[0]
    core["symbol"] = core.name.str.split("|").str[1]
    M = load_promoter_matrix(len(core))

    df = core[["chrom", "start", "end", "ensembl", "symbol"]].copy()
    for m in MARKS:
        nt_tot = cats[(m, "shNT")][1]
        bc_tot = cats[(m, "shBCAT1")][1]
        nt = M[(m, "shNT")] / nt_tot * 1e6
        bc = M[(m, "shBCAT1")] / bc_tot * 1e6 * SPIKEIN[m]
        df[f"{m}_shNT_cpm"] = nt
        df[f"{m}_shBCAT1_cpm"] = bc
        eps = 1e-2
        df[f"{m}_delta"] = np.log2((bc + eps) / (nt + eps))
        df[f"{m}_delta_nospike"] = np.log2((bc / SPIKEIN[m] + eps) / (nt + eps))
    df.to_csv(os.path.join(DATA, "promoter_delta_matrix.csv"), index=False)
    print(f"\n[ok] promoter delta matrix: {len(df):,} promoters x {len(MARKS)} marks")

    # ---------------- join expression ----------------
    ex = pd.read_csv(EXPR)
    ex["ensembl"] = ex.ensembl.str.split(".").str[0]
    j = df.merge(ex, on="ensembl", how="inner", suffixes=("", "_ex"))
    # one promoter per gene: keep the window with most total signal
    j["_sig"] = sum(j[f"{m}_shNT_cpm"] + j[f"{m}_shBCAT1_cpm"] for m in MARKS)
    j = j.sort_values("_sig", ascending=False).drop_duplicates("ensembl")
    print(f"[join] {len(j):,} genes with both promoter signal and expression")

    tab, rows = {}, ["mark,spikein,n_genes,rho_wholemRNA,rho_plasmidsaurus,"
                     "mean_delta_down,mean_delta_unchanged,mean_delta_up,"
                     "n_down,n_up,rho_wholemRNA_nospike"]
    print("\n=== promoter Δ vs expression ===")
    for m in [x for x in MARKS if x != "IgG"]:
        d = j[f"{m}_delta"].to_numpy()
        dns = j[f"{m}_delta_nospike"].to_numpy()
        e = {}
        for plat in ("wholemRNA", "plasmidsaurus"):
            ok = j[f"lfc_{plat}"].notna().to_numpy() & np.isfinite(d)
            e[plat] = sp(d[ok], j[f"lfc_{plat}"].to_numpy()[ok])
        okw = j["lfc_wholemRNA"].notna().to_numpy() & np.isfinite(dns)
        rho_ns = sp(dns[okw], j["lfc_wholemRNA"].to_numpy()[okw])

        sig = j["padj_wholemRNA"].to_numpy() < 0.05
        lfc = j["lfc_wholemRNA"].to_numpy()
        grp = {"down": sig & (lfc < 0), "up": sig & (lfc > 0),
               "unchanged": ~sig & np.isfinite(lfc)}
        st = {}
        for g, mask in grp.items():
            v = d[mask & np.isfinite(d)]
            st[f"mean_delta_{g}"] = float(np.mean(v))
            st[f"sem_delta_{g}"] = float(np.std(v) / max(np.sqrt(len(v)), 1))
        tab[m] = {"rho_wholemRNA": e["wholemRNA"],
                  "rho_plasmidsaurus": e["plasmidsaurus"], **st}
        print(f"  {m:8s} ×{SPIKEIN[m]:.2f}  ρ(whole)={e['wholemRNA']:+.3f}  "
              f"ρ(plasmid)={e['plasmidsaurus']:+.3f}  |  mean Δ  "
              f"down {st['mean_delta_down']:+.3f}  unch {st['mean_delta_unchanged']:+.3f}  "
              f"up {st['mean_delta_up']:+.3f}")
        rows.append(f"{m},{SPIKEIN[m]},{len(j)},{e['wholemRNA']:.4f},"
                    f"{e['plasmidsaurus']:.4f},{st['mean_delta_down']:.4f},"
                    f"{st['mean_delta_unchanged']:.4f},{st['mean_delta_up']:.4f},"
                    f"{int(grp['down'].sum())},{int(grp['up'].sum())},{rho_ns:.4f}")
    with open(os.path.join(DATA, "promoter_delta_vs_expression.csv"), "w") as fh:
        fh.write("\n".join(rows) + "\n")
    j.to_csv(os.path.join(DATA, "promoter_delta_with_expression.csv"), index=False)
    fig_expression(tab, os.path.join(PLOTS, "promoter_delta_vs_expression.png"))
    print("[ok] wrote", DATA)


if __name__ == "__main__":
    main()
