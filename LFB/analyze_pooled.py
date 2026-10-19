#!/usr/bin/env python3
"""
WT vs Bcat1-KO comparison of ROI-based LFB optical density, by region.

Input : deconvolved/LFB data.xlsx  (genotype, mouse, Region, replicate, mean)

Produces two matched figures:
    lfb_barplot.png / .svg          all 9 mice
    lfb_barplot_noWT1.png / .svg    WT-1 excluded (the one statistically flagged
                                    outlier; it is HIGH, so removing it is the
                                    conservative direction)
plus lfb_pooled_stats.csv / lfb_pooled_stats_noWT1.csv

Bars are mean +/- SEM of MOUSE means; dots are individual mice (5 ROIs each).
Displayed p = Welch t-test on mouse means, Holm-corrected across the four regions.
Mann-Whitney, a mixed model on all 180 ROIs, and the naive (pseudoreplicated)
ROI-level t-test are also written to the CSV for reference.
"""
import os
import json
import warnings
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.formula.api as smf
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")
# --- vector text settings -------------------------------------------------
# svg.fonttype "none"  -> SVG keeps real <text> elements referencing the font by
#                         name, so strings stay editable (not converted to paths).
# pdf.fonttype 42      -> embed TrueType (Type 42) rather than Type 3, so the PDF
#                         opens as editable type boxes in Illustrator.
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42

B = os.path.join(os.path.dirname(os.path.abspath(__file__)), "deconvolved")
REGIONS = ["CORTEX", "OB", "CB", "PONS"]
LABELS = {"CORTEX": "Cortex", "OB": "Olfactory bulb",
          "CB": "Cerebellum", "PONS": "Pons"}

pal = json.load(open(os.path.join(os.path.dirname(B), "..",
                                  "BCAT1_figure_palette.json")))["genotypes"]
WT_C, KO_C = pal["WT"]["hex"], pal["KO"]["hex"]

raw = pd.read_excel(os.path.join(B, "LFB data.xlsx")).rename(columns={"mean": "lfb"})
raw["id"] = raw.genotype + "-" + raw.mouse.astype(str)


def holm(p):
    p = np.asarray(p, float); n = len(p); order = np.argsort(p)
    adj = np.empty(n); run = 0.0
    for i, idx in enumerate(order):
        run = max(run, (n - i) * p[idx]); adj[idx] = min(run, 1.0)
    return adj


def fmt_p(p):
    if p < 1e-4:
        m, e = f"{p:.1e}".split("e")
        return f"p = {m}×10$^{{{int(e)}}}$"
    if p < 0.001:
        return f"p = {p:.4f}"
    return f"p = {p:.3f}"


def analyse(d):
    mouse = d.groupby(["genotype", "id", "Region"], as_index=False)["lfb"].mean()
    rows = []
    for r in REGIONS:
        mm = mouse[mouse.Region == r]
        w, k = mm[mm.genotype == "WT"].lfb.values, mm[mm.genotype == "KO"].lfb.values
        tt = stats.ttest_ind(w, k, equal_var=False)
        uu = stats.mannwhitneyu(w, k, alternative="two-sided")
        sp = np.sqrt(((len(w)-1)*w.var(ddof=1) + (len(k)-1)*k.var(ddof=1))
                     / (len(w)+len(k)-2))
        J = 1 - 3/(4*(len(w)+len(k)) - 9)
        g = (w.mean() - k.mean()) / sp * J

        sub = d[d.Region == r].copy()
        md = smf.mixedlm("lfb ~ genotype", sub, groups=sub["id"]).fit(reml=False)
        coef = [c for c in md.params.index if c.startswith("genotype")][0]
        mdr = smf.mixedlm("lfb ~ genotype", sub, groups=sub["id"]).fit(reml=True)
        gv, rv = float(mdr.cov_re.iloc[0, 0]), float(mdr.scale)
        icc = gv / (gv + rv) if (gv + rv) > 0 else 0.0

        roi = d[d.Region == r]
        p_naive = stats.ttest_ind(roi[roi.genotype == "WT"].lfb,
                                  roi[roi.genotype == "KO"].lfb, equal_var=False).pvalue

        rows.append(dict(
            region=r,
            WT_mean=round(w.mean(), 4), WT_sd=round(w.std(ddof=1), 4),
            WT_sem=round(w.std(ddof=1)/np.sqrt(len(w)), 4), n_WT=len(w),
            KO_mean=round(k.mean(), 4), KO_sd=round(k.std(ddof=1), 4),
            KO_sem=round(k.std(ddof=1)/np.sqrt(len(k)), 4), n_KO=len(k),
            KO_over_WT=round(k.mean()/w.mean(), 3), hedges_g=round(g, 2),
            p_mouse_welch=round(tt.pvalue, 5),
            p_mouse_mannwhitney=round(uu.pvalue, 5),
            ICC_between_mouse=round(icc, 3),
            p_mixedmodel=float("%.3g" % md.pvalues[coef]),
            lmm_trustworthy=("no - random effect at boundary" if icc < 0.05 else "yes"),
            p_naive_ROI_PSEUDOREP=float("%.2e" % p_naive)))

    res = pd.DataFrame(rows)
    for col in ["p_mouse_welch", "p_mouse_mannwhitney", "p_mixedmodel"]:
        res[col.replace("p_", "pholm_")] = np.round(holm(res[col].values), 5)
    res["PRIMARY_pholm"] = res["pholm_mouse_welch"]
    return res, mouse


def figure(res, mouse, path_stem, subtitle):
    fig, ax = plt.subplots(figsize=(8.2, 5.2))
    rng = np.random.default_rng(3)
    x = np.arange(len(REGIONS)); wdt = 0.36

    for j, (gt, c) in enumerate((("WT", WT_C), ("KO", KO_C))):
        mu = [res[res.region == r][f"{gt}_mean"].iloc[0] for r in REGIONS]
        se = [res[res.region == r][f"{gt}_sem"].iloc[0] for r in REGIONS]
        pos = x + (j - 0.5) * wdt
        ax.bar(pos, mu, wdt*0.92, color=c, alpha=.85, zorder=2,
               label="WT" if gt == "WT" else "Bcat1-KO")
        ax.errorbar(pos, mu, yerr=se, fmt="none", ecolor="0.25",
                    capsize=4, lw=1.3, zorder=5)
        for i, r in enumerate(REGIONS):
            v = mouse[(mouse.Region == r) & (mouse.genotype == gt)].lfb.values
            ax.scatter(pos[i] + rng.uniform(-.09, .09, len(v)), v, s=27,
                       color="white", edgecolor="0.25", linewidth=.85, zorder=4)

    for i, r in enumerate(REGIONS):
        row = res[res.region == r].iloc[0]
        y = max(mouse[mouse.Region == r].lfb.max(),
                row.WT_mean + row.WT_sem) * 1.06
        ax.plot([x[i]-wdt/2, x[i]+wdt/2], [y, y], color="0.3", lw=1.1)
        ax.text(x[i], y*1.03, fmt_p(row.PRIMARY_pholm), ha="center", fontsize=9.5)

    ax.set_xticks(x)
    ax.set_xticklabels([LABELS[r] for r in REGIONS], fontsize=10.5)
    ax.set_ylabel("Mean LFB optical density", fontsize=10.5)
    ax.set_title(subtitle, fontsize=9.5)
    ax.legend(frameon=False, fontsize=10, loc="upper center",
              bbox_to_anchor=(0.5, -0.08), ncol=2, handlelength=1.4)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_ylim(0, YMAX)
    fig.tight_layout()
    fig.savefig(path_stem + ".png", dpi=220, bbox_inches="tight")
    fig.savefig(path_stem + ".svg", bbox_inches="tight")
    plt.close(fig)



CM = 1/2.54


def figure_compact(res, mouse, path_stem, w_cm=4.0, h_cm=2.0):
    """Publication-panel version: exact physical size, no title, scaled type.

    NOTE: at 4 x 2 cm the type lands at ~3.5-4.5 pt, below the 5-7 pt minimum most
    journals specify. The SVG is vector with editable text, so the intended workflow
    is to place it at this size and bump the font up in Illustrator if needed.
    """
    fig = plt.figure(figsize=(w_cm*CM, h_cm*CM))
    ax = fig.add_axes([0.135, 0.215, 0.855, 0.745])   # exact margins, no bbox_inches
    rng = np.random.default_rng(3)
    x = np.arange(len(REGIONS)); wdt = 0.38
    FS_TICK, FS_LAB, FS_P, FS_LEG = 3.4, 3.6, 3.2, 3.2

    for j, (gt, c) in enumerate((("WT", WT_C), ("KO", KO_C))):
        mu = [res[res.region == r][f"{gt}_mean"].iloc[0] for r in REGIONS]
        se = [res[res.region == r][f"{gt}_sem"].iloc[0] for r in REGIONS]
        pos = x + (j - 0.5) * wdt
        ax.bar(pos, mu, wdt*0.9, color=c, linewidth=0, zorder=2,
               label="WT" if gt == "WT" else "Bcat1-KO")
        ax.errorbar(pos, mu, yerr=se, fmt="none", ecolor="0.2",
                    capsize=0.9, capthick=0.4, elinewidth=0.4, zorder=5)
        for i, r in enumerate(REGIONS):
            v = mouse[(mouse.Region == r) & (mouse.genotype == gt)].lfb.values
            ax.scatter(pos[i] + rng.uniform(-.10, .10, len(v)), v, s=1.6,
                       color="white", edgecolor="0.25", linewidth=0.22, zorder=4)

    for i, r in enumerate(REGIONS):
        row = res[res.region == r].iloc[0]
        y = max(mouse[mouse.Region == r].lfb.max(), row.WT_mean + row.WT_sem) * 1.10
        ax.plot([x[i]-wdt/2, x[i]+wdt/2], [y, y], color="0.3", lw=0.4)
        ax.text(x[i], y*1.045, fmt_p(row.PRIMARY_pholm).replace(" ", ""),
                ha="center", fontsize=FS_P)

    ax.set_xticks(x)
    ax.set_xticklabels([LABELS[r].replace(" ", "\n") for r in REGIONS], fontsize=FS_TICK)
    ax.set_yticks([0, 0.1, 0.2, 0.3, 0.4])
    ax.set_yticklabels(["0", "0.1", "0.2", "0.3", "0.4"], fontsize=FS_TICK)
    ax.set_ylabel("LFB optical density", fontsize=FS_LAB, labelpad=1.5)
    ax.tick_params(axis="both", length=1.3, width=0.4, pad=1.0)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_linewidth(0.4)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_ylim(0, YMAX)
    ax.legend(frameon=False, fontsize=FS_LEG, loc="upper left",
              bbox_to_anchor=(-0.015, 1.05), handlelength=0.7, handleheight=0.5,
              handletextpad=0.3, labelspacing=0.12, borderpad=0.0)

    fig.savefig(path_stem + ".svg")
    fig.savefig(path_stem + ".pdf")
    fig.savefig(path_stem + ".png", dpi=1200)
    plt.close(fig)
    return fig


# shared y-axis limit so the two panels are directly comparable
YMAX = raw.groupby(["genotype", "id", "Region"])["lfb"].mean().max() * 1.22

CONFIGS = [
    (None, "lfb_barplot", "lfb_pooled_stats.csv",
     "LFB myelin signal by region — all mice (WT n=5, KO n=4)"),
    (["WT-1"], "lfb_barplot_noWT1", "lfb_pooled_stats_noWT1.csv",
     "LFB myelin signal by region — WT-1 excluded (WT n=4, KO n=4)"),
]

for excl, stem, csvname, head in CONFIGS:
    d = raw if not excl else raw[~raw.id.isin(excl)].copy()
    res, mouse = analyse(d)
    res.to_csv(os.path.join(B, csvname), index=False)
    sub = (head + "\n"
           "bars = mean $\\pm$ SEM of mouse means; dots = individual mice (5 ROIs each)\n"
           "Welch t-test on mouse means, Holm-corrected across regions")
    figure(res, mouse, os.path.join(B, stem), sub)
    if excl == ["WT-1"]:
        figure_compact(res, mouse, os.path.join(B, stem + "_4x2cm"))
    print("=== " + head + " ===")
    print(res[["region", "WT_mean", "WT_sem", "KO_mean", "KO_sem", "KO_over_WT",
               "hedges_g", "p_mouse_welch", "PRIMARY_pholm"]].to_string(index=False))
    print()
