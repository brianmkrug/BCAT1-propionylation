#!/usr/bin/env python3
"""Figures 7B and S7A -- Kaplan-Meier survival on control vs BCAA-restricted diet.

Input:  diet data/diet_survival.xlsx
          sheet '2HIT'  Diet, time, status                (H3.3K27M; Trp53-KO allograft)
          sheet '4HIT'  Diet, time, status, sex           (+ shATRX, PDGFRA o/e)
Output: diet_survival_kaplan_meier.{png,svg,pdf}          Figure 7B  (2HIT above 4HIT)
        diet_survival_kaplan_meier_by_sex.{png,svg,pdf}   Figure S7A (4HIT, female/male)
        diet_survival_stats.csv

status: 1 = event (endpoint reached), 0 = censored.
Kaplan-Meier estimator with Greenwood-variance 95% confidence bands and a two-sided
Mantel-Cox (log-rank) test. Implemented directly so the script has no survival-analysis
dependency; verified against the p values printed on the figures.

Run from anywhere:  python "diet data/plot_diet_survival.py"
"""
import os, json
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats


def repo_root(start=None):
    d = os.path.abspath(start or os.path.dirname(__file__))
    while not os.path.exists(os.path.join(d, "BCAT1_figure_palette.json")):
        p = os.path.dirname(d)
        if p == d:
            raise RuntimeError("repository root not found")
        d = p
    return d


ROOT = repo_root(); HERE = os.path.join(ROOT, "diet data")
PAL = json.load(open(os.path.join(ROOT, "BCAT1_figure_palette.json")))
C = {"Control": PAL["conditions"]["control"]["hex"],
     "BCAAr":   PAL["conditions"]["BCAAr"]["hex"]}

plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans"],
                     "font.size": 9, "axes.linewidth": 0.9, "svg.fonttype": "none",
                     "pdf.fonttype": 42, "ps.fonttype": 42})


def km(time, status):
    """Kaplan-Meier survival with Greenwood 95% CI. Returns step arrays."""
    time = np.asarray(time, float); status = np.asarray(status, int)
    ts = np.unique(time[status == 1])
    S, var, out_t, out_s, lo, hi = 1.0, 0.0, [0.0], [1.0], [1.0], [1.0]
    for t in ts:
        n = (time >= t).sum(); d = ((time == t) & (status == 1)).sum()
        S *= (1 - d / n)
        var += d / (n * (n - d)) if n > d else 0.0
        se = S * np.sqrt(var)
        out_t.append(t); out_s.append(S)
        lo.append(max(0.0, S - 1.96 * se)); hi.append(min(1.0, S + 1.96 * se))
    return (np.array(out_t), np.array(out_s), np.array(lo), np.array(hi))


def logrank(t1, s1, t2, s2):
    """Two-sided Mantel-Cox log-rank test."""
    t = np.r_[t1, t2]; s = np.r_[s1, s2]
    g = np.r_[np.zeros(len(t1), int), np.ones(len(t2), int)]
    O1 = E1 = V = 0.0
    for tt in np.unique(t[s == 1]):
        n, n1 = (t >= tt).sum(), ((t >= tt) & (g == 0)).sum()
        d, d1 = ((t == tt) & (s == 1)).sum(), ((t == tt) & (s == 1) & (g == 0)).sum()
        if n < 2:
            continue
        O1 += d1; E1 += d * n1 / n
        V += d * (n1 / n) * (1 - n1 / n) * (n - d) / (n - 1)
    chi2 = (O1 - E1) ** 2 / V if V > 0 else 0.0
    return chi2, stats.chi2.sf(chi2, 1)


def panel(ax, df, title):
    for diet in ("Control", "BCAAr"):
        d = df[df.Diet == diet]
        t, S, lo, hi = km(d.time, d.status)
        ax.step(t, S, where="post", color=C[diet], lw=1.6,
                label=f"{diet} (n={len(d)})", zorder=3)
        ax.fill_between(t, lo, hi, step="post", color=C[diet], alpha=0.18,
                        linewidth=0, zorder=2)
    a = df[df.Diet == "Control"]; b = df[df.Diet == "BCAAr"]
    chi2, p = logrank(a.time, a.status, b.time, b.status)
    ax.text(0.03, 0.10, "Log-rank\n" + ("p < 0.0001" if p < 1e-4 else f"p = {p:.4f}"),
            transform=ax.transAxes, fontsize=8, va="bottom", linespacing=1.35)
    ax.set_title(title, fontsize=9)
    ax.set_ylim(-0.03, 1.05); ax.set_xlim(0, max(df.time) * 1.06)
    ax.set_ylabel("Survival probability")
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(width=0.9, length=3)
    ax.legend(frameon=False, fontsize=8, loc="lower left",
              bbox_to_anchor=(0.02, 0.22), handlelength=1.4)
    return chi2, p, len(a), len(b), a.time.median(), b.time.median()


rows = []
h2 = pd.read_excel(os.path.join(HERE, "diet_survival.xlsx"), sheet_name="2HIT").dropna(how="all")
h4 = pd.read_excel(os.path.join(HERE, "diet_survival.xlsx"), sheet_name="4HIT").dropna(how="all")
for d in (h2, h4):
    d["time"] = pd.to_numeric(d["time"]); d["status"] = pd.to_numeric(d["status"])

# ---- Figure 7B
fig, axes = plt.subplots(2, 1, figsize=(3.5, 5.0), sharex=True,
                         gridspec_kw=dict(hspace=0.30))
for ax, df, ttl, tag in [(axes[0], h2, "H3.3K27M, Trp53 KO allograft", "2HIT"),
                         (axes[1], h4, "H3.3K27M, Trp53 KO, shATRX,\nPDGFRA o/e allograft", "4HIT")]:
    chi2, p, n1, n2, m1, m2 = panel(ax, df, ttl)
    rows.append(dict(figure="Figure 7B", cohort=tag, stratum="all", n_control=n1,
                     n_BCAAr=n2, median_control=m1, median_BCAAr=m2, chi2=chi2, p=p))
axes[1].set_xlabel("Days after tumor implantation")
for e, kw in (("png", dict(dpi=400)), ("svg", {}), ("pdf", {})):
    fig.savefig(os.path.join(HERE, f"diet_survival_kaplan_meier.{e}"),
                bbox_inches="tight", **kw)
plt.close(fig); print("wrote diet_survival_kaplan_meier.{png,svg,pdf}")

# ---- Figure S7A
fig, axes = plt.subplots(2, 1, figsize=(3.5, 5.0), sharex=True,
                         gridspec_kw=dict(hspace=0.30))
for ax, sex in zip(axes, ("Female", "Male")):
    sub = h4[h4.sex == sex]
    chi2, p, n1, n2, m1, m2 = panel(ax, sub, f"4HIT: {sex}")
    rows.append(dict(figure="Figure S7A", cohort="4HIT", stratum=sex, n_control=n1,
                     n_BCAAr=n2, median_control=m1, median_BCAAr=m2, chi2=chi2, p=p))
axes[1].set_xlabel("Days after tumor implantation")
for e, kw in (("png", dict(dpi=400)), ("svg", {}), ("pdf", {})):
    fig.savefig(os.path.join(HERE, f"diet_survival_kaplan_meier_by_sex.{e}"),
                bbox_inches="tight", **kw)
plt.close(fig); print("wrote diet_survival_kaplan_meier_by_sex.{png,svg,pdf}")

st = pd.DataFrame(rows)
st["test"] = "two-sided Mantel-Cox (log-rank) test, control vs BCAAr"
st.to_csv(os.path.join(HERE, "diet_survival_stats.csv"), index=False)
print(); print(st.to_string(index=False))
