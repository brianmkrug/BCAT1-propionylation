"""Kaplan-Meier survival, HSJ019 K27M vs KO, with 95% CI shading + log-rank p."""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

D = REPO_ROOT + "/scXenograft/datasets"
OUT = REPO_ROOT + "/scXenograft"; FIG = f"{OUT}/figures"
df = pd.read_excel(f"{D}/xHSJ019 survival data.xlsx").iloc[:, :4]
df.columns = ["genotype", "sample", "time", "status"]
df = df.dropna(subset=["time", "status"])
df["genotype"] = df["genotype"].str.upper().str.strip()
K = df[df.genotype == "K27M"]; O = df[df.genotype == "KO"]
print("n K27M =", len(K), "events", int(K.status.sum()),
      "| n KO =", len(O), "events", int(O.status.sum()))

lr = logrank_test(K.time, O.time, K.status, O.status)
p = lr.p_value
med = {}
GT = {"K27M": ("#C1272D", K), "KO": ("#2B6CB0", O)}
fig, ax = plt.subplots(figsize=(6.2, 5))
for name, (col, sub) in GT.items():
    kmf = KaplanMeierFitter()
    kmf.fit(sub.time, sub.status, label=f"{name} (n={len(sub)})")
    kmf.plot_survival_function(ax=ax, ci_show=True, color=col, linewidth=2.0, ci_alpha=0.40)
    med[name] = kmf.median_survival_time_
ax.set_ylim(0, 1.02); ax.set_xlim(0, df.time.max()*1.02)
ax.set_xlabel("time (days)", fontsize=12); ax.set_ylabel("survival probability", fontsize=12)
ax.set_title("HSJ019 xenograft survival — K27M vs K27M-KO", fontsize=13)
ptxt = f"log-rank p = {p:.4f}" if p >= 1e-4 else f"log-rank p = {p:.2e}"
ax.text(0.03, 0.06, ptxt + f"\nmedian: K27M {med['K27M']:.0f} d, KO "
        + ("n.r." if not np.isfinite(med['KO']) else f"{med['KO']:.0f} d"),
        transform=ax.transAxes, fontsize=11, va="bottom")
ax.legend(frameon=False, fontsize=11, loc="upper right")
for sp in ["top", "right"]: ax.spines[sp].set_visible(False)
fig.tight_layout(); fig.savefig(f"{FIG}/hsj019_km_survival.png", dpi=170, bbox_inches="tight")
print(f"log-rank p = {p:.5g}; median K27M {med['K27M']}, KO {med['KO']}")
