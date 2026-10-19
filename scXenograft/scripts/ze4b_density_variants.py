"""Re-emit the two density variants (full-range binned; KDE filled) as vectors."""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import numpy as np, pandas as pd
from scipy.stats import gaussian_kde
import matplotlib.pyplot as plt
W = REPO_ROOT + "/scXenograft/work"
FIG = REPO_ROOT + "/scXenograft/figures"
e = pd.read_csv(f"{W}/pooled5_dpt.csv.gz", index_col=0)
K = e[e.genotype == "K27M"].dpt.values; O = e[e.genotype == "K27M-KO"].dpt.values
RED, BLUE = "#C1272D", "#2B6CB0"

# full-range binned (80 bins)
edges = np.linspace(0, 1, 81)
fig, ax = plt.subplots(figsize=(13, 4.6))
ax.hist(K, bins=edges, density=True, color=RED, alpha=0.40)
ax.hist(O, bins=edges, density=True, color=BLUE, alpha=0.40)
ax.hist(K, bins=edges, density=True, histtype="step", color=RED, lw=2.0, label=f"K27M (n={len(K)})")
ax.hist(O, bins=edges, density=True, histtype="step", color=BLUE, lw=2.0, label=f"KO r1+r2+r3 (n={len(O)})")
for v, c in [(np.median(K), RED), (np.median(O), BLUE)]: ax.axvline(v, color=c, ls=":", lw=1.3)
ax.set_xlim(0, 1); ax.set_ylim(0, None); ax.set_xlabel("pooled pseudotime"); ax.set_ylabel("density")
ax.legend(frameon=False, fontsize=10, loc="upper right")
ax.set_title("Cell density along pooled pseudotime by genotype (80 bins; overlap = purple; dotted = median)", fontsize=12)
for sp in ["top", "right"]: ax.spines[sp].set_visible(False)
fig.tight_layout(); fig.savefig(f"{FIG}/pooled5_density_binned.png", dpi=170, bbox_inches="tight"); plt.close(fig)

# KDE filled
xs = np.linspace(0, 1, 400)
kK = gaussian_kde(K, bw_method=0.15)(xs); kO = gaussian_kde(O, bw_method=0.15)(xs)
fig, ax = plt.subplots(figsize=(13, 4.6))
ax.fill_between(xs, kK, color=RED, alpha=0.42, lw=0); ax.fill_between(xs, kO, color=BLUE, alpha=0.42, lw=0)
ax.plot(xs, kK, color=RED, lw=2.2, label=f"K27M (n={len(K)})")
ax.plot(xs, kO, color=BLUE, lw=2.2, label=f"KO r1+r2+r3 (n={len(O)})")
for v, c in [(np.median(K), RED), (np.median(O), BLUE)]: ax.axvline(v, color=c, ls=":", lw=1.3)
ax.set_xlim(0, 1); ax.set_ylim(0, None); ax.set_xlabel("pooled pseudotime"); ax.set_ylabel("density")
ax.legend(frameon=False, fontsize=10, loc="upper right")
ax.set_title("Cell density along pooled pseudotime by genotype (KDE; overlap = purple; dotted = median)", fontsize=12)
for sp in ["top", "right"]: ax.spines[sp].set_visible(False)
fig.tight_layout(); fig.savefig(f"{FIG}/pooled5_density_filled.png", dpi=170, bbox_inches="tight"); plt.close(fig)
print("done")
