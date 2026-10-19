"""Filled pseudotime-density per genotype (all 5 libraries). K27M light red,
KO light blue, overlap -> purple via alpha. Wide panel."""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import numpy as np, pandas as pd
from scipy.stats import gaussian_kde
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

W = REPO_ROOT + "/scXenograft/work"
FIG = REPO_ROOT + "/scXenograft/figures"
e = pd.read_csv(f"{W}/pooled5_dpt.csv.gz", index_col=0)
K = e[e.genotype == "K27M"].dpt.values
O = e[e.genotype == "K27M-KO"].dpt.values
NBINS = 80
edges = np.linspace(0, 1, NBINS+1)
fig, ax = plt.subplots(figsize=(6.5, 4.6))   # ~50% narrower
RED, BLUE = "#C1272D", "#2B6CB0"
# filled bars (alpha overlap -> purple) + step outline
ax.hist(K, bins=edges, density=True, color=RED, alpha=0.40, zorder=1)
ax.hist(O, bins=edges, density=True, color=BLUE, alpha=0.40, zorder=1)
ax.hist(K, bins=edges, density=True, histtype="step", color=RED, lw=2.0, zorder=3,
        label=f"K27M (n={len(K)})")
ax.hist(O, bins=edges, density=True, histtype="step", color=BLUE, lw=2.0, zorder=3,
        label=f"KO r1+r2+r3 (n={len(O)})")
# median markers
for v, c in [(np.median(K), RED), (np.median(O), BLUE)]:
    ax.axvline(v, color=c, ls=":", lw=1.3, zorder=2)
XCUT = (0.2, 0.8)   # crop view (histogram computed on full data)
ax.set_xlim(*XCUT); ax.set_ylim(0, None)
ax.set_xlabel("pooled pseudotime", fontsize=11); ax.set_ylabel("density", fontsize=11)
ax.legend(frameon=False, fontsize=10, loc="upper right")
ax.set_title(f"Cell density along pooled pseudotime by genotype "
             f"({NBINS} bins; x-axis cut to {XCUT[0]}–{XCUT[1]}; overlap = purple; dotted = median)",
             fontsize=11)
for sp in ["top", "right"]: ax.spines[sp].set_visible(False)
fig.tight_layout()
fig.savefig(f"{FIG}/pooled5_density_binned_xcut.png", dpi=170, bbox_inches="tight")
print("written; medians K27M %.3f  KO %.3f" % (np.median(K), np.median(O)))
