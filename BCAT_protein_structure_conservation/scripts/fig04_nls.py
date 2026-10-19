#!/usr/bin/env python3
"""Figure 4 - is there an NLS in BCAT1 (and where BCAT2 differs)?
Basic-residue / charge landscape, live predictor output, and a direct
side-by-side with textbook NLSs."""
import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
from bcat_common import *

P = json.load(open(f"{D}/predictor_results.json"))
R = json.load(open(f"{D}/nls_report.json"))
S = json.load(open(f"{D}/nls_seqs.json"))
SF = json.load(open(f"{D}/struct_features.json"))

BASIC, ACID = set("KR"), set("DE")
NT = {"BCAT1": 25, "BCAT2": 36}
COL = {"BCAT1": C1, "BCAT2": C2}
BLUE, RED = "#2166ac", "#c0392b"


def win(s, w=20):
    n = len(s); h = w // 2
    f = np.zeros(n); c = np.zeros(n)
    for i in range(n):
        a, b = max(0, i - h), min(n, i + h + 1)
        seg = s[a:b]
        f[i] = sum(x in BASIC for x in seg) / len(seg)
        c[i] = sum(x in BASIC for x in seg) - sum(x in ACID for x in seg)
    return f, c


fig = plt.figure(figsize=(7.2, 7.0))
gs = fig.add_gridspec(4, 2, height_ratios=[1.0, 1.0, 1.05, 1.25],
                      hspace=0.72, wspace=0.28,
                      left=0.085, right=0.975, top=0.945, bottom=0.055)

# ------------------------------------------------ a,b: charge landscape
for r, g in enumerate(["BCAT1", "BCAT2"]):
    s = S[g]; n = len(s); x = np.arange(1, n + 1)
    f20, c20 = win(s, 20)
    ax = fig.add_subplot(gs[r, :])
    ax.set_title("ab"[r], loc="left", fontweight="bold")
    ax.axvspan(1, NT[g], color=COL[g], alpha=0.13, lw=0, zorder=0)
    ax.axhline(0, color="#999", lw=0.5)
    ax.fill_between(x, 0, c20, where=c20 >= 0, color=BLUE, alpha=0.55, lw=0, step="mid")
    ax.fill_between(x, 0, c20, where=c20 <= 0, color=RED, alpha=0.55, lw=0, step="mid")
    ax.plot(x, c20, color="#2b2b33", lw=0.6)
    # individual basic / acidic residues as rug marks
    for i, ch in enumerate(s):
        if ch in BASIC:
            ax.plot([i + 1, i + 1], [-9.6, -8.6], color=BLUE, lw=0.5, solid_capstyle="butt")
        elif ch in ACID:
            ax.plot([i + 1, i + 1], [-11.0, -10.0], color=RED, lw=0.5, solid_capstyle="butt")
    ax.text(-0.012, -9.1, "K/R", transform=ax.get_yaxis_transform(), ha="right",
            va="center", fontsize=5.8, color=BLUE)
    ax.text(-0.012, -10.5, "D/E", transform=ax.get_yaxis_transform(), ha="right",
            va="center", fontsize=5.8, color=RED)
    # cNLS Mapper bipartite hits (best score per start position)
    raw = sorted(((pos, sc) for pos, _, sc in P["cNLS_Mapper"][g]["bipartite"]))
    best = {}
    for pos, sc in raw:
        key = next((k for k in best if abs(k - pos) <= 12), pos)
        best[key] = max(best.get(key, 0), sc)
    for pos, sc in best.items():
        ax.plot([pos], [8.6], marker="v", ms=3.0 + sc * 0.45,
                color="#7f7f88" if sc < 3 else ACCENT, clip_on=False, zorder=6)
        if sc >= 3:
            ax.text(pos, 10.2, f"{sc}", ha="center", va="bottom", fontsize=6,
                    color=ACCENT, fontweight="bold", clip_on=False)
    ax.set_ylim(-11.4, 8.4); ax.set_xlim(0, 400)
    ax.set_ylabel("net charge\n(20-aa window)")
    ax.set_xlabel(f"{g} residue")
    ax.spines["bottom"].set_bounds(0, n)
    ax.spines["left"].set_bounds(-8, 8)
    ax.set_yticks([-8, -4, 0, 4, 8])
    ax.text(0.998, 1.05, g, transform=ax.transAxes, ha="right", va="bottom",
            fontsize=8, fontweight="bold", color=COL[g])
    ax.text(NT[g] + 4, -7.2, f"N-term 1–{NT[g]}\nnet charge "
            f"{sum(c in BASIC for c in s[:NT[g]]) - sum(c in ACID for c in s[:NT[g]]):+d}",
            fontsize=6, color=COL[g], va="center", ha="left")
    if r == 0:
        ax.text(0.30, 1.05, "▼ cNLS Mapper bipartite hit (size ∝ score; orange = score ≥ 3)",
                transform=ax.transAxes, ha="left", va="bottom", fontsize=6, color=GREY)
    ax.text(0.0, -0.44, "PSORT II: pat4 none · pat7 none · bipartite none · "
            f"NLS score {P['PSORT_II'][g]['NLS_score']} · predicted cytoplasmic"
            f"   |   NucPred {P['NucPred'][g]}",
            transform=ax.transAxes, ha="left", va="top", fontsize=5.8, color=GREY)

# ------------------------------------------------ c: predictor benchmark
ORDER = ["SV40_largeT", "c-Myc", "Nucleoplasmin", "NPM1_human",
         "BCAT2", "BCAT1", "ALDOA_human", "CA2_human"]
NICE = {"SV40_largeT": "SV40 large T", "c-Myc": "c-Myc", "Nucleoplasmin": "nucleoplasmin",
        "NPM1_human": "NPM1", "BCAT2": "BCAT2", "BCAT1": "BCAT1",
        "ALDOA_human": "aldolase A", "CA2_human": "CA2"}
NPKEY = {"Nucleoplasmin": "Nucleoplasmin_P05221"}
ROLE = {k: R[k]["role"] for k in ORDER}
BARC = {"positive": "#4d4d4d", "negative": "#c9ccd1", "target": None}

ax = fig.add_subplot(gs[2, 0]); ax.set_title("c", loc="left", fontweight="bold")
vals = [P["NucPred"][NPKEY.get(k, k)] for k in ORDER]
cols = [COL[k] if k in COL else BARC[ROLE[k]] for k in ORDER]
ax.barh(range(len(ORDER))[::-1], vals, color=cols, height=0.68, lw=0)
ax.axvline(0.8, color=ACCENT, lw=0.8, ls=(0, (3, 2)))
ax.text(0.8, len(ORDER) - 0.2, "0.8 cut-off", color=ACCENT, fontsize=5.8,
        ha="center", va="bottom")
ax.set_yticks(range(len(ORDER))[::-1]); ax.set_yticklabels([NICE[k] for k in ORDER], fontsize=6)
ax.set_xlim(0, 1.0); ax.set_xlabel("NucPred score")
for i, v in zip(range(len(ORDER))[::-1], vals):
    ax.text(v + 0.02, i, f"{v:.2f}", va="center", fontsize=5.8)

# max %K+R in any 20-aa window
ax = fig.add_subplot(gs[2, 1]); ax.set_title("d", loc="left", fontweight="bold")
vals2 = [100 * R[k]["max_basic_20"] for k in ORDER]
ax.barh(range(len(ORDER))[::-1], vals2, color=cols, height=0.68, lw=0)
ax.set_yticks(range(len(ORDER))[::-1]); ax.set_yticklabels([])
ax.set_xlabel("max % K+R in any 20-aa window")
ax.set_xlim(0, 55)
for i, k in zip(range(len(ORDER))[::-1], ORDER):
    nhit = len(R[k]["pat4"]) + len(R[k]["pat7"]) + len(R[k]["bipartite"])
    ax.text(vals2[ORDER.index(k)] + 1, i, f"{vals2[ORDER.index(k)]:.0f}%   "
            f"pattern hits: {nhit}", va="center", fontsize=5.8)

# ------------------------------------------------ e: sequence comparison
ax = fig.add_subplot(gs[3, :]); ax.set_title("e", loc="left", fontweight="bold")
ROWS = [
    ("SV40 large T 122–134", "STPPKKKRKVED", "classical monopartite (score 16)", "#4d4d4d"),
    ("c-Myc 320–328", "PAAKRVKLD", "classical monopartite", "#4d4d4d"),
    ("nucleoplasmin 155–170", "KRPAATKKAGQAKKKK", "classical bipartite", "#4d4d4d"),
    ("", "", "", None),
    ("BCAT1 1–30", S["BCAT1"][:30], "no NLS predicted by any method", C1),
    ("BCAT2 1–30", S["BCAT2"][:30], "cNLS Mapper bipartite, score 5.8 (25–51)", C2),
]
maxlen = max(len(r[1]) for r in ROWS)
ax.set_xlim(-0.5, maxlen + 12); ax.set_ylim(-0.6, len(ROWS) - 0.4)
ax.invert_yaxis(); ax.axis("off")
for y, (lab, seq, note, col) in enumerate(ROWS):
    if not seq:
        continue
    ax.text(-1.0, y, lab, ha="right", va="center", fontsize=6.2,
            color=col, fontweight="bold" if col in (C1, C2) else "normal")
    for i, ch in enumerate(seq):
        fc = BLUE if ch in BASIC else (RED if ch in ACID else None)
        if fc:
            ax.add_patch(Rectangle((i - 0.45, y - 0.36), 0.9, 0.72, fc=fc,
                                   ec="none", alpha=0.85))
        ax.text(i, y, ch, ha="center", va="center",
                fontdict={"family": "DejaVu Sans Mono", "size": 6.8},
                color="w" if fc else "#16161a")
    net = sum(c in BASIC for c in seq) - sum(c in ACID for c in seq)
    ax.text(maxlen + 0.8, y, f"net {net:+d}   {note}", ha="left", va="center",
            fontsize=6, color=col if col in (C1, C2) else GREY)
# mark the BCAT1 cysteines
b1 = S["BCAT1"][:30]
for i, ch in enumerate(b1):
    if ch == "C":
        ax.plot([i], [3.62], marker="v", ms=2.8, color="#d95f02", clip_on=False)
ax.text(14, 3.50, "C4 / C8 / C12 — cysteine-rich and acidic, not basic", fontsize=5.8,
        color="#d95f02", ha="left", va="bottom")

fig.suptitle("Is the BCAT1 N-terminal extension a nuclear localization signal? "
             "(live output from cNLS Mapper, PSORT II and NucPred, 2026-09-29)",
             y=0.995, fontsize=8.3)
save(fig, "Fig4_BCAT_NLS_prediction")
