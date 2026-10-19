#!/usr/bin/env python3
"""
make_SUDIPGXIII_motif_logo_table.py

Regenerate the SU-DIPG-XIII K27M-vs-KO promoter-motif sequence-logo table
(Fig. "SUDIPGXIII_motif_logo_table") directly from the DNAFORM / RECLU CAGE
motif-search output.

For the S.K27M-S.KO comparison it:
  1. reads the enriched-motif tables (combined up- and down-regulated),
  2. re-computes Benjamini-Hochberg FDR (q) across the full set of discovered
     motifs (both directions) using the pipeline's binomial enrichment p-values,
  3. selects the top 5 K27M-enriched motifs that carry a known TF (TomTom) and
     the top 5 KO-enriched motifs (ranked by -log10 p),
  4. draws information-content (bits) sequence logos from the de-novo PWMs, and
  5. renders a table (Dir | Motif | logo | Known TF | -log10 p | q(FDR)) with
     per-column auto-fit text, saved as PNG + editable PDF/SVG.

Usage:
    python make_SUDIPGXIII_motif_logo_table.py [CAGE_DATA_DIR] [OUT_DIR]

    CAGE_DATA_DIR  path to the DNAFORM "CAGE data" folder
                   (default: ./CAGE data)
    OUT_DIR        where to write the figure (default: current directory)

Requires: numpy, pandas, matplotlib, logomaker  (pip install logomaker)
"""
import os, sys, csv, math, re
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import logomaker

# text stays editable in vector output (Illustrator-friendly)
matplotlib.rcParams["svg.fonttype"] = "none"   # SVG -> <text> elements
matplotlib.rcParams["pdf.fonttype"] = 42        # PDF -> embedded CID TrueType
matplotlib.rcParams["ps.fonttype"]  = 42

# ---------------------------------------------------------------- paths
CAGE_DIR = sys.argv[1] if len(sys.argv) > 1 else "CAGE data"
OUT_DIR  = sys.argv[2] if len(sys.argv) > 2 else "."
RUN = os.path.join(CAGE_DIR, "All_data", "Result", "RECLU", "Motif_Search",
                   "S.K27M-S.KO", "MOTIF_Run_1")
MTEXT = os.path.join(RUN, "MOTIF_TEXT")
MCLUST = os.path.join(RUN, "MOTIF_CLUSTERING")

# ---------------------------------------------------------------- parsers
def parse_results(path, direction):
    """Enriched-motif table -> list of dicts (motif, consensus, fg, bg, p, tf, dir)."""
    out = []
    for ln in open(path):
        ln = ln.rstrip("\n")
        if not ln or ln.startswith("MotifNo") or ln.startswith("\t"):
            continue
        f = ln.split("\t")
        if len(f) < 5:
            continue
        out.append(dict(motif=f[0], consensus=f[1], fg=int(f[2]), bg=int(f[3]),
                        p=float(f[4]), tf=(f[5] if len(f) > 5 else "NA"),
                        direction=direction))
    return out

def parse_pwm(path):
    """MEME letter-probability matrices -> {motif_name: np.ndarray(w,4)}."""
    d, name, rows, grab = {}, None, [], False
    for ln in open(path):
        s = ln.strip()
        if s.startswith("MOTIF"):
            if name and rows: d[name] = np.array(rows)
            name, rows, grab = s.split()[1], [], False
        elif s.startswith("letter-probability"):
            grab = True
        elif grab:
            p = s.split()
            if len(p) == 4:
                try: rows.append([float(x) for x in p])
                except ValueError: grab = False
            else:
                grab = False
    if name and rows: d[name] = np.array(rows)
    return d

def tfnames(s):
    s = (s or "").strip()
    if s.upper() == "NA" or s == "": return "—"
    return ", ".join(re.sub(r"\(.*?\)", "", x).strip() for x in s.split(","))

def n_tested(folder):
    """Number of discovered/tested motifs = tomtom_out_* dirs (BH denominator)."""
    try:
        return sum(1 for n in os.listdir(folder) if n.startswith("tomtom_out_"))
    except FileNotFoundError:
        return 0

# ---------------------------------------------------------------- load + FDR
res  = (parse_results(os.path.join(MTEXT, "combined_downregulated_result.txt"), "K27M")
        + parse_results(os.path.join(MTEXT, "combined_upregulated_result.txt"),  "KO"))
M = (n_tested(os.path.join(MCLUST, "CombinedDownClusters"))
     + n_tested(os.path.join(MCLUST, "CombinedUpClusters"))) or 94   # 69 + 25

res.sort(key=lambda r: r["p"])                         # Benjamini-Hochberg
for i, r in enumerate(res, 1):
    r["q"] = min(1.0, r["p"] * M / i)
qmin = 1.0
for r in reversed(res):                                # enforce monotonicity
    qmin = min(qmin, r["q"]); r["q"] = qmin

nlp = lambda r: -math.log10(r["p"])
K = [r for r in res if r["direction"] == "K27M" and tfnames(r["tf"]) != "—"]
O = [r for r in res if r["direction"] == "KO"]
K = sorted(K, key=nlp, reverse=True)[:5]               # top 5 K27M-up with known TF
O = sorted(O, key=nlp, reverse=True)[:5]               # top 5 KO-up
sel = [("K27M↑", "#f6dada", r) for r in K] + [("KO↑", "#d8e6f2", r) for r in O]

pwm_down = parse_pwm(os.path.join(MCLUST, "CombinedDownClusters", "All_motif_matrix.txt"))
pwm_up   = parse_pwm(os.path.join(MCLUST, "CombinedUpClusters",   "All_motif_matrix.txt"))
def info_logo(direction, motif):
    pwm = (pwm_down if direction == "K27M" else pwm_up)[motif]
    df = pd.DataFrame(pwm, columns=list("ACGT"))
    return logomaker.transform_matrix(df, from_type="probability", to_type="information")

# ---------------------------------------------------------------- render
cols = ["Dir", "Motif", "Consensus (bits)", "Known TF (TomTom)", "−log₁₀p", "q (FDR)"]
wr   = [0.70, 1.30, 2.05, 4.30, 0.95, 0.95]
nrow = len(sel)
fig = plt.figure(figsize=(12.2, 7.6))
gs = GridSpec(nrow + 1, 6, width_ratios=wr, height_ratios=[0.60] + [1] * nrow,
              hspace=0.16, wspace=0.035, left=0.008, right=0.992, top=0.99, bottom=0.02)

reg = {}
def cell(r, c, text, fc, ha="left", weight="normal", color="black", reg_key=None):
    ax = fig.add_subplot(gs[r, c]); ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_visible(False)
    ax.set_facecolor(fc); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    t = ax.text(0.5 if ha == "center" else 0.06, 0.5, text, ha=ha, va="center",
                color=color, fontweight=weight, fontsize=12, transform=ax.transAxes)
    if reg_key is not None: reg.setdefault(reg_key, []).append((ax, t))

for c, name in enumerate(cols):
    cell(0, c, name, "#333", ha="center", weight="bold", color="white", reg_key=("hdr", c))
for i, (dirlab, tint, r) in enumerate(sel, start=1):
    cell(i, 0, dirlab, tint, ha="center", weight="bold", reg_key="dir")
    cell(i, 1, r["motif"], tint, reg_key="motif")
    axl = fig.add_subplot(gs[i, 2]); axl.set_facecolor("white")
    logomaker.Logo(info_logo(r["direction"], r["motif"]), ax=axl,
                   color_scheme="classic", show_spines=False)
    axl.set_xticks([]); axl.set_yticks([]); axl.set_ylim(0, 2)
    cell(i, 3, tfnames(r["tf"]), tint, reg_key="tf")
    cell(i, 4, f"{nlp(r):.2f}", tint, ha="center", reg_key="stat")
    cell(i, 5, f"{r['q']:.1e}", tint, ha="center", reg_key="stat")

# per-column auto-fit: largest font that fits every cell of the column
fig.canvas.draw(); rend = fig.canvas.get_renderer()
def fitsize(ax, t, cap, wfrac=0.90, hfrac=0.72, lo=5):
    fs = cap
    while fs > lo:
        t.set_fontsize(fs)
        e, b = t.get_window_extent(renderer=rend), ax.get_window_extent()
        if e.width <= b.width * wfrac and e.height <= b.height * hfrac: break
        fs -= 0.5
    return fs
CAP = {"dir": 20, "motif": 20, "tf": 26, "stat": 22}
for key, cells in reg.items():
    if isinstance(key, tuple):                         # header cell (fit individually)
        ax, t = cells[0]; t.set_fontsize(fitsize(ax, t, 13))
    else:
        fs = min(fitsize(ax, t, CAP[key]) for ax, t in cells)
        for ax, t in cells: t.set_fontsize(fs)

os.makedirs(OUT_DIR, exist_ok=True)
stem = os.path.join(OUT_DIR, "SUDIPGXIII_motif_logo_table")
for ext in ("png", "pdf", "svg"):
    fig.savefig(f"{stem}.{ext}", dpi=220, bbox_inches="tight")
print(f"wrote {stem}.png/.pdf/.svg  (BH denominator M={M}; "
      f"{len(K)} K27M-up + {len(O)} KO-up motifs)")
