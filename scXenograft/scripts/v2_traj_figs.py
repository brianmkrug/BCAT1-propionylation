"""Trajectory figures: PAGA graph, pseudotime UMAP, meta-score and gene trends
along pseudotime, with per-library curves as the replicate check."""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx

W = REPO_ROOT + "/scXenograft/work"
FIG = REPO_ROOT + "/scXenograft/figures"
CLS = ["RGC", "GLIP", "OPC", "OL", "ASTRO", "EPEN", "NEURON"]
LCOL = dict(zip(sorted(CLS), plt.cm.Set2(np.linspace(0, 1, len(CLS)))))
SETS = [("k27m", "HSJ019 K27M (r1+r2)"), ("ko13", "HSJ019 K27M-KO (r1+r3)")]

D = {w: pd.read_csv(f"{W}/hsj019_{w}_dpt.csv.gz", index_col=0) for w, _ in SETS}
P = {w: np.load(f"{W}/hsj019_{w}_paga_pos.npz", allow_pickle=True) for w, _ in SETS}

# ---------- Fig 1: PAGA graphs ----------
fig, axes = plt.subplots(1, 2, figsize=(14, 6.5))
for ax, (w, ttl) in zip(axes, SETS):
    conn, names = P[w]["conn"], list(P[w]["names"])
    d = D[w]
    sizes = d.leiden_harmony.astype(str).value_counts()
    G = nx.Graph()
    for i, n in enumerate(names):
        G.add_node(n)
    for i in range(len(names)):
        for j in range(i+1, len(names)):
            if conn[i, j] > 0.05:
                G.add_edge(names[i], names[j], weight=conn[i, j])
    pos = nx.spring_layout(G, weight="weight", seed=1, k=1.1)
    ew = [G[u][v]["weight"] for u, v in G.edges()]
    nx.draw_networkx_edges(G, pos, ax=ax, width=[3*w_/max(ew) for w_ in ew],
                           edge_color="#999999", alpha=.7)
    ns = [sizes.get(n.split(":")[0], 10)/6 for n in G.nodes()]
    nc = [LCOL[n.split(":")[1]] for n in G.nodes()]
    nx.draw_networkx_nodes(G, pos, ax=ax, node_size=ns, node_color=nc,
                           edgecolors="k", linewidths=.5)
    nx.draw_networkx_labels(G, pos, ax=ax, font_size=6.5)
    ax.set_title(f"{ttl}\nPAGA (edges > 0.05)", fontsize=10); ax.axis("off")
h = [plt.Line2D([], [], marker="o", ls="", color=LCOL[c], label=c) for c in sorted(CLS)]
axes[1].legend(handles=h, frameon=False, fontsize=8, loc="upper right")
fig.tight_layout(); fig.savefig(f"{FIG}/hsj019_paga_graphs.png", dpi=170); plt.close(fig)

# ---------- Fig 2: pseudotime on UMAP + trends ----------
fig, axes = plt.subplots(2, 3, figsize=(16, 9.5))
for row, (w, ttl) in enumerate(SETS):
    d = D[w]
    ax = axes[row, 0]
    o = np.argsort(d.dpt_pseudotime.values)
    sm = ax.scatter(d.UMAP1_harmony.values[o], d.UMAP2_harmony.values[o],
                    c=d.dpt_pseudotime.values[o], s=4, cmap="viridis",
                    linewidths=0, rasterized=True)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_visible(False)
    ax.set_title(f"{ttl}\nDPT pseudotime", fontsize=10)
    fig.colorbar(sm, ax=ax, fraction=.045)

    # binned meta-score trends
    ax = axes[row, 1]
    bins = np.linspace(0, 1, 26)
    ctr = (bins[:-1] + bins[1:]) / 2
    idx = np.digitize(d.dpt_pseudotime, bins) - 1
    for c in ["RGC", "OPC", "OL", "ASTRO", "PROLIF_sig"]:
        m = [d[c][idx == i].mean() if (idx == i).sum() > 15 else np.nan
             for i in range(25)]
        ax.plot(ctr, m, lw=2, label=c,
                color=LCOL.get(c, "#444444"), ls="--" if c == "PROLIF_sig" else "-")
    ax.axhline(0, lw=.6, color="k"); ax.set_xlabel("DPT pseudotime")
    ax.set_ylabel("mean meta-score"); ax.legend(fontsize=7, frameon=False)
    ax.set_title("lineage scores along pseudotime", fontsize=10)
    for sp in ["top", "right"]: ax.spines[sp].set_visible(False)

    # BCAT1 + markers, per library
    ax = axes[row, 2]
    libs = sorted(d["sample"].unique())
    lc = dict(zip(libs, [plt.cm.Dark2(0.0), plt.cm.Dark2(0.6)]))
    for g, ls in [("BCAT1", "-"), ("PDGFRA", ":"), ("MBP", "--")]:
        for l in libs:
            s = d[d["sample"] == l]
            ii = np.digitize(s.dpt_pseudotime, bins) - 1
            m = [s[g][ii == i].mean() if (ii == i).sum() > 15 else np.nan
                 for i in range(25)]
            ax.plot(ctr, m, ls=ls, lw=1.6, color=lc[l],
                    label=f"{g} — {l.split('_',1)[1]}")
    ax.set_xlabel("DPT pseudotime"); ax.set_ylabel("mean lognorm expression")
    ax.legend(fontsize=6, frameon=False, ncol=2)
    ax.set_title("BCAT1 / PDGFRA / MBP per library", fontsize=10)
    for sp in ["top", "right"]: ax.spines[sp].set_visible(False)
fig.tight_layout(); fig.savefig(f"{FIG}/hsj019_pseudotime_trends.png", dpi=170)
plt.close(fig)

# ---------- Fig 3: pseudotime density per library ----------
fig, axes = plt.subplots(1, 2, figsize=(12, 4.2))
for ax, (w, ttl) in zip(axes, SETS):
    d = D[w]
    for l in sorted(d["sample"].unique()):
        s = d[d["sample"] == l]
        h_, edges = np.histogram(s.dpt_pseudotime, bins=30, range=(0, 1), density=True)
        ax.plot((edges[:-1]+edges[1:])/2, h_, lw=2, label=f"{l.split('_',1)[1]} (n={len(s)})")
    ax.set_xlabel("DPT pseudotime"); ax.set_ylabel("density")
    ax.set_title(ttl, fontsize=10); ax.legend(fontsize=8, frameon=False)
    for sp in ["top", "right"]: ax.spines[sp].set_visible(False)
fig.suptitle("Pseudotime distribution per library (replicate consistency check)", fontsize=12)
fig.tight_layout(); fig.savefig(f"{FIG}/hsj019_pseudotime_density.png", dpi=170)
plt.close(fig)
print("figures written")

for w, ttl in SETS:
    d = D[w]
    print(f"\n{ttl}: median pseudotime per library")
    print(d.groupby("sample").dpt_pseudotime.agg(["median", "mean", "count"]).round(3).to_string())
