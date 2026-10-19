#!/usr/bin/env python3
"""
Global integrated heatmap with PER-DATASET z-scoring.
Each metabolite is z-scored within cell culture and within xenograft SEPARATELY
(removes the CC-vs-xeno batch offset so genotype/cell-line structure dominates).
The two blocks are drawn side by side with a white gap. Outlier HSJ019_K27M_2 removed.
No sample names. Rows clustered on the combined per-dataset-z matrix; columns clustered
within each dataset block. RdBu_r. Top annotation strips + super-pathway row bar + legends.
"""
import os, re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42

import matplotlib.colors as mcolors
from matplotlib.patches import Patch
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(HERE)
XLSX = os.path.join(BASE, "norm_cell_data copy.xlsx")
FIG = os.path.join(HERE, "figures")
META = ["Cell line", "Genotype", "Source"]
DROP_LABEL = ("HSJ019", "K27M", 2)

DATASET_C = {"Cell culture": "#762a83", "Xenograft": "#e08214"}
GENO_C = {"K27M": "#d1495b", "KO": "#2e86ab", "hfNSC (WT)": "#5f8f4e"}
LINE_C = {"BT245": "#66c2a5", "DIPGXIII": "#fc8d62", "HSJ019": "#8da0cb", "hfNSC": "#e78ac3"}
SP_ORDER = ["Amino Acid", "Peptide", "Carbohydrate", "Energy", "Lipid",
            "Nucleotide", "Cofactors and Vitamins", "Xenobiotics",
            "Partially Characterized Molecules"]
SP_C = dict(zip(SP_ORDER, plt.get_cmap("tab10").colors[:len(SP_ORDER)]))
GAP = 1  # blank columns between blocks


def norm(s): return re.sub(r"\s+", " ", str(s).strip().lower())


def load():
    cc = pd.read_excel(XLSX, sheet_name="cell culture")
    xe = pd.read_excel(XLSX, sheet_name="xenograft")
    # drop xeno outlier
    cnt, keep = {}, []
    for _, r in xe.iterrows():
        k = (r["Cell line"], r["Genotype"]); cnt[k] = cnt.get(k, 0) + 1
        keep.append(not (r["Cell line"] == DROP_LABEL[0] and
                         r["Genotype"] == DROP_LABEL[1] and cnt[k] == DROP_LABEL[2]))
    xe = xe.loc[keep].reset_index(drop=True)
    ca = pd.read_excel(os.path.join(BASE, "cell_annotation.xlsx"))
    xa = pd.read_excel(os.path.join(BASE, "xenograft_annotation.xlsx"))
    sp_map, sub_map = {}, {}
    for a in (ca, xa):
        for _, r in a.iterrows():
            sp_map.setdefault(r["Metabolite"], r["Super pathway"])
            sub_map.setdefault(r["Metabolite"], r["Sub pathway"])
    return cc, xe, sp_map, sub_map


def detected(df):
    cols = [c for c in df.columns if c not in META]
    return set(m for m in cols if (pd.to_numeric(df[m], errors="coerce").fillna(0) > 0).any())


def block(df, union, dsname):
    """Return per-dataset z-scored matrix (metabolites x samples) + sample meta."""
    cols, meta_rows, counts = {}, [], {}
    for i, row in df.iterrows():
        line, geno = row["Cell line"], row["Genotype"]
        counts[(line, geno)] = counts.get((line, geno), 0) + 1
        sid = f"{dsname[:2]}|{line}_{geno}_{counts[(line, geno)]}"
        cols[sid] = pd.to_numeric(df.loc[i, union], errors="coerce").values
        meta_rows.append((sid, dsname, "hfNSC (WT)" if geno == "WT" else geno, line))
    M = pd.DataFrame(cols, index=union)
    Mz = M.replace(0, np.nan)
    half = Mz.min(axis=1, skipna=True) / 2.0
    Ml = np.log10(Mz.T.fillna(half).T)
    Z = Ml.sub(Ml.mean(1), axis=0).div(Ml.std(1, ddof=1).replace(0, np.nan), axis=0)
    meta = pd.DataFrame(meta_rows, columns=["sid", "Dataset", "Genotype", "Cell line"]).set_index("sid")
    return Z, meta


def main():
    cc, xe, sp_map, sub_map = load()
    union = sorted(detected(cc) & detected(xe))
    BCAA_SUB = "Leucine, Isoleucine and Valine Metabolism"
    Zcc, mcc = block(cc, union, "Cell culture")
    Zxe, mxe = block(xe, union, "Xenograft")
    Zcc = Zcc.fillna(0); Zxe = Zxe.fillna(0)

    # cluster columns within each block
    occ = leaves_list(linkage(Zcc.T.values, method="ward"))
    oxe = leaves_list(linkage(Zxe.T.values, method="ward"))
    Zcc, mcc = Zcc.iloc[:, occ], mcc.iloc[occ]
    Zxe, mxe = Zxe.iloc[:, oxe], mxe.iloc[oxe]

    # cluster rows on combined per-dataset-z matrix
    comb = np.hstack([Zcc.values, Zxe.values])
    Lrow = linkage(comb, method="ward")
    orow = leaves_list(Lrow)

    ncc, nxe = Zcc.shape[1], Zxe.shape[1]
    ntot = ncc + GAP + nxe
    nrow = len(union)
    # assemble data matrix with NaN gap, rows in cluster order (origin lower -> reverse)
    ro = orow[::-1]
    data = np.full((nrow, ntot), np.nan)
    data[:, :ncc] = Zcc.values[ro]
    data[:, ncc + GAP:] = Zxe.values[ro]

    sp = pd.Series({m: sp_map.get(m, "Partially Characterized Molecules") for m in union}).reindex(union)
    sp_ord = sp.values[ro]
    bcaa = pd.Series({m: (sub_map.get(m) == BCAA_SUB) for m in union}).reindex(union)
    bcaa_ord = bcaa.values[ro]
    meta_all = pd.concat([mcc, mxe])
    col_sids = list(Zcc.columns) + [None] * GAP + list(Zxe.columns)

    # ---- figure layout: fixed 7 cm x 10 cm canvas, minimum font 6 pt ----
    from matplotlib.colors import ListedColormap
    CM = 1 / 2.54
    FS, TFS = 6, 6.5                                   # min body / title pt
    plt.rcParams.update({"svg.fonttype": "none",
                         "lines.linewidth": 0.4,      # dendrogram links
                         "axes.linewidth": 0.4})      # annotation/heatmap box outlines
    SP_ABBR = {"Amino Acid": "Amino Acid", "Peptide": "Peptide",
               "Carbohydrate": "Carbohydrate", "Energy": "Energy", "Lipid": "Lipid",
               "Nucleotide": "Nucleotide", "Cofactors and Vitamins": "Cofactors/Vit",
               "Xenobiotics": "Xenobiotics",
               "Partially Characterized Molecules": "Partially char."}

    fig = plt.figure(figsize=(8.55 * CM, 10 * CM))     # canvas grows to hold the wider heatmap
    hx0, hx1, hy0, hy1 = 0.151, 0.617, 0.045, 0.78    # heatmap region (a further +20% wide; legends kept at right)
    HW, HH = hx1 - hx0, hy1 - hy0
    ax = fig.add_axes([hx0, hy0, HW, HH])
    cmap = plt.get_cmap("RdBu_r").copy(); cmap.set_bad("white")
    ax.imshow(data, aspect="auto", cmap=cmap, vmin=-2.5, vmax=2.5,
              origin="lower", extent=[0, ntot, 0, nrow], interpolation="nearest")
    ax.set_xticks([]); ax.set_yticks([])

    # ---- row super-pathway strip ----
    axsp = fig.add_axes([hx0 - 0.023, hy0, 0.0147, HH])
    sp_rgb = np.array([[mcolors.to_rgb(SP_C[s]) for s in sp_ord]]).transpose(1, 0, 2)
    axsp.imshow(sp_rgb, aspect="auto", origin="lower", extent=[0, 1, 0, nrow])
    axsp.set_xticks([]); axsp.set_yticks([])

    # ---- BCAA-pathway row marker strip (right of heatmap) ----
    axbc = fig.add_axes([hx1 + 0.005, hy0, 0.0152, HH])   # BCAA marker column 50% wider
    bccol = np.where(bcaa_ord, 0.0, np.nan).reshape(nrow, 1)
    bcmap = ListedColormap(["#111111"]); bcmap.set_bad("white")
    axbc.imshow(bccol, aspect="auto", origin="lower", extent=[0, 1, 0, nrow],
                cmap=bcmap, vmin=0, vmax=1, interpolation="nearest")
    axbc.set_xticks([]); axbc.set_yticks([])
    axbc.text(0.5, 1.004, "BCAA", transform=axbc.transAxes, ha="center",
              va="bottom", fontsize=FS, rotation=90)
    for s in axbc.spines.values():
        s.set_linewidth(0.3)

    # ---- row dendrogram ----
    axd = fig.add_axes([hx0 - 0.098, hy0, 0.0615, HH])
    dendrogram(Lrow, orientation="left", ax=axd, no_labels=True, color_threshold=0,
               link_color_func=lambda k: "#555555")
    for coll in axd.collections:
        coll.set_linewidth(0.35)
    axd.set_xticks([]); axd.set_yticks([])
    for s in axd.spines.values():
        s.set_visible(False)
    axd.invert_yaxis()

    # ---- top annotation strips ----
    def color_row(values, cmap_dict):
        return np.array([[(1, 1, 1) if v is None else mcolors.to_rgb(cmap_dict[v])
                          for v in values]])
    ann = [("Cell line", [None if s is None else meta_all.loc[s, "Cell line"] for s in col_sids], LINE_C),
           ("Genotype", [None if s is None else meta_all.loc[s, "Genotype"] for s in col_sids], GENO_C),
           ("Dataset", [None if s is None else meta_all.loc[s, "Dataset"] for s in col_sids], DATASET_C)]
    strip_h = 0.013
    for i, (name, vals, cmapd) in enumerate(ann):
        ay = hy1 + 0.004 + i * (strip_h + 0.004)
        axa = fig.add_axes([hx0, ay, HW, strip_h])
        axa.imshow(color_row(vals, cmapd), aspect="auto", extent=[0, ntot, 0, 1],
                   interpolation="nearest")
        axa.set_xticks([]); axa.set_yticks([])
        axa.text(-0.02, 0.5, name, transform=axa.transAxes, va="center", ha="right", fontsize=FS)

    # ---- legends (right column, 6 pt, compact) ----
    LX = 0.658
    def leg(handles, title, y):
        lg = fig.legend(handles=handles, title=title, loc="upper left",
                        bbox_to_anchor=(LX, y), frameon=False, fontsize=FS,
                        title_fontsize=TFS, handlelength=0.9, handletextpad=0.4,
                        labelspacing=0.25, alignment="left")
        fig.add_artist(lg)
    leg([Patch(facecolor=c, label=k) for k, c in DATASET_C.items()], "Dataset", 0.90)
    leg([Patch(facecolor=c, label=k) for k, c in GENO_C.items()], "Genotype", 0.79)
    leg([Patch(facecolor=c, label=k) for k, c in LINE_C.items()], "Cell line", 0.645)
    leg([Patch(facecolor=SP_C[k], label=SP_ABBR[k]) for k in SP_ORDER], "Super pathway", 0.475)
    leg([Patch(facecolor="#111111", label="BCAA pathway")], "Marker", 0.13)

    # ---- colorbar: top-left gap under the title, label on top ----
    axcb = fig.add_axes([0.045, 0.905, 0.14, 0.013])
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=mcolors.Normalize(-2.5, 2.5), cmap=cmap),
                      cax=axcb, orientation="horizontal")
    cb.set_label("z-score (log10)", fontsize=FS, labelpad=2)
    cb.ax.xaxis.set_label_position("top")
    cb.ax.tick_params(labelsize=FS, width=0.3, length=1.5)

    for ext in ("png", "svg", "pdf"):
        out = os.path.join(FIG, f"global_heatmap_splitz.{ext}")
        fig.savefig(out, dpi=300)          # no tight bbox -> exact 7.85x10 cm canvas
        print("wrote", out)
    plt.close(fig)
    print("ncc", ncc, "nxe", nxe)


if __name__ == "__main__":
    main()
