"""Stage 4: score Jessa 2022 (braindex) atlas signatures on merged objects.

Scoring = Seurat AddModuleScore / sc.tl.score_genes logic (mean signature
expression minus mean of expression-bin-matched control genes), vectorised so
all 290 signatures are computed in one sparse matrix product.

Usage: python3 s4_score.py bt245|hsj019
"""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import sys, re, json
import numpy as np, pandas as pd, scipy.sparse as sp, scanpy as sc

W = REPO_ROOT + "/scXenograft/work"
GMT = REPO_ROOT + "/xHSJ019_K27M_vs_KO_analysis/gsea_braindex/gene_sets.gmt"
tag = sys.argv[1]
rng = np.random.default_rng(0)
N_CTRL, N_BINS = 50, 25

A = sc.read_h5ad(f"{W}/{tag}_norm.h5ad")
X = A.X.tocsc() if sp.issparse(A.X) else sp.csc_matrix(A.X)   # cells x genes, lognorm
genes = pd.Index(A.var_names)
n_cells, n_genes = X.shape

sigs = {}
for ln in open(GMT):
    p = ln.rstrip("\n").split("\t")
    g = [x for x in p[2:] if x]
    sigs[p[0]] = g

# expression bins for control matching
avg = np.asarray(X.mean(0)).ravel()
order = np.argsort(avg)
bins = np.zeros(n_genes, dtype=int)
bins[order] = np.minimum((np.arange(n_genes) / n_genes * N_BINS).astype(int), N_BINS - 1)
bin_members = {b: np.where(bins == b)[0] for b in range(N_BINS)}

names, rows_s, cols_s, vals_s, rows_c, cols_c, vals_c, cover = [], [], [], [], [], [], [], []
for i, (nm, gl) in enumerate(sigs.items()):
    idx = genes.get_indexer(gl)
    idx = idx[idx >= 0]
    if len(idx) < 5:
        continue
    names.append(nm); cover.append((nm, len(gl), len(idx)))
    k = len(idx)
    rows_s += [i]*k; cols_s += list(idx); vals_s += [1.0/k]*k
    ctrl = []
    for b, cnt in zip(*np.unique(bins[idx], return_counts=True)):
        pool = np.setdiff1d(bin_members[b], idx, assume_unique=False)
        take = min(len(pool), cnt * N_CTRL)
        if take:
            ctrl.append(rng.choice(pool, take, replace=False))
    ctrl = np.concatenate(ctrl)
    kc = len(ctrl)
    rows_c += [i]*kc; cols_c += list(ctrl); vals_c += [1.0/kc]*kc

n_sig = len(sigs)
S = sp.csr_matrix((vals_s, (rows_s, cols_s)), shape=(n_sig, n_genes))
C = sp.csr_matrix((vals_c, (rows_c, cols_c)), shape=(n_sig, n_genes))
sc_mat = np.asarray((X @ S.T).todense()) - np.asarray((X @ C.T).todense())
keep = np.array([i for i, nm in enumerate(sigs) if nm in set(names)])
scores = pd.DataFrame(sc_mat[:, keep], index=A.obs_names, columns=names)
scores.to_csv(f"{W}/{tag}_sigscores.csv.gz")
pd.DataFrame(cover, columns=["signature", "n_genes", "n_found"]).to_csv(
    f"{W}/{tag}_sig_coverage.csv", index=False)
print(f"{tag}: scored {len(names)}/{n_sig} signatures over {n_cells} cells")
print("median gene recovery:", int(np.median([c[2] for c in cover])),
      "of", int(np.median([c[1] for c in cover])))

# ---- lineage meta-scores: average z-scored signature scores within class ----
def cls(nm):
    s = re.sub(r"^[FP]-[a-z0-9]+_", "", nm)
    if re.match(r"^(V|D|G)?RGC", s): return "RGC"
    if s.startswith("OPC"): return "OPC"
    if s.startswith(("NFOL", "MOL")): return "OL"
    if s.startswith(("ASTR", "ASTP", "OPAS")): return "ASTRO"
    if s.startswith("GLIP"): return "GLIP"
    if s.startswith("EPEN"): return "EPEN"
    if s.startswith(("CEXN", "PEXN", "NEUR", "CINHN", "PINH", "MGINH", "GABAN",
                     "EXIP", "INIP", "UBCN", "SERN", "ACHN", "THLN", "MFN",
                     "CJRZN", "SMSN")): return "NEURON"
    if s.startswith(("PERI", "ENDO", "VSMC", "VLM", "MNG", "MGL", "MAC",
                     "CPLX", "SCHW", "OTIC", "CHEM", "ASEP", "FLRPT")): return "NON_NEURAL"
    return "OTHER"

z = (scores - scores.mean()) / scores.std().replace(0, np.nan)
meta = pd.DataFrame({c: z[[n for n in names if cls(n) == c]].mean(1)
                     for c in ["RGC", "GLIP", "OPC", "OL", "ASTRO", "EPEN",
                               "NEURON", "NON_NEURAL"]})
meta["top_class"] = meta.idxmax(1)
# proliferating-signature score (any '-P' suffixed cluster)
meta["PROLIF_sig"] = z[[n for n in names if n.endswith("-P")]].mean(1)
meta.to_csv(f"{W}/{tag}_metascores.csv.gz")
print("\nmeta-class sizes (n signatures):",
      {c: len([n for n in names if cls(n) == c]) for c in meta.columns[:8]})
