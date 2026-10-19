"""Depth-robust alternative to ssGSEA for the same metabolic panel.

Control-bin scoring (Seurat AddModuleScore / sc.tl.score_genes): mean expression
of the set minus mean of expression-bin-matched control genes. The control term
absorbs the per-cell detection-depth offset that ssGSEA is confounded by.
Same vectorised sparse-product implementation used for the Jessa signatures.
"""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import numpy as np, pandas as pd, scipy.sparse as sp, scanpy as sc
from scipy.stats import spearmanr

W = REPO_ROOT + "/scXenograft/work"
N_CTRL, N_BINS = 50, 25
rng = np.random.default_rng(0)

sets = {}
for ln in open(f"{W}/metabolic_panel.gmt"):
    f = ln.rstrip("\n").split("\t")
    sets[f[0]] = [x for x in f[2:] if x]

for which in ["k27m", "ko13"]:
    e = pd.read_csv(f"{W}/hsj019_{which}_dpt.csv.gz", index_col=0)
    A = sc.read_h5ad(f"{W}/hsj019_norm.h5ad")[e.index]
    X = A.X.tocsc() if sp.issparse(A.X) else sp.csc_matrix(A.X)
    genes = pd.Index(A.var_names); n_genes = X.shape[1]

    avg = np.asarray(X.mean(0)).ravel()
    order = np.argsort(avg)
    bins = np.zeros(n_genes, dtype=int)
    bins[order] = np.minimum((np.arange(n_genes)/n_genes*N_BINS).astype(int), N_BINS-1)
    members = {b: np.where(bins == b)[0] for b in range(N_BINS)}

    names, rs, cs, vs, rc, cc_, vc = [], [], [], [], [], [], []
    for i, (nm, gl) in enumerate(sets.items()):
        idx = genes.get_indexer(gl); idx = np.unique(idx[idx >= 0])
        if len(idx) < 5: continue
        names.append(nm); k = len(idx)
        rs += [len(names)-1]*k; cs += list(idx); vs += [1.0/k]*k
        ctrl = []
        for b, cnt in zip(*np.unique(bins[idx], return_counts=True)):
            pool = np.setdiff1d(members[b], idx)
            take = min(len(pool), cnt*N_CTRL)
            if take: ctrl.append(rng.choice(pool, take, replace=False))
        ctrl = np.concatenate(ctrl); kc = len(ctrl)
        rc += [len(names)-1]*kc; cc_ += list(ctrl); vc += [1.0/kc]*kc

    S = sp.csr_matrix((vs, (rs, cs)), shape=(len(names), n_genes))
    C = sp.csr_matrix((vc, (rc, cc_)), shape=(len(names), n_genes))
    sco = np.asarray((X @ S.T).todense()) - np.asarray((X @ C.T).todense())
    df = pd.DataFrame(sco, index=e.index, columns=names)
    df["CHOL_composite"] = df[[c for c in names if c.startswith("CHOL_")]].mean(1)
    df["FA_composite"] = df[[c for c in names if c.startswith("FA_")]].mean(1)
    df.to_csv(f"{W}/hsj019_{which}_metab_ctrlscore.csv.gz")

    print(f"\n=== {which}: control-bin scores ({len(names)} sets) ===")
    print(f"{'signature':<42} {'rho n_genes':>11} {'rho DPT':>9} {'partial DPT':>12}")
    ng = e.n_genes.values.astype(float); pt = e.dpt_pseudotime.values
    from scipy.stats import rankdata
    def partial(x, y, z):
        ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
        rx, ry, rz = (rankdata(v[ok]) for v in (x, y, z))
        def res(a, b):
            B = np.c_[np.ones_like(b), b]
            return a - B @ np.linalg.lstsq(B, a, rcond=None)[0]
        return np.corrcoef(res(rx, rz), res(ry, rz))[0, 1]
    for c in ["CHOL_composite", "FA_composite",
              "CHOL_Hallmark_cholesterol_homeostasis",
              "CHOL_Reactome_cholesterol_biosynthesis", "CHOL_WP_mevalonate_pathway",
              "FA_Hallmark_fatty_acid_metabolism", "FA_Reactome_mito_beta_oxidation",
              "FA_KEGGc_unsaturated_FA_biosynthesis",
              "CTX_KEGGc_BCAA_degradation", "CTX_KEGGc_propanoate_metabolism"]:
        if c not in df.columns: continue
        v = df[c].values
        print(f"{c[:42]:<42} {spearmanr(ng, v).statistic:>+11.3f} "
              f"{spearmanr(pt, v).statistic:>+9.3f} {partial(pt, v, ng):>+12.3f}")
    cc = df[names].corr(method="spearman")
    iu = np.triu_indices_from(cc, k=1)
    print(f"median pairwise rho between signatures: {np.nanmedian(cc.values[iu]):+.3f} "
          f"(ssGSEA was +0.65 / +0.82)")
