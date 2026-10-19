"""Score the metabolic panel (same 108-set list as the 4-lib pooled analysis)
on the 5-library pooled cells; control-bin + depth-residualise. -> pooled5_metab_resid"""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import numpy as np, pandas as pd, scipy.sparse as sp, scanpy as sc
B = REPO_ROOT + ""; W = f"{B}/scXenograft/work"
N_CTRL, N_BINS = 50, 25; rng = np.random.default_rng(0)

sets = {}
for ln in open(f"{W}/metabolic_panel.gmt"):
    f = ln.rstrip("\n").split("\t"); sets[f[0]] = [x for x in f[2:] if x]
for ln in open(f"{B}/xHSJ019_K27M_vs_KO_analysis/KEGG_metabolism_expanded.gmt"):
    f = ln.rstrip("\n").split("\t"); g = [x for x in f[2:] if x]
    if g: sets["KEGG_"+f[0]] = g
hall = {}
for ln in open(f"{B}/signatures/h.all.v7.5.symbols.gmt"):
    f = ln.rstrip("\n").split("\t"); hall[f[0]] = [x for x in f[2:] if x]
for h in ["HALLMARK_CHOLESTEROL_HOMEOSTASIS","HALLMARK_FATTY_ACID_METABOLISM",
          "HALLMARK_OXIDATIVE_PHOSPHORYLATION","HALLMARK_GLYCOLYSIS","HALLMARK_ADIPOGENESIS",
          "HALLMARK_PEROXISOME","HALLMARK_MTORC1_SIGNALING","HALLMARK_HYPOXIA",
          "HALLMARK_BILE_ACID_METABOLISM","HALLMARK_XENOBIOTIC_METABOLISM","HALLMARK_HEME_METABOLISM"]:
    sets["HM_"+h.replace("HALLMARK_","")] = hall[h]

emb = pd.read_csv(f"{W}/pooled5_dpt.csv.gz", index_col=0)
A = sc.read_h5ad(f"{W}/hsj019_norm.h5ad")[emb.index]
X = A.X.tocsc() if sp.issparse(A.X) else sp.csc_matrix(A.X)
genes = pd.Index(A.var_names); NG = X.shape[1]
avg = np.asarray(X.mean(0)).ravel(); order = np.argsort(avg)
bg = np.zeros(NG, int); bg[order] = np.minimum((np.arange(NG)/NG*N_BINS).astype(int), N_BINS-1)
members = {b: np.where(bg == b)[0] for b in range(N_BINS)}

names, rs, cs, vs, rc, cc_, vc = [], [], [], [], [], [], []
for nm, gl in sets.items():
    idx = genes.get_indexer(gl); idx = np.unique(idx[idx >= 0])
    if len(idx) < 5: continue
    i = len(names); names.append(nm); k = len(idx)
    rs += [i]*k; cs += list(idx); vs += [1.0/k]*k
    ctrl = []
    for b, cnt in zip(*np.unique(bg[idx], return_counts=True)):
        pool = np.setdiff1d(members[b], idx); take = min(len(pool), cnt*N_CTRL)
        if take: ctrl.append(rng.choice(pool, take, replace=False))
    ctrl = np.concatenate(ctrl); kc = len(ctrl)
    rc += [i]*kc; cc_ += list(ctrl); vc += [1.0/kc]*kc
S = sp.csr_matrix((vs,(rs,cs)), shape=(len(names),NG)); C = sp.csr_matrix((vc,(rc,cc_)), shape=(len(names),NG))
sco = np.asarray((X @ S.T).todense()) - np.asarray((X @ C.T).todense())
D = pd.DataFrame(sco, index=emb.index, columns=names)
ng = emb.n_genes.values.astype(float); M = np.c_[np.ones_like(ng), ng, ng**2]
D = pd.DataFrame(D.values - M @ np.linalg.lstsq(M, D.values, rcond=None)[0], index=emb.index, columns=names)
D.to_csv(f"{W}/pooled5_metab_resid.csv.gz")
print(f"pooled5_metab_resid: {D.shape} ({D.shape[1]} pathways over {D.shape[0]} cells)")
