"""Bin-conditional genotype comparison of metabolic pathways on the pooled axis.

Scores: control-bin (depth-robust) on the POOLED cells -> common scale for both
genotypes; then residualised on n_genes. NOT z-scored within genotype.

For each pathway:
  naive_delta   marginal K27M-minus-KO mean (ignores pseudotime) -- the confounded
                comparison, shown for contrast
  binadj_delta  mean over pseudotime bins of the within-bin genotype difference,
                where each genotype's bin value is the mean of its libraries
                (so libraries, not cells, are the replication unit) -- controls
                for compositional bias
  wilcox_p      Wilcoxon signed-rank of the per-bin difference across bins
  lib_consist   min over the 2x2 library pairings of sign-agreement of binadj_delta
  shape_r       correlation of the two genotypes' binned profiles (do they move
                the same way along the shared axis)
"""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import sys
import numpy as np, pandas as pd, scipy.sparse as sp, scanpy as sc
from scipy.stats import wilcoxon, spearmanr

B = REPO_ROOT + ""
W = f"{B}/scXenograft/work"; OUT = f"{B}/scXenograft"
N_CTRL, N_BINS_CTRL, NB = 50, 25, 20
rng = np.random.default_rng(0)

# panel
sets = {}
for ln in open(f"{W}/metabolic_panel.gmt"):
    f = ln.rstrip("\n").split("\t"); sets[f[0]] = [x for x in f[2:] if x]
def add(p, pre):
    for ln in open(p):
        f = ln.rstrip("\n").split("\t"); g = [x for x in f[2:] if x]
        if g: sets[pre+f[0]] = g
add(f"{B}/xHSJ019_K27M_vs_KO_analysis/KEGG_metabolism_expanded.gmt", "KEGG_")
hall = {}
for ln in open(f"{B}/signatures/h.all.v7.5.symbols.gmt"):
    f = ln.rstrip("\n").split("\t"); hall[f[0]] = [x for x in f[2:] if x]
for h in ["HALLMARK_CHOLESTEROL_HOMEOSTASIS","HALLMARK_FATTY_ACID_METABOLISM",
          "HALLMARK_OXIDATIVE_PHOSPHORYLATION","HALLMARK_GLYCOLYSIS","HALLMARK_ADIPOGENESIS",
          "HALLMARK_PEROXISOME","HALLMARK_MTORC1_SIGNALING","HALLMARK_HYPOXIA",
          "HALLMARK_BILE_ACID_METABOLISM","HALLMARK_XENOBIOTIC_METABOLISM","HALLMARK_HEME_METABOLISM"]:
    sets["HM_"+h.replace("HALLMARK_","")] = hall[h]

emb = pd.read_csv(f"{W}/pooled_dpt.csv.gz", index_col=0)
if "bin" not in emb.columns:
    emb["bin"] = pd.qcut(emb.dpt, NB, labels=False, duplicates="drop")
A = sc.read_h5ad(f"{W}/hsj019_norm.h5ad")[emb.index]
X = A.X.tocsc() if sp.issparse(A.X) else sp.csc_matrix(A.X)
genes = pd.Index(A.var_names); NG = X.shape[1]
avg = np.asarray(X.mean(0)).ravel(); order = np.argsort(avg)
binsg = np.zeros(NG, int); binsg[order] = np.minimum((np.arange(NG)/NG*N_BINS_CTRL).astype(int), N_BINS_CTRL-1)
members = {b: np.where(binsg == b)[0] for b in range(N_BINS_CTRL)}

names, rs, cs, vs, rc, cc_, vc = [], [], [], [], [], [], []
for nm, gl in sets.items():
    idx = genes.get_indexer(gl); idx = np.unique(idx[idx >= 0])
    if len(idx) < 8: continue
    i = len(names); names.append(nm); k = len(idx)
    rs += [i]*k; cs += list(idx); vs += [1.0/k]*k
    ctrl = []
    for b, cnt in zip(*np.unique(binsg[idx], return_counts=True)):
        pool = np.setdiff1d(members[b], idx); take = min(len(pool), cnt*N_CTRL)
        if take: ctrl.append(rng.choice(pool, take, replace=False))
    ctrl = np.concatenate(ctrl); kc = len(ctrl)
    rc += [i]*kc; cc_ += list(ctrl); vc += [1.0/kc]*kc
S = sp.csr_matrix((vs, (rs, cs)), shape=(len(names), NG))
C = sp.csr_matrix((vc, (rc, cc_)), shape=(len(names), NG))
sco = np.asarray((X @ S.T).todense()) - np.asarray((X @ C.T).todense())
D = pd.DataFrame(sco, index=emb.index, columns=names)
# depth residualise (common model, both genotypes)
ng = emb.n_genes.values.astype(float); M = np.c_[np.ones_like(ng), ng, ng**2]
D = pd.DataFrame(D.values - M @ np.linalg.lstsq(M, D.values, rcond=None)[0],
                 index=emb.index, columns=names)
D.to_csv(f"{W}/pooled_metab_resid.csv.gz")

gt = emb.genotype.values
binv = emb.bin.values
samp = emb["sample"].values
k_libs = sorted(set(samp[gt == "K27M"])); o_libs = sorted(set(samp[gt == "K27M-KO"]))
bins = sorted(emb.bin.dropna().unique())

rows = []
for c in names:
    v = D[c].values
    naive = v[gt == "K27M"].mean() - v[gt == "K27M-KO"].mean()
    # per (library, bin) means
    lib_bin = {}
    for l in k_libs + o_libs:
        for b in bins:
            m = (samp == l) & (binv == b)
            lib_bin[(l, b)] = v[m].mean() if m.sum() >= 10 else np.nan
    kprof = np.array([np.nanmean([lib_bin[(l, b)] for l in k_libs]) for b in bins])
    oprof = np.array([np.nanmean([lib_bin[(l, b)] for l in o_libs]) for b in bins])
    delta = kprof - oprof
    ok = np.isfinite(delta)
    binadj = np.nanmean(delta)
    try: wp = wilcoxon(delta[ok]).pvalue
    except Exception: wp = np.nan
    shape = np.corrcoef(kprof[ok], oprof[ok])[0, 1] if ok.sum() >= 6 else np.nan
    # library-pairing consistency: sign of mean delta for each K-lib vs each O-lib
    signs = []
    for kl in k_libs:
        for ol in o_libs:
            dl = np.array([lib_bin[(kl, b)] - lib_bin[(ol, b)] for b in bins])
            signs.append(np.sign(np.nanmean(dl)))
    lib_consist = 1.0 if len(set(signs)) == 1 else (sum(s == np.sign(binadj) for s in signs)/len(signs))
    rows.append(dict(pathway=c, naive_delta=naive, binadj_delta=binadj,
                     wilcox_p=wp, frac_bins_Khigher=float(np.nanmean(delta > 0)),
                     lib_consist=lib_consist, shape_r=shape))
R = pd.DataFrame(rows).set_index("pathway")
from statsmodels.stats.multitest import multipletests
R["wilcox_fdr"] = multipletests(R.wilcox_p.fillna(1), method="fdr_bh")[1]
R["shift_naive_to_binadj"] = R.naive_delta - R.binadj_delta
R = R.sort_values("binadj_delta")
R.round(4).to_csv(f"{OUT}/pooled_pathway_genotype_binconditional.csv")

print("How much does conditioning on pseudotime change the genotype comparison?")
print(f"  corr(naive_delta, binadj_delta) = "
      f"{spearmanr(R.naive_delta, R.binadj_delta).statistic:+.3f}")
print(f"  pathways that FLIP sign after conditioning: "
      f"{int((np.sign(R.naive_delta) != np.sign(R.binadj_delta)).sum())} of {len(R)}\n")
sig = R[(R.wilcox_fdr < 0.05) & (R.lib_consist == 1.0)]
print(f"bin-adjusted genotype differences, FDR<0.05 & library-consistent: {len(sig)}\n")
print("-- higher in KO at matched pseudotime --")
print(sig.head(12)[["naive_delta","binadj_delta","wilcox_fdr","frac_bins_Khigher","shape_r"]].round(3).to_string())
print("\n-- higher in K27M at matched pseudotime --")
print(sig.tail(12)[["naive_delta","binadj_delta","wilcox_fdr","frac_bins_Khigher","shape_r"]].round(3).to_string())
