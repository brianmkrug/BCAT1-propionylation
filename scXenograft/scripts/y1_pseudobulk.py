"""Library-level pseudobulk from the HSJ019 single-cell data.

Unit of replication = the 10x library, matching the bulk design. KO uses r1+r3
(r2 excluded earlier). This gives n=2 vs n=2 -- minimal for DESeq2, see memo.

Also emits a per-library x cell-state pseudobulk for descriptive use only
(pseudo-replication; NOT valid as DESeq2 replicates).
"""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import numpy as np, pandas as pd, scanpy as sc, scipy.sparse as sp

W = REPO_ROOT + "/scXenograft/work"
OUT = REPO_ROOT + "/scXenograft"

A = sc.read_h5ad(f"{W}/hsj019_norm.h5ad")
keep = ((A.obs.genotype == "K27M") |
        ((A.obs.genotype == "K27M-KO") & (A.obs.replicate != "r2")))
A = A[keep].copy()
print("cells used:", A.n_obs)
print(A.obs.groupby(["genotype", "sample"], observed=True).size().to_string())

C = A.layers["counts"]
C = C.tocsr() if sp.issparse(C) else sp.csr_matrix(C)

# ---- library-level pseudobulk ----
labs = A.obs["sample"].astype(str).values
libs = sorted(set(labs))
pb = {}
for l in libs:
    m = labs == l
    pb[l] = np.asarray(C[m].sum(0)).ravel()
pb = pd.DataFrame(pb, index=A.var_names).astype(int)
short = {l: ("K27M-" if "_k27m_" in l else "KO-") + l.rsplit("_", 1)[-1] for l in libs}
pb.columns = [short[c] for c in pb.columns]
pb = pb[[c for c in pb.columns if c.startswith("K27M")] +
        [c for c in pb.columns if c.startswith("KO")]]
pb.index.name = "gene"
pb.to_csv(f"{OUT}/sc_pseudobulk_counts_by_library.csv")
print("\nlibrary pseudobulk:", pb.shape)
print(pb.sum(0).to_frame("total_counts").assign(
    n_cells=[int((labs == [k for k, v in short.items() if v == c][0]).sum())
             for c in pb.columns],
    genes_detected=(pb > 0).sum(0)).to_string())

# ---- library x cell-state pseudobulk (descriptive only) ----
ann = {}
for w in ["k27m", "ko13"]:
    a = pd.read_csv(f"{W}/hsj019_{w}_dpt.csv.gz", index_col=0)
    ann[w] = a["label"]
lab = pd.concat(ann.values())
lab = lab.reindex(A.obs_names)
grp = A.obs["sample"].astype(str) + "|" + lab.fillna("NA").astype(str)
rows = {}
for g in sorted(set(grp)):
    m = (grp == g).values
    if m.sum() >= 30:
        rows[g] = np.asarray(C[m].sum(0)).ravel()
pbs = pd.DataFrame(rows, index=A.var_names).astype(int)
pbs.index.name = "gene"
pbs.to_csv(f"{OUT}/sc_pseudobulk_counts_by_library_state.csv.gz")
print(f"\nlibrary x state pseudobulk: {pbs.shape[1]} groups (>=30 cells)")
