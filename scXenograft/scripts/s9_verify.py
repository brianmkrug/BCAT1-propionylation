"""Stage 9: replicate-level verification. Are any genotype effects consistent
across every HSJ019 replicate? Uses per-library means (pseudobulk of scores),
which is the level at which the design actually has replication."""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import numpy as np, pandas as pd, scanpy as sc
W = REPO_ROOT + "/scXenograft/work"
OUT = REPO_ROOT + "/scXenograft"
CLS = ["RGC", "GLIP", "OPC", "OL", "ASTRO", "EPEN", "NEURON", "PROLIF_sig"]

rows = []
for tag in ["bt245", "hsj019"]:
    d = pd.read_csv(f"{W}/{tag}_annotated.csv.gz", index_col=0)
    A = sc.read_h5ad(f"{W}/{tag}_norm.h5ad")[d.index]
    for g in ["BCAT1", "OLIG1", "OLIG2", "PDGFRA", "MBP", "PLP1", "GFAP", "APOE", "MKI67"]:
        if g in A.var_names:
            d[g] = np.asarray(A[:, g].X.todense()).ravel()
    pb = d.groupby(["genotype", "sample"])[CLS + [g for g in
         ["BCAT1","OLIG1","OLIG2","PDGFRA","MBP","PLP1","GFAP","APOE","MKI67"]
         if g in d.columns]].mean()
    pb["line"] = tag
    rows.append(pb)
    print(f"\n=== {tag}: per-library means ===")
    print(pb.round(2).to_string())

    if tag == "hsj019":
        k = d[d.genotype == "K27M"].groupby("sample")
        o = d[d.genotype == "K27M-KO"].groupby("sample")
        print("\nconsistency check (does every KO library differ from every K27M "
              "library in the same direction?)")
        for c in CLS + ["BCAT1", "OLIG1", "PDGFRA", "APOE"]:
            if c not in d.columns: continue
            kv, ov = k[c].mean().values, o[c].mean().values
            sep = "K27M > KO" if kv.min() > ov.max() else (
                  "KO > K27M" if ov.min() > kv.max() else "OVERLAPS")
            print(f"  {c:<12} K27M {np.round(kv,2)}  KO {np.round(ov,2)}   -> {sep}")

pd.concat(rows).to_csv(f"{OUT}/scx_pseudobulk_scores_by_library.csv")
print("\nwrote scx_pseudobulk_scores_by_library.csv")
