"""Stage 5: label Leiden clusters from z-scored Jessa atlas signatures;
cross-check against canonical markers.  Usage: python3 s5_label.py bt245|hsj019

Labelling uses only the neural/glial classes (RGC, GLIP, OPC, OL, ASTRO, EPEN,
NEURON).  NON_NEURAL is retained as a contamination/QC score, not as a candidate
label -- these xenograft cells are human tumour, so a pericyte/microglia call
would indicate a problem, not a cell state.
"""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import sys
import numpy as np, pandas as pd, scanpy as sc

W = REPO_ROOT + "/scXenograft/work"
tag = sys.argv[1]
emb = pd.read_csv(f"{W}/{tag}_embed.csv.gz", index_col=0)
meta = pd.read_csv(f"{W}/{tag}_metascores.csv.gz", index_col=0).loc[emb.index]
sigs = pd.read_csv(f"{W}/{tag}_sigscores.csv.gz", index_col=0).loc[emb.index]
zs = (sigs - sigs.mean()) / sigs.std().replace(0, np.nan)

CAND = ["RGC", "GLIP", "OPC", "OL", "ASTRO", "EPEN", "NEURON"]
cm = meta.groupby(emb.leiden)[CAND + ["NON_NEURAL", "PROLIF_sig"]].mean()
lab = cm[CAND].idxmax(1)
zc = zs.groupby(emb.leiden).mean()
top3 = zc.apply(lambda r: ", ".join(r.nlargest(3).index), axis=1)

MK = ["SOX2","HES1","VIM","FABP7","PTPRZ1","OLIG1","OLIG2","PDGFRA","CSPG4","DLL3",
      "SOX10","BCAS1","MBP","PLP1","MOG","CNP","GFAP","AQP4","S100B","APOE","CLU",
      "MKI67","TOP2A","BCAT1","EGFR"]
A = sc.read_h5ad(f"{W}/{tag}_norm.h5ad")[emb.index]
mk = [g for g in MK if g in A.var_names]
E = pd.DataFrame(np.asarray(A[:, mk].X.todense()), index=emb.index, columns=mk)
me = E.groupby(emb.leiden).mean()

out = pd.DataFrame({"n": emb.leiden.value_counts().sort_index(), "label": lab,
                    "top3_atlas_sig": top3})
out = out.join(cm.round(2)).join(me.round(2))
out = out.join(pd.crosstab(emb.leiden, emb.genotype, normalize="index").mul(100).round(1).add_prefix("pct_"))
out.to_csv(f"{W}/{tag}_cluster_annotation.csv")

pd.set_option("display.width", 320)
print(f"=== {tag} ===")
print(out[["n","label","top3_atlas_sig"] + CAND + ["NON_NEURAL","PROLIF_sig"]
          + [c for c in out.columns if c.startswith("pct_")]].to_string())
print()
show = [g for g in ["SOX2","FABP7","PTPRZ1","OLIG1","OLIG2","PDGFRA","DLL3","SOX10",
                    "BCAS1","MBP","PLP1","MOG","GFAP","AQP4","APOE","MKI67","BCAT1"] if g in me.columns]
print(out[["label"] + show].to_string())

emb["label"] = emb.leiden.map(lab)
emb = emb.join(meta[CAND + ["NON_NEURAL", "PROLIF_sig"]])
emb.to_csv(f"{W}/{tag}_annotated.csv.gz")
print("\nlabel x genotype (% within genotype):")
print((pd.crosstab(emb.label, emb.genotype, normalize="columns")*100).round(1).to_string())
print("\nlabel x sample (% within sample):")
print((pd.crosstab(emb.label, emb["sample"], normalize="columns")*100).round(1).to_string())
