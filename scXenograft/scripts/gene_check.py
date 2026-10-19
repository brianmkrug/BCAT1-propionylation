import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import gzip, os
import numpy as np, pandas as pd, scipy.io as sio
D = REPO_ROOT + "/scXenograft/datasets"
OUT = REPO_ROOT + "/scXenograft"
samples = sorted(set(f.split(".")[0] for f in os.listdir(D) if f.endswith(".mtx.gz")))
genes = ["BCAT1","H3F3A","SOX2","VIM","HES1","PTPRZ1","OLIG1","OLIG2","PDGFRA","CSPG4",
         "SOX10","MBP","PLP1","CNP","GFAP","AQP4","S100B","APOE","MKI67","EGFR","DLL3","HOPX","FABP7"]
pc = pd.read_csv(f"{OUT}/scx_percell_qc.csv.gz")
bcsets = {}
rows = []
for s in samples:
    feat = pd.read_csv(f"{D}/{s}.features.tsv.gz", sep="\t", header=None, names=["id","symbol","type"])
    sym = feat.symbol.values
    idx = {g: np.where(sym == "hg19_"+g)[0] for g in genes}
    with gzip.open(f"{D}/{s}.matrix.mtx.gz") as fh:
        M = sio.mmread(fh).tocsr()
    d = pc[pc["sample"] == s]
    bcsets[s] = set(d.barcode)
    hum = d.is_human.values
    hg = feat.id.str.startswith("hg19_").values
    tot = np.asarray(M[hg].sum(0)).ravel()
    r = {"sample": s}
    for g, ii in idx.items():
        if len(ii) == 0:
            r[g] = np.nan; continue
        v = np.asarray(M[ii].sum(0)).ravel()[hum]
        cpm = 1e4*v/np.maximum(tot[hum],1)
        r[g+"_pctdet"] = round(100*float((v>0).mean()),1)
        r[g+"_meanCP10K"] = round(float(cpm.mean()),2)
    rows.append(r)
    print(s, "done", flush=True)
df = pd.DataFrame(rows)
df.to_csv(f"{OUT}/scx_marker_check.csv", index=False)
pd.set_option("display.width",250)
print(df[[c for c in df.columns if c=="sample" or c.endswith("_pctdet")]].to_string(index=False))
print()
print(df[[c for c in df.columns if c=="sample" or c.endswith("_meanCP10K")]].to_string(index=False))
print("\nbarcode overlap hsj019 k27m_r1 vs k27m_r2:",
      len(bcsets['GSM6209151_hsj019_k27m_r1'] & bcsets['GSM6209152_hsj019_k27m_r2']))
