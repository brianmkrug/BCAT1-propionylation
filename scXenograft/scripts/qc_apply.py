import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import numpy as np, pandas as pd
OUT = REPO_ROOT + "/scXenograft"
pc = pd.read_csv(f"{OUT}/scx_percell_qc.csv.gz")
res = []
for s, d in pc.groupby("sample"):
    h = d[d.is_human].copy()
    # paper filters
    keep = (h.pct_mt <= 50) & (h.hg19_genes >= 500) & (h.hg19_umi >= 500)
    h = h[keep]
    for col in ["hg19_genes", "hg19_umi"]:
        m, sd = h[col].mean(), h[col].std()
        h = h[(h[col] > m - 2*sd) & (h[col] < m + 2*sd)]
    res.append(dict(sample=s, barcodes=len(d), human=int(d.is_human.sum()),
                    post_QC=len(h),
                    med_genes=int(h.hg19_genes.median()), med_umi=int(h.hg19_umi.median()),
                    med_pctMT=round(h.pct_mt.median(),1),
                    med_mm10_frac=round(float((1-d[d.is_human].frac_hg19).median()),4)))
r = pd.DataFrame(res)
r["line"] = r["sample"].str.extract(r"_(bt245|hsj019)_")
r["genotype"] = np.where(r["sample"].str.contains("_k27m_"), "K27M", "K27M-KO")
print(r.to_string(index=False))
r.to_csv(f"{OUT}/scx_postQC_summary.csv", index=False)
print("\nBy line/genotype:")
print(r.groupby(["line","genotype"]).post_QC.agg(["sum","count"]).to_string())
