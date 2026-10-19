import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import gzip, os, sys
import numpy as np, pandas as pd
import scipy.io as sio
import scipy.sparse as sp

D = REPO_ROOT + "/scXenograft/datasets"
OUT = REPO_ROOT + "/scXenograft"
samples = sorted(set(f.split(".")[0] for f in os.listdir(D) if f.endswith(".mtx.gz")))

rows = []
percell = []
for s in samples:
    feat = pd.read_csv(f"{D}/{s}.features.tsv.gz", sep="\t", header=None,
                       names=["id", "symbol", "type"])
    bc = pd.read_csv(f"{D}/{s}.barcodes.tsv.gz", sep="\t", header=None)[0].values
    with gzip.open(f"{D}/{s}.matrix.mtx.gz") as fh:
        M = sio.mmread(fh).tocsr()   # genes x cells
    hg = feat.id.str.startswith("hg19_").values
    mm = feat.id.str.startswith("mm10_").values
    hgM = M[hg]; mmM = M[mm]
    hg_umi = np.asarray(hgM.sum(0)).ravel()
    mm_umi = np.asarray(mmM.sum(0)).ravel()
    tot = hg_umi + mm_umi
    fhum = np.where(tot > 0, hg_umi / np.maximum(tot, 1), 0)
    hg_genes = np.asarray((hgM > 0).sum(0)).ravel()
    # mito on human side
    sym = feat.symbol[hg].str.replace("hg19_", "", regex=False).values
    mito = np.char.startswith(sym.astype(str), "MT-")
    mt_umi = np.asarray(hgM[mito].sum(0)).ravel()
    pmt = np.where(hg_umi > 0, 100 * mt_umi / np.maximum(hg_umi, 1), 0)

    human = fhum > 0.75
    df = pd.DataFrame(dict(sample=s, barcode=bc, hg19_umi=hg_umi, mm10_umi=mm_umi,
                           frac_hg19=fhum, hg19_genes=hg_genes, pct_mt=pmt,
                           is_human=human))
    percell.append(df)
    h = df[human]
    rows.append(dict(
        sample=s, barcodes=len(bc), n_features=len(feat),
        n_hg19_feat=int(hg.sum()), n_mm10_feat=int(mm.sum()),
        human_cells=int(human.sum()), pct_human=round(100*human.mean(), 1),
        med_hg19_UMI=int(np.median(h.hg19_umi)), med_hg19_genes=int(np.median(h.hg19_genes)),
        med_pct_mt=round(float(np.median(h.pct_mt)), 2),
        genes_detected_human=int(((hgM[:, human] > 0).sum(1) > 0).sum()),
    ))
    print(rows[-1], flush=True)

pd.DataFrame(rows).to_csv(f"{OUT}/scx_sample_summary.csv", index=False)
pc = pd.concat(percell)
pc.to_csv(f"{OUT}/scx_percell_qc.csv.gz", index=False, compression="gzip")
print("done")
