"""Stage 2: concatenate per-line, normalize per methods (CP10K + natural log),
HVG, PCA. Usage: python3 s2_merge.py BT245|HSJ019
"""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import sys, os, glob
import numpy as np, pandas as pd, scanpy as sc

W = REPO_ROOT + "/scXenograft/work"
line = sys.argv[1]
tag = "bt245" if line == "BT245" else "hsj019"
files = sorted(glob.glob(f"{W}/*_{tag}_*.h5ad"))
print("merging:", [os.path.basename(f) for f in files], flush=True)

adatas = [sc.read_h5ad(f) for f in files]
A = sc.concat(adatas, join="outer", index_unique=None)
A.obs_names_make_unique()
A.var["ensembl"] = adatas[0].var["ensembl"].reindex(A.var_names).values
del adatas
print("merged:", A.shape, flush=True)

sc.pp.filter_genes(A, min_cells=3)
A.layers["counts"] = A.X.copy()

# methods: scale to 10,000 UMIs/cell, natural log
sc.pp.normalize_total(A, target_sum=1e4)
sc.pp.log1p(A)                      # natural log
A.raw = A

sc.pp.highly_variable_genes(A, n_top_genes=2000, flavor="seurat", batch_key="sample")
print("HVGs:", int(A.var.highly_variable.sum()), flush=True)

# cell cycle (Tirosh/Regev, human symbols)
S = """MCM5 PCNA TYMS FEN1 MCM2 MCM4 RRM1 UNG GINS2 MCM6 CDCA7 DTL PRIM1 UHRF1 MLF1IP
HELLS RFC2 RPA2 NASP RAD51AP1 GMNN WDR76 SLBP CCNE2 UBR7 POLD3 MSH2 ATAD2 RAD51 RRM2
CDC45 CDC6 EXO1 TIPIN DSCC1 BLM CASP8AP2 USP1 CLSPN POLA1 CHAF1B BRIP1 E2F8""".split()
G2M = """HMGB2 CDK1 NUSAP1 UBE2C BIRC5 TPX2 TOP2A NDC80 CKS2 NUF2 CKS1B MKI67 TMPO
CENPF TACC3 FAM64A SMC4 CCNB2 CKAP2L CKAP2 AURKB BUB1 KIF11 ANP32E TUBB4B GTSE1 KIF20B
HJURP CDCA3 HN1 CDC20 TTK CDC25C KIF2C RANGAP1 NCAPD2 DLGAP5 CDCA2 CDCA8 ECT2 KIF23
HMMR AURKA PSRC1 ANLN LBR CKAP5 CENPE CTCF NEK2 G2E3 GAS2L3 CBX5 CENPA""".split()
S = [g for g in S if g in A.var_names]; G2M = [g for g in G2M if g in A.var_names]
sc.tl.score_genes_cell_cycle(A, s_genes=S, g2m_genes=G2M)

sc.pp.scale(A, max_value=10)
sc.tl.pca(A, n_comps=50, svd_solver="arpack")
A.X = A.raw.X.copy()          # keep lognorm in .X; scaled matrix not retained
A.write_h5ad(f"{W}/{tag}_pca.h5ad")
print("wrote", f"{tag}_pca.h5ad", A.shape, flush=True)
print(A.obs.groupby(["genotype", "replicate"], observed=True).size().to_string())
