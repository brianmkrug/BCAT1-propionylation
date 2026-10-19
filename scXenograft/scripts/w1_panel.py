"""Assemble the metabolic signature panel: cholesterol/mevalonate + fatty acid,
from the curated human KEGG-78 list, MSigDB Hallmark, and Reactome/WP.
BCAA and propanoate sets are carried as directly relevant context.
"""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import json, os
import pandas as pd

B = REPO_ROOT + ""
S = f"{B}/signatures"
OUT = REPO_ROOT + "/scXenograft/work"

def read_gmt(p, prefix=""):
    d = {}
    for ln in open(p):
        f = ln.rstrip("\n").split("\t")
        g = [x for x in f[2:] if x]
        if g:
            d[prefix + f[0]] = g
    return d

kegg = read_gmt(f"{B}/xHSJ019_K27M_vs_KO_analysis/KEGG_metabolism_expanded.gmt", "KEGGc_")
hall = read_gmt(f"{S}/h.all.v7.5.symbols.gmt")
c2 = read_gmt(f"{S}/c2.all.v7.5.symbols.gmt")

PANEL = {}
# --- cholesterol / mevalonate ---
CHOL = {
    "CHOL_Hallmark_cholesterol_homeostasis": hall.get("HALLMARK_CHOLESTEROL_HOMEOSTASIS"),
    "CHOL_Reactome_cholesterol_biosynthesis": c2.get("REACTOME_CHOLESTEROL_BIOSYNTHESIS"),
    "CHOL_WP_mevalonate_pathway": c2.get("WP_MEVALONATE_PATHWAY"),
    "CHOL_Reactome_SREBP_regulation": c2.get("REACTOME_REGULATION_OF_CHOLESTEROL_BIOSYNTHESIS_BY_SREBP_SREBF"),
    "CHOL_KEGGc_steroid_biosynthesis": kegg.get("KEGGc_Steroid biosynthesis"),
    "CHOL_KEGGc_terpenoid_backbone": kegg.get("KEGGc_Terpenoid backbone biosynthesis"),
    "CHOL_KEGGc_cholesterol_metabolism": kegg.get("KEGGc_Cholesterol metabolism"),
}
# --- fatty acid ---
FA = {
    "FA_Hallmark_fatty_acid_metabolism": hall.get("HALLMARK_FATTY_ACID_METABOLISM"),
    "FA_Reactome_fatty_acid_metabolism": c2.get("REACTOME_FATTY_ACID_METABOLISM"),
    "FA_Reactome_mito_beta_oxidation": c2.get("REACTOME_MITOCHONDRIAL_FATTY_ACID_BETA_OXIDATION"),
    "FA_WP_fatty_acid_biosynthesis": c2.get("WP_FATTY_ACID_BIOSYNTHESIS"),
    "FA_KEGGc_fatty_acid_biosynthesis": kegg.get("KEGGc_Fatty acid biosynthesis"),
    "FA_KEGGc_fatty_acid_degradation": kegg.get("KEGGc_Fatty acid degradation"),
    "FA_KEGGc_fatty_acid_elongation": kegg.get("KEGGc_Fatty acid elongation"),
    "FA_KEGGc_unsaturated_FA_biosynthesis": kegg.get("KEGGc_Biosynthesis of unsaturated fatty acids"),
}
# --- context: BCAT1-adjacent and general metabolic state ---
CTX = {
    "CTX_KEGGc_BCAA_degradation": kegg.get("KEGGc_Valine, leucine and isoleucine degradation"),
    "CTX_KEGGc_propanoate_metabolism": kegg.get("KEGGc_Propanoate metabolism"),
    "CTX_KEGGc_oxidative_phosphorylation": kegg.get("KEGGc_Oxidative phosphorylation"),
    "CTX_KEGGc_glycolysis": kegg.get("KEGGc_Glycolysis / Gluconeogenesis"),
    "CTX_Hallmark_mTORC1_signaling": hall.get("HALLMARK_MTORC1_SIGNALING"),
    "CTX_Hallmark_hypoxia": hall.get("HALLMARK_HYPOXIA"),
    "CTX_KEGGc_peroxisome": kegg.get("KEGGc_Peroxisome"),
    "CTX_KEGGc_sphingolipid_metabolism": kegg.get("KEGGc_Sphingolipid metabolism"),
    "CTX_KEGGc_glycerophospholipid": kegg.get("KEGGc_Glycerophospholipid metabolism"),
}
for d in (CHOL, FA, CTX):
    for k, v in d.items():
        if v is None:
            print("MISSING:", k)
        else:
            PANEL[k] = sorted(set(v))

with open(f"{OUT}/metabolic_panel.gmt", "w") as fh:
    for k, v in PANEL.items():
        fh.write("\t".join([k, "na"] + v) + "\n")
print(f"\npanel: {len(PANEL)} signatures")
print(pd.Series({k: len(v) for k, v in PANEL.items()}).to_string())
