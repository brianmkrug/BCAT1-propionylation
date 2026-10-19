import numpy as np,pandas as pd,json,scipy.sparse as sp
import statsmodels.api as sm
exec(open('scripts/metab_score.py').read().split("S=pscore(pb)")[0])
gr=json.load(open('work/metab_groups.json')); OPC='OPC (PDGFRA+ OLIG2+)'
M=sp.load_npz('work/brain_counts_all.npz').tocsr(); tot=np.asarray(M.sum(1)).ravel()
L=sp.diags(1e4/tot)@M; L.data=np.log1p(L.data); L=L.tocsc()
cells=ab.grp.isin([OPC]+gr['neural']).values
L=L[cells]; a=ab[cells].reset_index(drop=True); tot=tot[cells]
mean=np.asarray(L.mean(0)).ravel(); expr=np.where(mean>0)[0]
bins=pd.qcut(pd.Series(mean[expr]).rank(method='first'),25,labels=False).values; binof=pd.Series(bins,index=expr)
rng=np.random.default_rng(0)
def score(idx):
    idx=[i for i in idx if i in binof.index]
    ctrl=set()
    for i in idx:
        pool=expr[bins==binof[i]]; ctrl|=set(rng.choice(pool,min(10,len(pool)),replace=False))
    ctrl=list(ctrl-set(idx))
    return np.asarray(L[:,idx].mean(1)).ravel()-np.asarray(L[:,ctrl].mean(1)).ravel()
S_genes='MCM5 PCNA TYMS FEN1 MCM2 MCM4 RRM1 UNG GINS2 MCM6 CDCA7 DTL PRIM1 UHRF1 HELLS RFC2 RPA2 NASP RAD51AP1 GMNN WDR76 SLBP CCNE2 UBR7 POLD3 MSH2 ATAD2 RAD51 RRM2 CDC45 CDC6 EXO1 TIPIN DSCC1 BLM CASP8AP2 USP1 CLSPN POLA1 CHAF1B BRIP1 E2F8'.split()
G2M_genes='HMGB2 CDK1 NUSAP1 UBE2C BIRC5 TPX2 TOP2A NDC80 CKS2 NUF2 CKS1B MKI67 TMPO CENPF TACC3 FAM64A SMC4 CCNB2 CKAP2L CKAP2 AURKB BUB1 KIF11 ANP32E TUBB4B GTSE1 KIF20B HJURP CDCA3 HN1 CDC20 TTK CDC25C KIF2C RANGAP1 NCAPD2 DLGAP5 CDCA2 CDCA8 ECT2 KIF23 HMMR AURKA PSRC1 ANLN LBR CKAP5 CENPE CTCF NEK2 G2E3 GAS2L3 CBX5 CENPA'.split()
a['S']=score([s2i[g] for g in S_genes if g in s2i]); a['G2M']=score([s2i[g] for g in G2M_genes if g in s2i])
a['opc']=(a.grp==OPC).astype(float); a['ldepth']=np.log10(tot)
a['cycling']=(a[['S','G2M']].max(1)>np.quantile(a[['S','G2M']].max(1),0.8))
print('cycling frac OPC',a[a.opc==1].cycling.mean().round(2),'others',a[a.opc==0].cycling.mean().round(2))
sc={k:score(v) for k,v in PI.items() if len(v)>=3}
SC=pd.DataFrame(sc); SC.to_parquet('work/metab_cell_scores.parquet'); a.to_parquet('work/metab_cell_meta.parquet')
X=pd.get_dummies(a[['opc','S','G2M','ldepth','sample']],columns=['sample'],drop_first=True).astype(float); X=sm.add_constant(X)
prog=a.marker_coarse.isin(['Neural progenitor']).values|(a.opc==1).values
res=[]
for k in SC:
    y=(SC[k]-SC[k].mean())/SC[k].std()
    f1=sm.OLS(y,X).fit()
    nc=~a.cycling.values; f2=sm.OLS(y[nc],X[nc]).fit()
    f3=sm.OLS(y[prog],X[prog]).fit()
    res.append((k,f1.params.opc,f1.tvalues.opc,f2.params.opc,f2.tvalues.opc,f3.params.opc,f3.tvalues.opc,f1.params.G2M))
R=pd.DataFrame(res,columns=['pathway','beta_all','t_all','beta_noncycling','t_noncycling','beta_vs_progenitors','t_vs_progenitors','G2M_coef']).set_index('pathway')
R.round(3).sort_values('beta_all',ascending=False).to_csv('work/metab_cell_ols.csv')
pd.set_option('display.width',200); print(R.round(2).sort_values('beta_all',ascending=False).to_string())
