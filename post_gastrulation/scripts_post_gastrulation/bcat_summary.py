import pandas as pd,numpy as np,scipy.sparse as sp
pd.set_option('display.width',300)
c=pd.read_parquet('work/zeng_cells_qc.parquet'); m=pd.read_parquet('work/zeng_marker_labels.parquet'); t=pd.read_parquet('work/zeng_transfer.parquet')
d=pd.concat([c,m,t],axis=1)
tissue={'W3-1':'whole body','W4-1':'whole body','W5-1':'whole body','W4-2':'head','W4-3':'head','W5-2':'head','W5-3':'head','W6-1':'head','W7-1':'head','W8-1':'head','W9-1':'brain','W9-2':'brain','W12-1':'brain'}
d['tissue']=d['sample'].map(tissue); d['week']=d['sample'].str.extract(r'W(\d+)')[0].astype(int)
# OPC gate
Z=sp.load_npz('work/zeng_sub.npz'); G=pd.read_csv('work/geneset_ensembl.csv').iloc[:,0].values
zg=pd.read_parquet('work/zeng_genes.parquet').set_index('ensembl'); sym=pd.Series(np.arange(len(G)),index=zg.loc[G,'symbol'].values)
def pos(g): return (Z[:,sym[g]].toarray().ravel()>0)
d['OPC_like']=pos('PDGFRA')&pos('OLIG2')&(d.marker_coarse=='Neural progenitor')
d.drop(columns=['barcode']).to_parquet('work/zeng_annotated_cells.parquet')
def summ(keys,mn=50):
    g=d.groupby(keys,observed=True); o=pd.DataFrame({'n_cells':g.size(),'BCAT1_pbCPM':g.BCAT1.sum()/g.total.sum()*1e6,'BCAT1_pct_expr':g.BCAT1.apply(lambda x:(x>0).mean()*100),'BCAT2_pbCPM':g.BCAT2.sum()/g.total.sum()*1e6,'MYC_pbCPM':g.MYC.sum()/g.total.sum()*1e6}).round(1)
    return o[o.n_cells>=mn]
A=summ('marker_coarse').sort_values('BCAT1_pbCPM',ascending=False); print(A.to_string())
B=summ(['marker_coarse','week']); print(B['BCAT1_pbCPM'].unstack().round(0).to_string())
C=summ('marker_fine').sort_values('BCAT1_pbCPM',ascending=False); print(C.to_string())
D=summ(['xu_developmental_system','tissue']); 
E=summ('OPC_like',1); print(E.to_string()); print(d[d.OPC_like].groupby('sample').size())
F=summ(['OPC_like','tissue'],1); print(F.to_string())
import os; os.makedirs('out',exist_ok=True)
A.to_csv('out/zeng2023_BCAT1_by_marker_coarse.csv'); B.to_csv('out/zeng2023_BCAT1_by_marker_coarse_week.csv'); C.to_csv('out/zeng2023_BCAT1_by_marker_fine.csv')
D.to_csv('out/zeng2023_BCAT1_by_xu_system_tissue.csv'); summ(['marker_fine','week']).to_csv('out/zeng2023_BCAT1_by_marker_fine_week.csv'); F.to_csv('out/zeng2023_BCAT1_OPC_like.csv')
