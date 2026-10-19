import numpy as np,pandas as pd,scipy.sparse as sp
Z=sp.load_npz('work/zeng_sub.npz'); G=pd.read_csv('work/geneset_ensembl.csv').iloc[:,0].values
zg=pd.read_parquet('work/zeng_genes.parquet').set_index('ensembl'); sym=zg.loc[G,'symbol'].values
c=pd.read_parquet('work/zeng_cells_qc.parquet'); cl=pd.read_parquet('work/zeng_leiden.parquet').leiden.astype(int).values
K=cl.max()+1; I=sp.csr_matrix((np.ones(len(cl)),(cl,np.arange(len(cl)))),shape=(K,len(cl)))
n=np.asarray(I.sum(1)).ravel(); mean=(I@Z).toarray()/n[:,None]; B=Z.copy(); B.data[:]=1; pct=(I@B).toarray()/n[:,None]
M=pd.DataFrame(mean,columns=sym); P=pd.DataFrame(pct,columns=sym)
M.to_csv('work/zeng_cluster_mean.csv'); P.to_csv('work/zeng_cluster_pct.csv')
comp=pd.crosstab(cl,c['sample'],normalize='index').round(2)
rest=(mean.sum(0)-mean)/(K-1)
for k in range(K):
    fc=mean[k]-rest[k]; top=[sym[i] for i in np.argsort(-fc)[:14]]
    s=comp.loc[k]; samp=', '.join(f'{i}:{v:.2f}' for i,v in s[s>0.08].items())
    print(f'C{k} n={n[k]:.0f} | {" ".join(top)} | {samp}')
