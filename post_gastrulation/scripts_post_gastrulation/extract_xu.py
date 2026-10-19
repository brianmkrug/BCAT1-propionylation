import numpy as np,pandas as pd,scipy.sparse as sp
G=pd.read_csv('work/geneset_ensembl.csv').iloc[:,0].values
g=pd.read_csv('xu/GSE157329_gene_annotate.txt.gz',sep='\t'); pos=pd.Series(np.arange(len(g)),index=g.gene_id)
remap=-np.ones(len(g),int); remap[pos.reindex(G).values]=np.arange(len(G))
c=pd.read_csv('xu_cells_BCAT.csv.gz'); tot=c.tot.values
R=[];C=[];V=[]
for ch in pd.read_csv('xu/GSE157329_raw_counts.mtx.gz',sep=' ',comment='%',skiprows=2,header=None,names=['g','c','v'],dtype={'g':np.int32,'c':np.int32,'v':np.float32},chunksize=20_000_000):
    gg=remap[ch.g.values-1]; k=gg>=0; ci=ch.c.values[k]-1
    R.append(ci); C.append(gg[k]); V.append(np.log1p(ch.v.values[k]*1e4/tot[ci]).astype(np.float32))
M=sp.csr_matrix((np.concatenate(V),(np.concatenate(R),np.concatenate(C))),shape=(len(c),len(G)))
sp.save_npz('work/xu_sub.npz',M); print('Xu',M.shape,M.nnz)
