import h5py,numpy as np,pandas as pd,scipy.sparse as sp
G=pd.read_csv('work/geneset_ensembl.csv').iloc[:,0].values
# Zeng
f=h5py.File('zeng/human_raw.h5ad','r'); X=f['X']; ip=X['indptr'][:]; n=len(ip)-1
ens=np.array([s.decode() for s in f['var'][:]['gene_ids']]); col=pd.Series(np.arange(len(ens)),index=ens)
gi=col.reindex(G).values; ng=len(ens); remap=-np.ones(ng,int); remap[gi]=np.arange(len(G))
c=pd.read_parquet('work/zeng_cells.parquet')
hi_t=c.groupby('sample').total.transform(lambda t:t.quantile(.995))
keep=((c.n_genes>=500)&(c.pct_mt<=20)&(c.total<=hi_t)).values; print('zeng keep',keep.sum(),'of',n)
parts=[]
for a in range(0,n,20000):
    b=min(n,a+20000); lo,h=ip[a],ip[b]
    m=sp.csr_matrix((X['data'][lo:h],X['indices'][lo:h],ip[a:b+1]-lo),shape=(b-a,ng))
    m=m[keep[a:b]]; t=c.total.values[a:b][keep[a:b]]
    m=sp.diags(1e4/t)@m[:,gi]; m.data=np.log1p(m.data).astype(np.float32); parts.append(m.astype(np.float32))
Z=sp.vstack(parts).tocsr(); sp.save_npz('work/zeng_sub.npz',Z); c[keep].reset_index(drop=True).to_parquet('work/zeng_cells_qc.parquet'); print('Z',Z.shape,Z.nnz)
