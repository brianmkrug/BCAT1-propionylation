# Pass Z1: stream Zeng CSR in chunks: per-cell QC, per-gene log-norm mean/var, BCAT1/2/MYC counts
import h5py, numpy as np, pandas as pd, scipy.sparse as sp
f=h5py.File('/home/claude/pg/zeng/human_raw.h5ad','r')
X=f['X']; ip=X['indptr'][:]; n=len(ip)-1
var=f['var'][:]; sym=np.array([s.decode() for s in var['index']]); ens=np.array([s.decode() for s in var['gene_ids']])
obs=f['obs'][:]; cats=[c.decode() for c in f['uns/week_stage_categories'][:]]
ng=len(sym); mt=np.array([s.startswith('MT-') for s in sym]); rib=np.array([s.startswith(('RPS','RPL')) for s in sym])
sel={g:np.where(sym==g)[0][0] for g in ['BCAT1','BCAT2','MYC']}
tot=np.zeros(n); ngen=np.zeros(n); mtc=np.zeros(n); gc={g:np.zeros(n) for g in sel}
s1=np.zeros(ng); s2=np.zeros(ng); nz=np.zeros(ng)
step=20000
for a in range(0,n,step):
    b=min(n,a+step); lo,hi=ip[a],ip[b]
    m=sp.csr_matrix((X['data'][lo:hi],X['indices'][lo:hi],ip[a:b+1]-lo),shape=(b-a,ng))
    t=np.asarray(m.sum(1)).ravel(); tot[a:b]=t; ngen[a:b]=np.diff(m.indptr); mtc[a:b]=np.add.reduceat(np.r_[m.data*mt[m.indices],0],m.indptr[:-1]) if False else np.asarray((m@sp.csr_matrix(mt.astype(float)).T).todense()).ravel()
    csc=m[:,list(sel.values())].toarray()
    for k,g in enumerate(sel): gc[g][a:b]=csc[:,k]
    m=sp.diags(1e4/np.maximum(t,1))@m; m.data=np.log1p(m.data)
    s1+=np.bincount(m.indices,weights=m.data,minlength=ng); s2+=np.bincount(m.indices,weights=m.data**2,minlength=ng); nz+=np.bincount(m.indices,minlength=ng)
    print(a,end=' ',flush=True)
cell=pd.DataFrame({'barcode':[s.decode() for s in obs['index']],'sample':[cats[i] for i in obs['week_stage']],'total':tot,'n_genes':ngen,'pct_mt':mtc/np.maximum(tot,1)*100,**gc})
cell.to_parquet('/home/claude/pg/work/zeng_cells.parquet')
mean=s1/n; v=s2/n-mean**2
pd.DataFrame({'symbol':sym,'ensembl':ens,'mean':mean,'var':v,'frac':nz/n}).to_parquet('/home/claude/pg/work/zeng_genes.parquet')
print('done')
