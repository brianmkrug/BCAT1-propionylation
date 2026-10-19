import h5py,numpy as np,pandas as pd,scipy.sparse as sp,json
P=json.load(open('metab/metabolic_panel.json')); T=json.load(open('metab/targeted_modules.json'))
f=h5py.File('zeng/human_raw.h5ad','r'); X=f['X']; ip=X['indptr'][:]; n=len(ip)-1
var=f['var'][:]; sym=np.array([s.decode() for s in var['index']]); ens=np.array([s.decode() for s in var['gene_ids']])
c=pd.read_parquet('work/zeng_cells.parquet'); hi=c.groupby('sample').total.transform(lambda t:t.quantile(.995))
keep=((c.n_genes>=500)&(c.pct_mt<=20)&(c.total<=hi)).values; orig=np.where(keep)[0]
a=pd.read_parquet('work/zeng_annotated_cells.parquet'); assert len(a)==len(orig)
brain=(a.tissue=='brain').values; bidx=orig[brain]; ab=a[brain].reset_index(drop=True)
ab['grp']=np.where(ab.OPC_like,'OPC (PDGFRA+ OLIG2+)',ab.marker_fine)
ab.to_parquet('work/brain_cells.parquet')
G=sorted(ab.grp.unique()); gi=pd.Series(range(len(G)),index=G); code=gi[ab.grp].values
ng=len(sym); pb=np.zeros((ng,len(G))); pbs={s:np.zeros((ng,len(G))) for s in ['W9-1','W9-2','W12-1']}
parts=[]
order=np.argsort(bidx); bs=bidx[order]  # sorted
rows=[]
for st in range(0,len(bs),5000):
    sel=bs[st:st+5000]
    lo,hi_=ip[sel[0]],ip[sel[-1]+1]
    m=sp.csr_matrix((X['data'][lo:hi_],X['indices'][lo:hi_],ip[sel[0]:sel[-1]+2]-lo),shape=(sel[-1]-sel[0]+1,ng))[sel-sel[0]]
    oi=order[st:st+5000]; cc=code[oi]
    I=sp.csr_matrix((np.ones(len(cc)),(cc,np.arange(len(cc)))),shape=(len(G),len(cc)))
    pb+=(I@m).T.toarray()
    for s in pbs:
        k=(ab['sample'].values[oi]==s)
        if k.any(): pbs[s]+=(I[:,k]@m[k]).T.toarray()
    parts.append(m.astype(np.float32)); rows.append(oi)
M=sp.vstack(parts).tocsr(); ro=np.concatenate(rows); inv=np.argsort(ro); M=M[inv]  # back to ab order
sp.save_npz('work/brain_counts_all.npz',M)
np.save('work/brain_pb.npy',pb); np.save('work/brain_pb_bysample.npy',np.stack([pbs[s] for s in ['W9-1','W9-2','W12-1']]))
pd.Series(G).to_csv('work/brain_groups.csv',index=False); pd.DataFrame({'symbol':sym,'ensembl':ens}).to_csv('work/zeng_var.csv',index=False)
print('brain cells',len(ab),'groups',len(G),'nnz',M.nnz)
