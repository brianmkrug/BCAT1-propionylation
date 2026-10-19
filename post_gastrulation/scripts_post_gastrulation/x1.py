# Pass X1: Xu pseudobulk by final_annotation + per-gene log-norm mean/var
import pandas as pd, numpy as np
c=pd.read_csv('/home/claude/pg/xu_cells_BCAT.csv.gz'); tot=c.tot.values
lab=(c['developmental system']+'|'+c['final_annotation']).astype('category'); lc=lab.cat.codes.values; L=len(lab.cat.categories)
ng=32351; pb=np.zeros((ng,L)); s1=np.zeros(ng); s2=np.zeros(ng); nzc=np.zeros((ng,L))
for ch in pd.read_csv('/home/claude/pg/xu/GSE157329_raw_counts.mtx.gz',sep=' ',comment='%',skiprows=2,header=None,names=['g','c','v'],dtype={'g':np.int32,'c':np.int32,'v':np.float32},chunksize=20_000_000):
    g=ch.g.values-1; ci=ch.c.values-1; v=ch.v.values
    idx=g*L+lc[ci]
    pb+=np.bincount(idx,weights=v,minlength=ng*L).reshape(ng,L)
    nzc+=np.bincount(idx,minlength=ng*L).reshape(ng,L)
    ln=np.log1p(v*1e4/tot[ci]); s1+=np.bincount(g,weights=ln,minlength=ng); s2+=np.bincount(g,weights=ln**2,minlength=ng)
    print('.',end='',flush=True)
n=len(c); np.save('/home/claude/pg/work/xu_pb_type.npy',pb); np.save('/home/claude/pg/work/xu_nz_type.npy',nzc)
pd.Series(lab.cat.categories).to_csv('/home/claude/pg/work/xu_types.csv',index=False)
g=pd.read_csv('/home/claude/pg/xu/GSE157329_gene_annotate.txt.gz',sep='\t')
g['mean']=s1/n; g['var']=s2/n-g['mean']**2; g.to_parquet('/home/claude/pg/work/xu_genes.parquet'); print('done')
