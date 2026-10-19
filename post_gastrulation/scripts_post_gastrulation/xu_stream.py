import pandas as pd, numpy as np, sys
c=pd.read_csv('GSE157329_cell_annotate.txt.gz',sep='\t')
sysl=c['developmental system'].astype('category'); sc=sysl.cat.codes.values; S=len(sysl.cat.categories)
ng,nc=32351,185140
tot=np.zeros(nc); pb=np.zeros((ng,S)); sel={15025:"MYC",19811:"BCAT1",30463:"BCAT2"}; gv={k:np.zeros(nc) for k in sel}
for ch in pd.read_csv('GSE157329_raw_counts.mtx.gz',sep=' ',comment='%',skiprows=2,header=None,names=['g','c','v'],dtype={'g':np.int32,'c':np.int32,'v':np.float32},chunksize=20_000_000):
    g=ch.g.values-1; ci=ch.c.values-1; v=ch.v.values
    np.add.at(tot,ci,v) if False else None
    tot+=np.bincount(ci,weights=v,minlength=nc)
    np.add.at(pb,(g,sc[ci]),v)
    for k in sel:
        m=g==k
        if m.any(): gv[k]+=np.bincount(ci[m],weights=v[m],minlength=nc)
    print('.',end='',flush=True)
np.save('/home/claude/pg/xu_pb_system.npy',pb)
out=pd.DataFrame({'tot':tot,**{sel[k]:gv[k] for k in sel}}); out=pd.concat([c,out],axis=1); out.to_csv('/home/claude/pg/xu_cells_BCAT.csv.gz',index=False)
print('\ndone', np.allclose(tot,c.total_UMIs))
