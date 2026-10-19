import numpy as np,pandas as pd,json
exec(open('scripts/metab_score.py').read().split("S=pscore(pb)")[0])  # reuse PI, G, pb, ab, coarse
gr=json.load(open('work/metab_groups.json')); OPC='OPC (PDGFRA+ OLIG2+)'; neural=gr['neural']
lc=pd.DataFrame(np.log1p(pb/pb.sum(0,keepdims=True)*1e6),columns=G)
opc=lc[OPC].values; oth=lc[neural].mean(1).values; d=opc-oth; avg=(opc+oth)/2
expr=avg>np.log1p(1)  # CPM>1 averaged
print('genome-wide median delta (expressed genes)',np.median(d[expr]).round(3),' mean',d[expr].mean().round(3))
bins=pd.qcut(avg[expr],20,labels=False); eidx=np.where(expr)[0]; bin_of=pd.Series(bins,index=eidx)
rng=np.random.default_rng(0); out=[]
for k,v in PI.items():
    v=[i for i in v if expr[i]]
    if len(v)<3: continue
    obs=d[v].mean(); vb=bin_of[v].values
    null=np.array([d[[rng.choice(eidx[bins==b]) for b in vb]].mean() for _ in range(2000)])
    out.append((k,len(v),obs,null.mean(),(obs-null.mean())/null.std(),(np.sum(null>=obs)+1)/2001,(np.sum(null<=obs)+1)/2001))
R=pd.DataFrame(out,columns=['pathway','n_expr_genes','delta','null_mean','z_vs_matched','p_up','p_down']).set_index('pathway')
from statsmodels.stats.multitest import multipletests
R['p2']=np.minimum(1,2*R[['p_up','p_down']].min(1)); R['fdr']=multipletests(R.p2,method='fdr_bh')[1]
R['excess_delta']=R.delta-R.null_mean
R.sort_values('z_vs_matched',ascending=False).round(3).to_csv('work/metab_null_results.csv')
pd.set_option('display.width',200); print(R.sort_values('z_vs_matched',ascending=False).round(3).to_string())
