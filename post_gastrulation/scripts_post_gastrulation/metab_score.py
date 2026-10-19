import numpy as np,pandas as pd,json,scipy.sparse as sp
P=json.load(open('metab/metabolic_panel.json')); T=json.load(open('metab/targeted_modules.json'))
P['MCT/lactate (targeted)']=T['MCT_lactate_transport']; P['Propionyl-CoA input Ile/Val (targeted, -BCAT1)']=[g for g in T['Propionyl_CoA_input_IleVal'] if g!='BCAT1']; P['Propionyl-CoA node (targeted)']=T['Propionyl_CoA_node']
var=pd.read_csv('work/zeng_var.csv'); sym=var.symbol.values; s2i={s:i for i,s in enumerate(sym)}; e2i={e:i for i,e in enumerate(var.ensembl)}
x=pd.read_csv('xu/GSE157329_gene_annotate.txt.gz',sep='\t').set_index('gene_short_name').gene_id
alias={'ADSS1':'ADSSL1','ADSS2':'ADSS','EPRS1':'EPRS','KARS1':'KARS','RARS1':'RARS','WARS1':'WARS','H1-0':'H1F0','H2AZ1':'H2AFZ','H4C3':'HIST1H4C','MARCHF2':'MARCH2','MARCHF6':'MARCH6','MARCHF8':'MARCH8','MTARC2':'MARC2','GFUS':'TSTA3','G6PC1':'G6PC','DGKK':'DGKK','UQCRHL':'UQCRHL'}
mtenc={'ATP6','ATP8','COX1','COX2','COX3','CYTB','ND1','ND2','ND3','ND4','ND4L','ND5','ND6'}
def idx(g):
    if g in mtenc or g=='MCT1': return None
    if g in s2i: return s2i[g]
    if g in x.index and x[g] in e2i: return e2i[x[g]]
    if alias.get(g) in s2i: return s2i[alias[g]]
    return None
PI={k:sorted({i for i in (idx(g) for g in v) if i is not None}) for k,v in P.items()}
cov=pd.Series({k:f'{len(PI[k])}/{len(v)}' for k,v in P.items()})
G=pd.read_csv('work/brain_groups.csv').iloc[:,0].values; pb=np.load('work/brain_pb.npy'); pbs=np.load('work/brain_pb_bysample.npy')
ab=pd.read_parquet('work/brain_cells.parquet'); ncell=ab.grp.value_counts().reindex(G).values
coarse=ab.groupby('grp').marker_coarse.first().reindex(G); coarse[G=='OPC (PDGFRA+ OLIG2+)']='OPC'
def pscore(pbm):
    lc=np.log1p(pbm/pbm.sum(0,keepdims=True)*1e6)
    return pd.DataFrame({k:lc[v].mean(0) for k,v in PI.items()},index=G).T
S=pscore(pb)
OPC='OPC (PDGFRA+ OLIG2+)'
neural=[g for g,c,n in zip(G,coarse,ncell) if c in ('Neural progenitor','Neuron') and n>=300]
nonneu=[g for g,c,n in zip(G,coarse,ncell) if c in ('Endothelium','Immune','Mesenchyme') and n>=300]
cols=[OPC]+neural
Z=S[cols].sub(S[cols].mean(1),axis=0).div(S[cols].std(1),axis=0)
delta=S[OPC]-S[neural].mean(1)
# rank of OPC among neural groups
rank=S[cols].rank(axis=1,ascending=False)[OPC]
# per-sample deltas
ds={}
for k,s in enumerate(['W9-1','W9-2','W12-1']):
    Ss=pscore(pbs[k]+1e-9); nn=[g for g in neural if (ab[(ab['sample']==s)].grp==g).sum()>=50]
    ds[s]=Ss[OPC]-Ss[nn].mean(1)
res=pd.DataFrame({'genes_used':cov,'OPC_logCPM':S[OPC],'other_neural_mean':S[neural].mean(1),'delta_OPC_vs_neural':delta,'OPC_z_among_neural':Z[OPC],'OPC_rank_of_%d'%len(cols):rank,
 'delta_W9-1(n=43)':ds['W9-1'],'delta_W9-2(n=28)':ds['W9-2'],'delta_W12-1(n=440)':ds['W12-1']})
res['delta_vs_nonneural_residents']=S[OPC]-S[nonneu].mean(1)
res.to_csv('work/metab_group_results.csv'); S.to_csv('work/metab_group_logcpm.csv'); Z.to_csv('work/metab_group_z_neural.csv')
json.dump({'neural':neural,'nonneural':nonneu},open('work/metab_groups.json','w'))
pd.set_option('display.width',250); print(res.sort_values('delta_OPC_vs_neural',ascending=False).round(3).to_string())
