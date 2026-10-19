import pandas as pd,numpy as np,matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
C1,C2='#1b9e77','#7570b3'
d=pd.read_parquet('work/zeng_annotated_cells.parquet'); d=d[d.marker_coarse!='Doublet/mixed']
g=d.groupby(['marker_coarse','week']); n=g.size()
m={k:(g[k].sum()/g.total.sum()*1e6).where(n>=50).unstack() for k in ['BCAT1','BCAT2']}
order=m['BCAT1'].mean(axis=1).sort_values(ascending=False).index
b=d[(d.tissue=='brain')&d.marker_coarse.isin(['Neural progenitor','Neuron'])].copy()
b['grp']=np.where(b.OPC_like,'OPC (PDGFRA+ OLIG2+)',b.marker_fine.str.replace(r'Cortical neurons \+ RG mix','Cortical neurons + RG',regex=True).str.replace('Late radial glia / glial progenitor','Late radial glia').str.replace('Dorsal neural tube / roof plate','Roof plate').str.replace(r', [^)]*\)',')',regex=True))
gb=b.groupby('grp'); bb=pd.DataFrame({'n':gb.size(),'B1':gb.BCAT1.sum()/gb.total.sum()*1e6,'B2':gb.BCAT2.sum()/gb.total.sum()*1e6}); bb=bb[bb.n>=300].sort_values('B1')
fig=plt.figure(figsize=(14,5.0)); gs=fig.add_gridspec(1,3,width_ratios=[1,1,1.25],wspace=0.15)
axs=[fig.add_subplot(gs[0,0]),fig.add_subplot(gs[0,1])]; gs2=gs[0,2].subgridspec(1,2,width_ratios=[1.1,1]); axs+= [None,fig.add_subplot(gs2[0,1])]
cmaps={'BCAT1':LinearSegmentedColormap.from_list('a',['#f4f8f6',C1,'#0b4a37']),'BCAT2':LinearSegmentedColormap.from_list('b',['#f5f5fa',C2,'#2e2a5c'])}
for k,(key,ax) in enumerate(zip(['BCAT1','BCAT2'],axs[:2])):
    M=m[key].reindex(order); im=ax.imshow(M.values,cmap=cmaps[key],vmin=0,vmax=150,aspect='auto')
    ax.set_xticks(range(M.shape[1]),[f'W{w}' for w in M.columns]); ax.set_yticks(range(len(order)),order if k==0 else [])
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v=M.values[i,j]
            if np.isfinite(v): ax.text(j,i,f'{v:.0f}',ha='center',va='center',fontsize=6.5,color='white' if v>90 else '#222')
    ax.set_title(f'{key} by lineage and week',loc='left',fontsize=9,color=C1 if key=='BCAT1' else C2)
    for s in ax.spines.values(): s.set_visible(False)
    ax.tick_params(length=0)
    cb=fig.colorbar(im,ax=ax,orientation='horizontal',fraction=0.035,pad=0.09); cb.set_label(f'{key} pseudobulk CPM'); cb.outline.set_visible(False)
fig.text(0.13,0.96,'Columns: post-conceptional week (W3–W8 whole body/head; W9–W12 brain only). Same colour scale (0–150) in both heatmaps.',fontsize=7.5,color='#444')
ax=axs[3]; y=np.arange(len(bb)); h=0.38
ax.barh(y+h/2,bb.B1,h,color=C1,label='BCAT1'); ax.barh(y-h/2,bb.B2,h,color=C2,label='BCAT2')
ax.set_yticks(y,[f'{i} (n={n:,})' for i,n in zip(bb.index,bb.n)],fontsize=7)
for t in ax.get_yticklabels():
    if 'OPC' in t.get_text(): t.set_fontweight('bold')
for i,(v1,v2) in enumerate(zip(bb.B1,bb.B2)): ax.text(max(v1,v2)+1.5,i,f'{v1/v2:.1f}×',va='center',fontsize=6.5,color='#444')
ax.set_xlabel('pseudobulk CPM'); ax.set_title('Brain W9–12 neural cell types \n(label = BCAT1/BCAT2 ratio)',loc='left',fontsize=9)
ax.xaxis.grid(True,color='#e6e6e6',lw=0.6); ax.set_axisbelow(True); ax.tick_params(axis='y',length=0); ax.legend(frameon=False,loc='lower right')
fig.text(0.01,-0.04,'Zeng et al. 2023 (GSE155121), 446,689 cells after QC, marker-based labels; OPC = PDGFRA+OLIG2+ cells in neural-progenitor clusters. Pseudobulk CPM = summed gene UMIs / summed UMIs ×1e6; groups <50 cells blank.',fontsize=6,color='#555')
fig.savefig('out/zeng2023_BCAT1_vs_BCAT2.png',dpi=200,bbox_inches='tight'); fig.savefig('out/zeng2023_BCAT1_vs_BCAT2.pdf',bbox_inches='tight')
