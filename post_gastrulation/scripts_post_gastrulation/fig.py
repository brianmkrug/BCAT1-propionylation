import pandas as pd,numpy as np,matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
d=pd.read_parquet('work/zeng_annotated_cells.parquet'); d=d[~d.marker_coarse.isin(['Doublet/mixed'])]
g=d.groupby(['marker_coarse','week']); n=g.size(); cpm=(g.BCAT1.sum()/g.total.sum()*1e6).where(n>=50).unstack()
order=cpm.mean(1).sort_values(ascending=False).index; cpm=cpm.loc[order]
b=d[(d.tissue=='brain')&d.marker_coarse.isin(['Neural progenitor','Neuron'])].copy(); b['grp']=np.where(b.OPC_like,'OPC (PDGFRA+ OLIG2+)',b.marker_fine.str.replace(r'Cortical neurons \+ RG mix','Cortical neurons + RG',regex=True).str.replace(r'Late radial glia / glial progenitor','Late radial glia',regex=True).str.replace(r'Dorsal neural tube / roof plate','Roof plate',regex=True).str.replace(r', [^)]*\)',')',regex=True))
gb=b.groupby('grp'); bb=pd.DataFrame({'n':gb.size(),'cpm':gb.BCAT1.sum()/gb.total.sum()*1e6}); bb=bb[(bb.n>=300)&~bb.index.isin(['Erythroid','Endothelium','Immune','Mesenchyme / fibroblast','Mixed neural-mesenchyme'])].sort_values('cpm')
fig,ax=plt.subplots(1,2,figsize=(10.5,4.3),gridspec_kw={'width_ratios':[1.15,1],'wspace':1.25})
cm=LinearSegmentedColormap.from_list('g',['#f2f7f5','#1b9e77','#0b4a37']); cm.set_bad('#ffffff')
im=ax[0].imshow(cpm.values,cmap=cm,vmin=0,vmax=150,aspect='auto')
ax[0].set_xticks(range(cpm.shape[1]),[f'W{w}' for w in cpm.columns]); ax[0].set_yticks(range(len(order)),order)
for i in range(cpm.shape[0]):
    for j in range(cpm.shape[1]):
        v=cpm.values[i,j]
        if np.isfinite(v): ax[0].text(j,i,f'{v:.0f}',ha='center',va='center',fontsize=6.5,color='white' if v>90 else '#222')
ax[0].set_xlabel('Post-conceptional week  (W3–W8 whole body/head; W9–W12 brain only)'); ax[0].set_title('BCAT1 by lineage and week (Zeng 2023)',loc='left',fontsize=9)
for s in ax[0].spines.values(): s.set_visible(False)
ax[0].tick_params(length=0)
cb=fig.colorbar(im,ax=ax[0],fraction=0.035,pad=0.03); cb.set_label('BCAT1 pseudobulk CPM'); cb.outline.set_visible(False)

col=['#6c3483' if 'OPC' in i else '#1b9e77' for i in bb.index]
ax[1].barh(range(len(bb)),bb.cpm,color=col,height=0.7)
ax[1].set_yticks(range(len(bb)),[f'{i} (n={n:,})' for i,n in zip(bb.index,bb.n)],fontsize=7)
ax[1].set_xlabel('BCAT1 pseudobulk CPM'); ax[1].set_title('Brain, weeks 9–12: neural cell types',loc='left',fontsize=9)
ax[1].xaxis.grid(True,color='#e6e6e6',lw=0.6); ax[1].set_axisbelow(True); ax[1].tick_params(axis='y',length=0)
for i,v in enumerate(bb.cpm): ax[1].text(v+1.5,i,f'{v:.0f}',va='center',fontsize=6.5,color='#333')
fig.text(0.01,0.005,'Zeng et al. 2023 (GSE155121), 446,689 cells after QC; labels from marker-based clustering, OPC = PDGFRA+OLIG2+ cells within neural-progenitor clusters. Pseudobulk CPM = summed BCAT1 UMIs / summed UMIs ×1e6; groups <50 cells blank.',fontsize=5.8,color='#555')
fig.savefig('out/zeng2023_BCAT1_overview.png',dpi=220,bbox_inches='tight'); fig.savefig('out/zeng2023_BCAT1_overview.pdf',bbox_inches='tight')
