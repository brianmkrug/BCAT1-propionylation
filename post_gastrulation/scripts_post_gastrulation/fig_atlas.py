import pandas as pd,numpy as np,matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':7.5,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
C1,C2,COPC='#1b9e77','#7570b3','#6c3483'
d=pd.read_parquet('work/zeng_annotated_cells.parquet'); d=d[d.marker_coarse!='Doublet/mixed'].copy()
OPC='OPC (PDGFRA+ OLIG2+)'
d['grp']=np.where(d.OPC_like&(d.tissue=='brain'),OPC,d.marker_fine)
d.loc[d.grp==OPC,'marker_coarse']='Oligodendrocyte lineage'
d['grp']=d.grp.str.replace(r'Cortical neurons \+ RG mix','Cortical neurons + RG',regex=True).str.replace('Late radial glia / glial progenitor','Late radial glia').str.replace('Dorsal neural tube / roof plate','Roof plate').str.replace('Immune (macrophage/microglia) + endothelial','Microglia / macrophage (+endoth.)',regex=False).str.replace(r', [^)]*\)',')',regex=True)
g=d.groupby(['grp','week']); n=g.size()
B={k:(g[k].sum()/g.total.sum()*1e6).where(n>=50).unstack() for k in ['BCAT1','BCAT2']}
lin=d.groupby('grp').marker_coarse.first()
lorder=['Oligodendrocyte lineage','Neural progenitor','Neuron','Neural crest / Schwann','Eye','Epithelium','Early mesoderm (W3)','Lateral plate mesoderm','Somite / muscle','Mesenchyme','Endothelium','Immune','Erythroid']
rows=[]
for L in lorder:
    r=B['BCAT1'].loc[lin[lin==L].index].mean(axis=1).sort_values(ascending=False).index.tolist(); rows+=r
rows=[r for r in rows if B['BCAT1'].loc[r].notna().any()]
# panel A data: W9-12 all cell types
late=d[d.week>=9]; gl=late.groupby('grp'); bb=pd.DataFrame({'n':gl.size(),'B1':gl.BCAT1.sum()/gl.total.sum()*1e6,'B2':gl.BCAT2.sum()/gl.total.sum()*1e6})
bb=bb[bb.n>=100].sort_values('B1'); bb['lin']=lin.reindex(bb.index)
fig=plt.figure(figsize=(17,9.2)); gs0=fig.add_gridspec(1,4,width_ratios=[1.0,0.62,1,1],wspace=0.06); gs=type('G',(),{})(); gs.__getitem__=None
class _G:
    def __getitem__(self,k): return gs0[0,{0:0,1:2,2:3}[k[1]]]
gs=_G()
ax=fig.add_subplot(gs[0,0]); y=np.arange(len(bb)); h=0.4
ax.barh(y+h/2,bb.B1,h,color=C1,label='BCAT1'); ax.barh(y-h/2,bb.B2,h,color=C2,label='BCAT2')
ax.set_yticks(y,[f'{i}  (n={k:,})' for i,k in zip(bb.index,bb.n)])
nonneural=~bb.lin.isin(['Neural progenitor','Neuron','Oligodendrocyte lineage'])
for t,nn in zip(ax.get_yticklabels(),nonneural):
    if OPC in t.get_text(): t.set_fontweight('bold'); t.set_color(COPC)
    elif nn: t.set_color('#8a6d3b'); t.set_style('italic')
for i,(v1,v2) in enumerate(zip(bb.B1,bb.B2)): ax.text(max(v1,v2)+1.5,i,f'{v1:.0f} / {v2:.0f}',va='center',fontsize=6.3,color='#444')
ax.set_xlabel('pseudobulk CPM'); ax.set_title('A  W9–12 brain dissections: all cell types\n(italic = non-neural residents; label = BCAT1 / BCAT2)',loc='left',fontsize=8.5)
ax.xaxis.grid(True,color='#e6e6e6',lw=0.6); ax.set_axisbelow(True); ax.tick_params(axis='y',length=0); ax.legend(frameon=False,loc='center right'); ax.set_xlim(0,112)
for k,(key,col,title) in enumerate([('BCAT1',C1,'B  BCAT1 across the full atlas (cell type × week)'),('BCAT2',C2,'C  BCAT2 across the full atlas')]):
    sg=gs[0,k+1].subgridspec(2,1,height_ratios=[1,0.025],hspace=0.06); a=fig.add_subplot(sg[0])
    cm=LinearSegmentedColormap.from_list(key,['#f6f7f7',col,'#1d1d3a' if key=='BCAT2' else '#0b4a37']); cm.set_bad('white')
    M=B[key].loc[rows]; im=a.imshow(M.values,cmap=cm,vmin=0,vmax=150,aspect='auto')
    a.set_xticks(range(M.shape[1]),[f'W{w}' for w in M.columns]); a.xaxis.tick_top()
    a.set_yticks(range(len(rows)),rows if k==0 else []); a.tick_params(length=0)
    if k==0:
        for t in a.get_yticklabels():
            if OPC in t.get_text(): t.set_fontweight('bold'); t.set_color(COPC)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v=M.values[i,j]
            if np.isfinite(v): a.text(j,i,f'{v:.0f}',ha='center',va='center',fontsize=5.6,color='white' if v>90 else '#222')
    i0=rows.index(OPC); a.add_patch(Rectangle((-0.5,i0-0.5),M.shape[1],1,fill=False,ec=COPC,lw=1.4))
    # lineage separators
    prev=None
    for i,r in enumerate(rows):
        if prev is not None and lin[r]!=prev: a.axhline(i-0.5,color='white',lw=2.2)
        prev=lin[r]
    a.axvline(5.5,color='#333',lw=0.8,ls=(0,(2,2)))
    for s in a.spines.values(): s.set_visible(False)
    a.set_title(title,loc='left',fontsize=8.5,pad=18)
    cb=fig.colorbar(im,cax=fig.add_subplot(sg[1]),orientation='horizontal'); cb.set_label(f'{key} pseudobulk CPM (0–150, same scale)'); cb.outline.set_visible(False)
fig.text(0.01,0.005,'Zeng et al. 2023 (GSE155121). W3–W8 = whole body/head dissections; W9–W12 (right of dashed line) = brain dissections only, so non-brain organs are absent there. '
 'Labels: marker-based clustering; OPC = PDGFRA+OLIG2+ cells within neural-progenitor clusters (brain samples). Cells: groups <50 cells blank (heatmaps), <100 cells omitted (bars).',fontsize=6.2,color='#555')
fig.savefig('out/zeng2023_BCAT_atlas_OPC.png',dpi=200,bbox_inches='tight'); fig.savefig('out/zeng2023_BCAT_atlas_OPC.pdf',bbox_inches='tight')
