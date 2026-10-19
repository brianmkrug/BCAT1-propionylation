import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt, matplotlib.cm as cm
plt.rcParams['svg.fonttype']='none'; plt.rcParams['pdf.fonttype']=42
__ORIG_SF=plt.savefig
def __sf_all(fname,*a,**k):
    import os as _os
    b,_=_os.path.splitext(str(fname)); k.pop('format',None); k.setdefault('bbox_inches','tight')
    __ORIG_SF(b+'.png',*a,**dict(k))
    __ORIG_SF(b+'.svg',*a,**dict(k,format='svg'))
    __ORIG_SF(b+'.pdf',*a,**dict(k,format='pdf'))
plt.savefig=__sf_all
from matplotlib.patches import Rectangle
J=REPO_ROOT + "/jessa2022_metabolic_ssgsea"
alldf=pd.read_csv(f"{J}/dotplot_panel_top20lipid.csv")
glord=['DMG-H3.1_2K27M-pons','DMG-H3.3K27M-pons','DMG-H3.3K27M-thalamus','HGG-H3.3G34-cortex','HGG-H3WT-cortex','Ependymoma-PFA']
order=glord+['MB-WNT','MB-SHH','MB-Group3','MB-Group4']
short=['DMG H3.1/2K27M pons','DMG H3.3K27M pons','DMG H3.3K27M thal','HGG G34 cortex','HGG H3WT cortex','Epend PFA','MB WNT','MB SHH','MB Group3','MB Group4']
groups=[
        ('Fatty acid','#fdae6b',['SCD5','ACSL4']),
        ('Glycero-/phospholipid','#a1d99b',['DGKI','DGKG','AGPAT4']),
        ('Sphingolipid','#9ecae1',['UGT8','NEU4']),
        ('Glycosphingolipid','#bcbddc',['ST8SIA1','ST3GAL5','ST6GALNAC3','B3GALT5']),
        ('Eicosanoid /\narachidonic','#fc9272',['PTGDS','PLA2G4A']),
        ('Sterol / bile acid','#fee391',['FDFT1','DHCR24']),
        ('GPI','#c9b38c',['PIGM'])]
gorder=[g for _,_,gg in groups for g in gg]
tgroups=[('H3K27M-DMG','#e0322f',[0,1,2]),('G34R/V','#2ca02c',[3]),('H3WT','#1f77b4',[4]),
         ('EPN','#8e44ad',[5]),('Medulloblastoma','#ff7f0e',[6,7,8,9])]
xi={s:i for i,s in enumerate(order)}; yi={g:i for i,g in enumerate(gorder)}
vmax=np.nanmax(alldf['mean_lognorm'].values); cmap=plt.cm.Reds
GENE_FS=13.0
fig,ax=plt.subplots(figsize=(15,13.5))
for _,r in alldf.iterrows():
    if r['gene'] not in yi: continue
    x=xi[r['subgroup']]; y=yi[r['gene']]
    if np.isnan(r['pct_expressing']): ax.scatter(x,y,s=42,facecolors='none',edgecolors='0.55',lw=0.9,marker='x'); continue
    ax.scatter(x,y,s=14+r['pct_expressing']*14,c=[cmap(r['mean_lognorm']/vmax)],edgecolors='k',lw=0.5,zorder=3)
ax.set_yticks([])
GENE_X=-0.65
texts=[ax.text(GENE_X,i,g,ha='right',va='center',fontsize=GENE_FS) for g,i in yi.items()]
ax.set_xticks(range(len(order))); ax.set_xticklabels(short,rotation=45,ha='right',fontsize=12)
ax.invert_yaxis()
BW=0.42   # pathway box width
LABELSPACE=2.0
GLAB={'Fatty acid':'Fatty\nacid','Glycero-/phospholipid':'Glycero-/\nphospholipid','Sphingolipid':'Sphingo-\nlipid','Glycosphingolipid':'Glyco-\nsphingolipid','Eicosanoid /\narachidonic':'Eicosanoid\n/ arachid.','Sterol / bile acid':'Sterol /\nbile acid','GPI':'GPI'}
left=-3.6
for _ in range(4):   # iterate: measure gene-label extent, tuck boxes right beside them
    ax.set_xlim(left,len(order)-0.4); ax.set_ylim(len(gorder)-0.3,-2.35)
    fig.canvas.draw()
    inv=ax.transData.inverted()
    minx=min(inv.transform((t.get_window_extent().x0,0))[0] for t in texts)
    box_r=minx-0.10; box_l=box_r-BW; left=box_l-0.12
labx=box_l-0.18; left=labx-LABELSPACE
ax.set_xlim(left,len(order)-0.4); ax.set_ylim(len(gorder)-0.3,-2.35)
for name,col,gg in groups:
    ys=[yi[g] for g in gg]; y0,y1=min(ys)-0.45,max(ys)+0.45
    ax.add_patch(Rectangle((box_l,y0),BW,y1-y0,facecolor=col,edgecolor='k',lw=0.7,clip_on=False,zorder=2))
    ax.text(labx,(y0+y1)/2,GLAB.get(name,name),ha='right',va='center',rotation=0,fontsize=13,zorder=3)
    ax.axhline(y1,color='0.85',lw=0.6,zorder=0)
# tumor color bar + color-coded subgroup names above it (no legend)
bb,bt=-1.30,-1.75
for name,col,cols in tgroups:
    x0=min(cols)-0.5+0.1; x1=max(cols)+0.5-0.1
    ax.add_patch(Rectangle((x0,bt),x1-x0,bb-bt,facecolor=col,edgecolor='k',lw=0.6,clip_on=False,zorder=4))
    ax.text((x0+x1)/2,bt-0.16,name,ha='center',va='bottom',fontsize=14,fontweight='bold',color=col,zorder=5)
sm=cm.ScalarMappable(cmap=cmap,norm=plt.Normalize(0,vmax)); sm.set_array([])
cb=fig.colorbar(sm,ax=ax,shrink=0.42,pad=0.01); cb.set_label('mean log-norm expr',fontsize=13); cb.ax.tick_params(labelsize=11)
sh=[plt.scatter([],[],s=14+p*14,c='0.5',edgecolors='k',lw=0.5,label=f'{p}%') for p in [5,25,50]]
ax.legend(handles=sh,title='% expressing',loc='upper left',bbox_to_anchor=(1.10,1.0),fontsize=12,title_fontsize=13,frameon=False)
for sp in ax.spines.values(): sp.set_visible(False)
ax.tick_params(axis='both',length=0)
plt.tight_layout()
plt.savefig(f"{J}/dotplot_panel_top20lipid_grouped_v2.png",dpi=150,bbox_inches='tight')
plt.savefig(f"{J}/dotplot_panel_top20lipid_grouped_v2.svg",format='svg',bbox_inches='tight')
print("saved revised panel")
