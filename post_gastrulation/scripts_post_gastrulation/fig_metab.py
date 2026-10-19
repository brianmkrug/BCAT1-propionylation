import numpy as np,pandas as pd,json,matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':7.5,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
exec(open('scripts/metab_score.py').read().split("S=pscore(pb)")[0])
gr=json.load(open('work/metab_groups.json')); OPC='OPC (PDGFRA+ OLIG2+)'; COPC='#6c3483'
D=pd.read_csv('work/metab_consensus.csv',index_col=0); S=pd.read_csv('work/metab_group_logcpm.csv',index_col=0)
import re
short=lambda s: 'OPC' if s==OPC else re.sub(r', [^)]*\)',')',s.replace('Cortical neurons + RG mix','Cortical neurons + RG').replace('Late radial glia / glial progenitor','Late radial glia').replace('Dorsal neural tube / roof plate','Roof plate'))
prog=[g for g in gr['neural'] if ab[ab.grp==g].marker_coarse.iloc[0]=='Neural progenitor']; neu=[g for g in gr['neural'] if g not in prog]
cols=[OPC]+sorted(prog,key=lambda g:-S.loc['Sphingolipid',g])+sorted(neu,key=lambda g:-S.loc['Sphingolipid',g])
up=['Sphingolipid','Unsat FA synthesis','Glycerophospholipid','Amino/nucleotide sugar','Lactate/MCT transport','MCT/lactate (targeted)','NAD/nicotinamide','Purine metabolism','Pyruvate metabolism','Glycolysis','HALLMARK mTORC1','Glucose transport']
mid=['Propionyl-CoA node (targeted)','BCAA degradation','One-carbon/folate','Pyrimidine metabolism','Mitochondrial biogenesis']
dn=['Inositol phosphate','Butanoate metab','Ala/Asp/Glu metab']
rows=up+mid+dn
Z=S.loc[rows,cols]; Z=Z.sub(Z.mean(axis=1),axis=0).div(Z.std(axis=1),axis=0)
# genes
genes={'Galactosphingolipid / myelin lipid':['UGT8','GAL3ST1','ASAH1','SGMS1','SCD5','ELOVL5','ELOVL1','PPAPDC1A','PLA2G4A','LPCAT2'],
 'UDP-sugar / glycogen':['PGM1','PGM2','UGP2','UGDH'],'Lactate / monocarboxylate':['SLC16A3','SLC16A7','SLC16A1','LDHA','LDHB'],
 'NAD / purinergic / cAMP':['NT5E','SIRT2','PNP','NAPRT1','PDE7B','PDE4B'],'Inositol':['SLC2A13','ISYNA1'],
 'Propionyl-CoA node':['PCCB','SUCLG1','SUCLG2'],'Proliferative 1C / nucleotide':['TYMS','DHFR','MTHFD1','RRM2'],'BCAT':['BCAT1','BCAT2']}
gl=[g for v in genes.values() for g in v]; gidx=[s2i[g] for g in gl]
cpm=pb/pb.sum(0,keepdims=True)*1e6; lc=pd.DataFrame(np.log1p(cpm[gidx]),index=gl,columns=G)[cols]
GZ=lc.sub(lc.mean(axis=1),axis=0).div(lc.std(axis=1),axis=0)
fig=plt.figure(figsize=(16,10)); gs=fig.add_gridspec(1,3,width_ratios=[1.25,0.55,1.15],wspace=0.55)
def heat(ax,M,title,ylabels):
    im=ax.imshow(M.values,cmap='RdBu_r',vmin=-2,vmax=2,aspect='auto')
    ax.set_xticks(range(M.shape[1]),[short(c) for c in M.columns],rotation=60,ha='right'); ax.set_yticks(range(M.shape[0]),ylabels)
    ax.get_xticklabels()[0].set_color(COPC); ax.get_xticklabels()[0].set_fontweight('bold')
    ax.axvline(0.5,color='white',lw=2.6); ax.axvline(len(prog)+0.5,color='white',lw=2.6)
    for s in ax.spines.values(): s.set_visible(False)
    ax.tick_params(length=0); ax.set_title(title,loc='left',fontsize=9); return im
ax=fig.add_subplot(gs[0,0]); im=heat(ax,Z,'A  Metabolic pathway scores, W9–12 brain neural cell types\n    (mean log1p pseudobulk CPM, row z-score; Velmeshev-matched panel)\n    columns: OPC | progenitors | neurons',rows)
for y in [len(up)-0.5,len(up)+len(mid)-0.5]: ax.axhline(y,color='white',lw=2.6)

cb=fig.colorbar(im,ax=ax,fraction=0.025,pad=0.02); cb.set_label('Row z-score'); cb.outline.set_visible(False)
# B: evidence
ax=fig.add_subplot(gs[0,1]); y=np.arange(len(rows))
ax.barh(y-0.2,D.loc[rows,'pb_excess'],0.38,color=COPC,label='Zeng W9–12: OPC − neural\n(offset-corrected Δlog-CPM)')
ax.barh(y+0.2,D.loc[rows,'vel_OPC_vs_progenitors'].fillna(0),0.38,color='#b39cc8',label='Velmeshev: OPC − progenitors\n(median-centred Δlog-CPM)')
ax.set_yticks(y,['']*len(rows)); ax.set_ylim(len(rows)-0.5,-0.5); ax.axvline(0,color='#333',lw=0.6)
for i,r in enumerate(rows):
    c=D.loc[r,'call']; m='●' if c=='OPC-high (robust)' else ('▼' if c=='OPC-low (robust)' else '')
    ax.text(ax.get_xlim()[1] if False else 0.95,i,m,va='center',fontsize=7,color=COPC)
ax.set_xlim(-0.4,1.0); ax.set_title('B  OPC-distinct Δ\n    ● robust up  ▼ robust down',loc='left',fontsize=9)
ax.legend(frameon=False,fontsize=6.3,loc='upper left',bbox_to_anchor=(-0.05,-0.02)); ax.xaxis.grid(True,color='#eee'); ax.set_axisbelow(True)
for yy in [len(up)-0.5,len(up)+len(mid)-0.5]: ax.axhline(yy,color='#ccc',lw=0.8)
# C: genes
ax=fig.add_subplot(gs[0,2]); im2=heat(ax,GZ,'C  Driver genes (log1p pseudobulk CPM, row z-score)',gl)
k=0
for name,v in genes.items():
    k+=len(v); ax.axhline(k-0.5,color='white',lw=2.6)
    ax.text(GZ.shape[1]-0.3,k-len(v)/2-0.5,' '+name,va='center',fontsize=6.8,color='#333')
fig.text(0.01,-0.14,'Zeng et al. 2023 brain samples W9-1, W9-2, W12-1 (OPC n=511: 43/28/440). Robust = matched-random-gene-set test on pseudobulk (z≥1.5) AND control-gene module score higher in OPC vs all neural cells and vs progenitors '
 '(OLS adjusted for S/G2M score, depth, sample; t≥5). Panel B Velmeshev values from human_data/fine_pathway_logcpm.csv (OPC − mean of ExN/IN/Glial progenitor). Glucose transport is driven by SLC2A13 (myo-inositol transporter).',fontsize=6,color='#555',wrap=True)
fig.savefig('out/zeng2023_OPC_metabolic_signatures.png',dpi=200,bbox_inches='tight'); fig.savefig('out/zeng2023_OPC_metabolic_signatures.pdf',bbox_inches='tight')
