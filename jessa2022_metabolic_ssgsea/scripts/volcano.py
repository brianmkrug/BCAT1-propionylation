import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import pandas as pd, numpy as np, json, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams['svg.fonttype']='none'; plt.rcParams['pdf.fonttype']=42
__ORIG_SF=plt.savefig
def __sf_all(fname,*a,**k):
    import os as _os
    b,_=_os.path.splitext(str(fname)); k.pop('format',None); k.setdefault('bbox_inches','tight')
    __ORIG_SF(b+'.png',*a,**dict(k))
    __ORIG_SF(b+'.svg',*a,**dict(k,format='svg'))
    __ORIG_SF(b+'.pdf',*a,**dict(k,format='pdf'))
plt.savefig=__sf_all
from matplotlib.lines import Line2D
from adjustText import adjust_text
J=REPO_ROOT + "/jessa2022_metabolic_ssgsea"
res=pd.read_csv(f"{J}/DE_H3K27M_vs_nonK27M_glioma.tsv",sep='\t',index_col=0)
gene_cat=json.load(open(f"{J}/gene_metabolic_category.json"))
res['cat']=[gene_cat.get(g) for g in res.index]
res['nlp']=-np.log10(res['padj'].clip(lower=res.loc[res['padj']>0,'padj'].min()))
res['lfc']=res['log2FoldChange']
cats=['Lipid','Amino acid','Carbohydrate','Nucleotide','Energy','Glycan','Cofactor/Vitamin','Xenobiotic']
colors={'Lipid':'#d62728','Amino acid':'#2ca02c','Carbohydrate':'#1f77b4','Nucleotide':'#9467bd',
        'Energy':'#ff7f0e','Glycan':'#17becf','Cofactor/Vitamin':'#8c564b','Xenobiotic':'#e377c2'}
LFC_T,P_T=1.0,0.05; YCAP=16
fig,ax=plt.subplots(figsize=(13.5,9.5))
bg=res[res['cat'].isna()]
ax.scatter(bg['lfc'],bg['nlp'].clip(upper=YCAP),s=6,c='#dcdcdc',alpha=0.5,edgecolors='none',rasterized=True,zorder=1)
for c in cats:
    d=res[res['cat']==c]
    ax.scatter(d['lfc'],d['nlp'].clip(upper=YCAP),s=30,c=colors[c],alpha=0.9,edgecolors='k',linewidths=0.3,zorder=3)
ax.axhline(-np.log10(P_T),color='0.4',ls='--',lw=0.9); ax.axvline(LFC_T,color='0.4',ls='--',lw=0.9); ax.axvline(-LFC_T,color='0.4',ls='--',lw=0.9)
ax.set_xlim(-6,6); ax.set_ylim(-1,YCAP+0.5)
# labels: top metabolic sig genes both directions + key lipid drivers
sig=res[(res['padj']<P_T)&(res['lfc'].abs()>LFC_T)&(res['cat'].notna())].copy()
lab=sig[sig['nlp']>4]
texts=[ax.text(r['lfc'],min(r['nlp'],YCAP),g,fontsize=13,fontweight='bold',color=colors[r['cat']]) for g,r in lab.iterrows()]
adjust_text(texts,ax=ax,arrowprops=dict(arrowstyle='-',color='0.6',lw=0.5),
            expand_points=(1.5,1.8),expand_text=(1.35,1.5),force_text=(0.8,1.1),force_static=(0.3,0.4),only_move={'text':'xy'},max_move=(40,60))
ax.set_xlabel('log2 fold change  (H3K27M-DMG  vs  non-DMG glioma)',fontsize=15)
ax.set_ylabel('-log10 adjusted p-value',fontsize=15)
ax.annotate('higher in K27M-DMG →',(5.9,-0.6),ha='right',fontsize=12,color='0.35')
ax.annotate('← higher in non-DMG glioma',(-5.9,-0.6),ha='left',fontsize=12,color='0.35')
handles=[]
for c in cats:
    n=int(((res['cat']==c)&(res['padj']<P_T)&(res['lfc'].abs()>LFC_T)).sum())
    up=int(((res['cat']==c)&(res['padj']<P_T)&(res['lfc']>LFC_T)).sum())
    handles.append(Line2D([0],[0],marker='o',color='w',markerfacecolor=colors[c],markeredgecolor='k',markersize=8,label=f'{c}  (n={n}, {up}↑)'))
handles.append(Line2D([0],[0],marker='o',color='w',markerfacecolor='#dcdcdc',markersize=7,label='non-metabolic'))
ax.legend(handles=handles,title='KEGG metabolic category (significant)',loc='upper left',bbox_to_anchor=(1.005,1.0),fontsize=12,title_fontsize=13,frameon=True,framealpha=0.95)
ax.tick_params(labelsize=12)
for sp in ['top','right']: ax.spines[sp].set_visible(False)
ax.set_title("Metabolic genes among H3K27M-DMG vs non-DMG glioma DEGs",fontsize=16)
plt.tight_layout()
plt.savefig(f"{J}/volcano_metabolic_K27M_vs_nonDMG.png",dpi=200,bbox_inches='tight')
plt.savefig(f"{J}/volcano_metabolic_K27M_vs_nonDMG.svg",format='svg',bbox_inches='tight')
print("saved")
