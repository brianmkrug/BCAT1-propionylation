import os, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42

from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator
from scipy import stats

HERE=os.path.dirname(os.path.abspath(__file__)); BASE=os.path.dirname(HERE)
XLSX=os.path.join(BASE,"norm_cell_data copy.xlsx"); FIG=os.path.join(HERE,"figures")
GENO_C={"K27M":"#d1495b","KO":"#2e86ab"}; LINES=["BT245","DIPGXIII","HSJ019"]
# reduced first, then oxidized / disulfide / conjugates+adducts
PANELS=[("GSH\n(reduced)","glutathione, reduced (GSH)"),
        ("GSSG\n(oxidized)","glutathione, oxidized (GSSG)"),
        ("Cys-GSH\ndisulfide","cysteine-glutathione disulfide"),
        ("Cyclic dGSH","cyclic dGSH"),
        ("GSNO","S-nitrosoglutathione (GSNO)"),
        ("S-lactoyl-GSH","S-lactoylglutathione"),
        ("S-methyl-GSH","S-methylglutathione"),
        ("HNE-GSH","4-hydroxy-nonenal-glutathione"),
        ("CoA-GSH","CoA-glutathione*"),
        ("S-DCE-GSH","S-(1,2-dicarboxyethyl)glutathione")]
FS=6
plt.rcParams.update({"font.size":FS,"axes.labelsize":FS,"xtick.labelsize":FS,"ytick.labelsize":FS,"legend.fontsize":FS,
    "axes.facecolor":"white","figure.facecolor":"white","axes.edgecolor":"black","axes.linewidth":0.35,"axes.grid":False})

cc=pd.read_excel(XLSX,sheet_name="cell culture")
cc=cc[cc["Cell line"].isin(LINES)].reset_index(drop=True)
def pfmt(p): return f"p={p:.0e}" if p<1e-2 else f"p={p:.2f}"

np.random.seed(7)
NROW,NCOL=2,5
fig,axes=plt.subplots(NROW,NCOL,figsize=(18/2.54,6.5/2.54)); axes=axes.ravel()
for idx,(ax,(title,met)) in enumerate(zip(axes,PANELS)):
    vals=pd.to_numeric(cc[met],errors="coerce"); xpos,width=np.arange(len(LINES)),0.36
    for gi,geno in enumerate(["K27M","KO"]):
        for li,line in enumerate(LINES):
            m=(cc["Cell line"]==line)&(cc["Genotype"]==geno); v=vals[m.values].dropna().values
            x=xpos[li]+(gi-0.5)*width
            if len(v):
                ax.bar(x,v.mean(),width=width,color=GENO_C[geno],alpha=.85,edgecolor="none",
                       yerr=(v.std(ddof=1)/np.sqrt(len(v)) if len(v)>1 else 0),capsize=1.2,error_kw=dict(lw=.5))
                ax.scatter(np.full(len(v),x)+np.random.uniform(-.05,.05,len(v)),v,s=3,c="white",edgecolor="k",lw=.35,zorder=5)
    ymax=np.nanmax(vals.values); nz=vals[vals>0]; hm=(nz.min()/2.0) if len(nz) else np.nan
    imp=np.log10(vals.replace(0,np.nan).fillna(hm))
    for li,line in enumerate(LINES):
        k=imp[((cc["Cell line"]==line)&(cc["Genotype"]=="K27M")).values].dropna()
        o=imp[((cc["Cell line"]==line)&(cc["Genotype"]=="KO")).values].dropna()
        if len(k)>=2 and len(o)>=2:
            _,p=stats.ttest_ind(k,o,equal_var=False)
            if np.isfinite(p):
                yb=ymax*1.03
                ax.plot([xpos[li]-width/2,xpos[li]+width/2],[yb,yb],color="k",lw=.4)
                ax.text(xpos[li],yb*1.004,pfmt(p),ha="center",va="bottom",fontsize=5)
    ax.set_xticks(xpos)
    if idx>=NCOL: ax.set_xticklabels(LINES,rotation=30,ha="right",fontsize=FS)
    else: ax.set_xticklabels([])
    ax.set_title(title,fontsize=FS,fontweight="normal",pad=1.5)
    ax.set_ylim(0,ymax*1.16); ax.tick_params(length=1.2,width=0.35)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=3))
    ax.ticklabel_format(axis="y",style="sci",scilimits=(0,0)); ax.yaxis.get_offset_text().set_fontsize(FS-2)
    if idx%NCOL==0: ax.set_ylabel("norm. abund.",fontsize=FS)
fig.canvas.draw()
for _ax in np.atleast_1d(axes).ravel()[:len(PANELS)]:
    _ot=_ax.yaxis.get_offset_text(); _s=_ot.get_text(); _ot.set_visible(False)
    if _s: _ax.text(-0.02,1.0,_s,transform=_ax.transAxes,ha="right",va="bottom",fontsize=FS-2.5)
fig.legend(handles=[Patch(facecolor=GENO_C["K27M"],edgecolor="k",lw=.4,label="K27M"),
                    Patch(facecolor=GENO_C["KO"],edgecolor="k",lw=.4,label="KO")],
           loc="lower center",ncol=2,frameon=False,fontsize=FS,handlelength=1,handleheight=1,
           columnspacing=1.2,bbox_to_anchor=(0.5,0.0))
fig.subplots_adjust(left=0.07,right=0.99,top=0.90,bottom=0.20,hspace=0.65,wspace=0.55)
out=os.path.join(FIG,"cellculture_glutathione_bars.png")
for ext in ("png","pdf","svg"):
    fig.savefig(out.replace(".png","."+ext),bbox_inches="tight",dpi=600)
print("wrote",out,"(+pdf/svg)")
