"""KO cell-type composition along pooled5 pseudotime, split by library, to show
what r2 contributes and whether r1/r3 patterns are preserved."""
import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))  # publication: repo-relative root
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

W = REPO_ROOT + "/scXenograft/work"
FIG = REPO_ROOT + "/scXenograft/figures"
CLS = ["RGC","GLIP","OPC","OL","ASTRO","EPEN","NEURON"]
COL = dict(zip(CLS, ["#7f7f7f","#9467bd","#e6a010","#f4d03f","#2ca089","#e8746a","#7fbf3f"]))
NBIN = 24
emb = pd.read_csv(f"{W}/pooled5_dpt.csv.gz", index_col=0)
emb["label"] = emb[CLS].idxmax(1)
lo, hi = emb.dpt.min(), emb.dpt.max()
edges = np.linspace(lo, hi, NBIN+1); ctr = (edges[:-1]+edges[1:])/2
emb["fb"] = np.clip(np.digitize(emb.dpt, edges)-1, 0, NBIN-1)
def sm(y):
    k=np.array([.25,.5,.25]); return np.convolve(np.r_[y[0],y,y[-1]],k,"same")[1:-1]

ko = emb[emb.genotype=="K27M-KO"]
libs = [("GSM6209153_hsj019_ko_r1","r1 (702)"),("GSM6209154_hsj019_ko_r2","r2 (3824)"),
        ("GSM6209155_hsj019_ko_r3","r3 (3702)")]
fig, axes = plt.subplots(2, 3, figsize=(16, 8), sharex=True)
for j,(lib,ttl) in enumerate(libs):
    sub = ko[ko["sample"]==lib]
    frac = np.zeros((len(CLS),NBIN)); n=np.zeros(NBIN)
    for b in range(NBIN):
        d=sub[sub.fb==b]; n[b]=len(d)
        for i,c in enumerate(CLS): frac[i,b]=(d.label==c).mean() if len(d) else 0
    ax=axes[0,j]; bottom=np.zeros(NBIN)
    for i,c in enumerate(CLS):
        y=sm(frac[i]); ax.fill_between(ctr,bottom,bottom+y,color=COL[c],lw=0,label=c); bottom+=y
    ax.set_ylim(0,1); ax.set_title(f"KO {ttl}\nhard-label composition",fontsize=10)
    axc=ax.twinx(); axc.plot(ctr,n,"k:",lw=1); axc.set_ylim(0,max(n.max(),1)*1.1); axc.tick_params(labelsize=6)
    ax2=axes[1,j]
    for i,c in enumerate(CLS): ax2.plot(ctr,sm(frac[i]),color=COL[c],lw=2)
    ax2.set_ylim(0,None); ax2.set_xlabel("pooled pseudotime"); ax2.set_title("per-type lines",fontsize=9)
    for a in (ax,ax2):
        for s in ["top","right"]: a.spines[s].set_visible(False)
axes[0,0].legend(loc="upper right",fontsize=6,frameon=False)
axes[0,0].set_ylabel("fraction"); axes[1,0].set_ylabel("fraction")
fig.suptitle("KO cell-type composition along pooled pseudotime, split by library "
             "(r2 included) — is the pattern consistent across replicates?",fontsize=12)
fig.tight_layout(); fig.savefig(f"{FIG}/celltype5_KO_by_library.png",dpi=160,bbox_inches="tight")
plt.close(fig)

# quantify r2 vs r1/r3 endpoint composition (late pseudotime, dpt>median)
late = ko[ko.dpt > ko.dpt.median()]
print("KO late-trajectory (dpt>median) label fraction by library:")
print((pd.crosstab(late["sample"].str.split("_").str[-1], late.label, normalize="index")*100).round(1).to_string())
print("\nfigure written")
