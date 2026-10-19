import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))  # publication: repo-relative root
import json, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

P=REPO_ROOT + '/'
pal = json.load(open(P+'BCAT1_figure_palette.json'))
C = {'shNT': pal['conditions']['KO']['hex'], 'shBCAT1': pal['genotypes']['KO']['hex']}

df = pd.read_excel(P+'manuscript bar plots/Clonogenic.xlsx', sheet_name='Sheet1')
V='Normalized cell index'
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['DejaVu Sans'],
                     'font.size':9,'axes.linewidth':0.9,'svg.fonttype':'none',
                     'pdf.fonttype':42,'ps.fonttype':42})
rng = np.random.default_rng(7)
def pfmt(p): return f"p = {p:.3f}" if p >= 0.001 else f"p = {p:.1e}"

lines = ['BT245','DIPGXIII','HSJ019']
genos = ['K27M','KO']
cons  = ['shNT','shBCAT1']
xs = {('K27M','shNT'):0, ('K27M','shBCAT1'):0.8, ('KO','shNT'):2.1, ('KO','shBCAT1'):2.9}

fig, axes = plt.subplots(1, 3, figsize=(6.6, 3.2), gridspec_kw={'wspace':0.42})
res = []

for ax, cl in zip(axes, lines):
    sub = df[df['Cell line']==cl]
    hi = sub[V].max()
    for (g, c), x in xs.items():
        v = sub[(sub.Genotype==g)&(sub.Condition==c)][V].values
        ax.bar(x, v.mean(), width=0.7, color=C[c], edgecolor='none',
               yerr=v.std(ddof=1)/np.sqrt(len(v)), capsize=3,
               error_kw=dict(elinewidth=0.9, capthick=0.9), zorder=2)
        ax.scatter(x+rng.uniform(-0.15,0.15,len(v)), v, s=13, facecolor='white',
                   edgecolor='black', linewidth=0.6, zorder=3)

    for g in genos:                       # paired t-test, shNT vs shBCAT1 within genotype
        a = sub[(sub.Genotype==g)&(sub.Condition=='shNT')][V].values
        b = sub[(sub.Genotype==g)&(sub.Condition=='shBCAT1')][V].values
        t, p = stats.ttest_rel(a, b)
        res.append((cl, g, a.mean(), b.mean(), len(a), t, p))
        x1, x2 = xs[(g,'shNT')], xs[(g,'shBCAT1')]
        y = hi*1.10
        ax.plot([x1,x1,x2,x2],[y,y*1.04,y*1.04,y], color='black', lw=0.9)
        ax.text((x1+x2)/2, y*1.06, pfmt(p), ha='center', va='bottom', fontsize=7.5)

    ax.set_ylim(0, hi*1.36); ax.set_xlim(-0.6, 3.5)
    ax.set_xticks([0.4, 2.5]); ax.set_xticklabels(['K27M','KO'], fontsize=9)
    ax.set_ylabel('Normalized cell index')
    ax.set_title(cl, fontsize=10, fontweight='bold', pad=2)
    ax.spines[['top','right']].set_visible(False)
    ax.tick_params(width=0.9, length=3)
    ax.tick_params(axis='x', length=0)

h = [plt.Rectangle((0,0),1,1, fc=C[c], ec='none') for c in cons]
fig.legend(h, cons, loc='upper center', bbox_to_anchor=(0.5,1.05), ncol=2,
           frameon=False, fontsize=8.5, handlelength=1.1, columnspacing=1.4)

out=P+'manuscript bar plots/Clonogenic_barplots'
fig.savefig(out+'.png', dpi=400, bbox_inches='tight')
fig.savefig(out+'.svg', bbox_inches='tight')
fig.savefig(out+'.pdf', bbox_inches='tight')

st = pd.DataFrame(res, columns=['cell_line','genotype','mean_shNT','mean_shBCAT1','n','t','p_two_tailed'])
st['test']='paired two-tailed t-test, shNT vs shBCAT1 within genotype (n=4 replicates, paired by replicate order)'
st.to_csv(P+'manuscript bar plots/Clonogenic_stats.csv', index=False)
print(st.to_string(index=False))
