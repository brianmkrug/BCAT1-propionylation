import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))  # publication: repo-relative root
import json, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

P=REPO_ROOT + '/'
pal = json.load(open(P+'BCAT1_figure_palette.json'))
C = {'WT': pal['genotypes']['WT']['hex'], 'KO': pal['genotypes']['KO']['hex']}
df = pd.read_excel(P+'manuscript bar plots/Source data.xlsx', sheet_name='S6phos brain WBs')
V = 'pS6/S6 signal'

plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['DejaVu Sans'],
                     'font.size':9,'axes.linewidth':0.9,'svg.fonttype':'none','pdf.fonttype':42,'ps.fonttype':42})
rng = np.random.default_rng(2)
def pfmt(p): return f"p = {p:.3f}" if p >= 0.001 else f"p = {p:.1e}"

fig, ax = plt.subplots(figsize=(3.2, 3.1))
groups = [('P6','WT'),('P6','KO'),('P90','WT'),('P90','KO')]
xs = [0, 0.8, 2.1, 2.9]
res = []
for x,(age,gt) in zip(xs, groups):
    v = df[(df.age==age)&(df.genotype==gt)][V].values
    ax.bar(x, v.mean(), width=0.7, color=C[gt], edgecolor='none',
           yerr=v.std(ddof=1)/np.sqrt(len(v)), capsize=3,
           error_kw=dict(elinewidth=0.9, capthick=0.9), zorder=2)
    ax.scatter(x+rng.uniform(-0.13,0.13,len(v)), v, s=14, facecolor='white',
               edgecolor='black', linewidth=0.7, zorder=3)

hi = max(df[V].max(),
         max((df[(df.age==a)&(df.genotype==g)][V].mean()
              + df[(df.age==a)&(df.genotype==g)][V].std(ddof=1)/np.sqrt(3)) for a,g in groups))
for (x1,x2),age in zip([(xs[0],xs[1]),(xs[2],xs[3])], ['P6','P90']):
    a = df[(df.age==age)&(df.genotype=='WT')][V].values
    b = df[(df.age==age)&(df.genotype=='KO')][V].values
    t,p = stats.ttest_ind(a,b,equal_var=False)
    res.append(('pS6/S6',age,a.mean(),b.mean(),t,p))
    y = hi*1.12
    ax.plot([x1,x1,x2,x2],[y,y*1.04,y*1.04,y], color='black', lw=0.9)
    ax.text((x1+x2)/2, y*1.07, pfmt(p), ha='center', va='bottom', fontsize=8.5)

ax.set_ylim(0, hi*1.34); ax.set_xlim(-0.6, 3.5)
ax.set_xticks(xs); ax.set_xticklabels([g for _,g in groups])
ax.set_ylabel('pS6 / total S6\n(normalized)')
ax.set_title('pS6 (brain)', fontsize=10, fontweight='bold', pad=2)
for x,lab in zip([(xs[0]+xs[1])/2,(xs[2]+xs[3])/2], ['P6','P90']):
    ax.text(x, -0.155, lab, transform=ax.get_xaxis_transform(), ha='center',
            va='top', fontsize=9.5, fontweight='bold')
    ax.plot([x-0.75,x+0.75],[-0.115]*2, transform=ax.get_xaxis_transform(),
            color='black', lw=0.9, clip_on=False)
ax.spines[['top','right']].set_visible(False)
ax.tick_params(width=0.9, length=3)

out=P+'manuscript bar plots/S6phos_brain_WB_barplot'
fig.savefig(out+'.png', dpi=400, bbox_inches='tight')
fig.savefig(out+'.svg', bbox_inches='tight')
fig.savefig(out+'.pdf', bbox_inches='tight')

st = pd.DataFrame(res, columns=['readout','age','mean_WT','mean_KO','t','p_two_tailed'])
st['test']="Welch's unpaired two-tailed t-test (n=3/group)"
st.to_csv(P+'manuscript bar plots/S6phos_brain_WB_stats.csv', index=False)
print(st.to_string(index=False))
# secondary: age effect within genotype
for gt in ['WT','KO']:
    a=df[(df.genotype==gt)&(df.age=='P6')][V]; b=df[(df.genotype==gt)&(df.age=='P90')][V]
    print(gt,'P6 vs P90 p=%.3f'%stats.ttest_ind(a,b,equal_var=False)[1])
