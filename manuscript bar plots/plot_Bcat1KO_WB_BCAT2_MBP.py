import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))  # publication: repo-relative root
import json, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

pal = json.load(open(REPO_ROOT + '/BCAT1_figure_palette.json'))
C = {'WT': pal['genotypes']['WT']['hex'], 'KO': pal['genotypes']['KO']['hex']}

df = pd.read_excel(REPO_ROOT + '/manuscript bar plots/Source data.xlsx', sheet_name='Bcat1-KO WB')

plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['DejaVu Sans'],
                     'font.size':9,'axes.linewidth':0.9,'svg.fonttype':'none','pdf.fonttype':42,'ps.fonttype':42})

def pfmt(p):
    return f"p = {p:.3f}" if p >= 0.001 else f"p = {p:.1e}"

fig, axes = plt.subplots(1, 2, figsize=(6.0, 3.1),
                         gridspec_kw={'width_ratios':[2.0,1.0],'wspace':0.42})
rng = np.random.default_rng(0)
res = []

# ---- Panel A: BCAT2, genotype x age ----
ax = axes[0]
groups = [('P6','WT'),('P6','KO'),('P90','WT'),('P90','KO')]
xs = [0, 0.8, 2.1, 2.9]
for x,(age,gt) in zip(xs, groups):
    v = df[(df.age==age)&(df.genotype==gt)]['BCAT2'].values
    ax.bar(x, v.mean(), width=0.7, color=C[gt], edgecolor='none',
           yerr=v.std(ddof=1)/np.sqrt(len(v)), capsize=3,
           error_kw=dict(elinewidth=0.9, capthick=0.9), zorder=2)
    ax.scatter(x+rng.uniform(-0.13,0.13,len(v)), v, s=14, facecolor='white',
               edgecolor='black', linewidth=0.7, zorder=3)

ymax = df['BCAT2'].max()
top = ymax*1.36
for (x1,x2),age in zip([(xs[0],xs[1]),(xs[2],xs[3])], ['P6','P90']):
    a = df[(df.age==age)&(df.genotype=='WT')]['BCAT2'].values
    b = df[(df.age==age)&(df.genotype=='KO')]['BCAT2'].values
    t,p = stats.ttest_ind(a,b,equal_var=False)
    res.append(('BCAT2',age,a.mean(),b.mean(),t,p))
    y = ymax*1.16
    ax.plot([x1,x1,x2,x2],[y,y*1.04,y*1.04,y],color='black',lw=0.9,clip_on=False)
    ax.text((x1+x2)/2, y*1.07, pfmt(p), ha='center', va='bottom', fontsize=8.5)

ax.set_ylim(0, top); ax.set_xlim(-0.6, 3.5)
ax.set_xticks(xs); ax.set_xticklabels([g for _,g in groups])
ax.set_ylabel('BCAT2 / loading control\n(normalized)')
ax.set_title('BCAT2', fontsize=10, fontweight='bold', pad=2)
for x,lab in zip([(xs[0]+xs[1])/2,(xs[2]+xs[3])/2], ['P6','P90']):
    ax.text(x, -0.155, lab, transform=ax.get_xaxis_transform(), ha='center',
            va='top', fontsize=9.5, fontweight='bold')
    ax.plot([x-0.75,x+0.75],[-0.115]*2, transform=ax.get_xaxis_transform(),
            color='black', lw=0.9, clip_on=False)

# ---- Panel B: MBP, P6 only ----
ax = axes[1]
xs2 = [0, 0.8]
vals = {}
for x,gt in zip(xs2, ['WT','KO']):
    v = df[(df.age=='P6')&(df.genotype==gt)]['MBP'].values; vals[gt]=v
    ax.bar(x, v.mean(), width=0.7, color=C[gt], edgecolor='none',
           yerr=v.std(ddof=1)/np.sqrt(len(v)), capsize=3,
           error_kw=dict(elinewidth=0.9, capthick=0.9), zorder=2)
    ax.scatter(x+rng.uniform(-0.13,0.13,len(v)), v, s=14, facecolor='white',
               edgecolor='black', linewidth=0.7, zorder=3)
t,p = stats.ttest_ind(vals['WT'],vals['KO'],equal_var=False)
res.append(('MBP','P6',vals['WT'].mean(),vals['KO'].mean(),t,p))
m = df[df.age=='P6']['MBP'].max(); y = m*1.16
ax.plot([xs2[0],xs2[0],xs2[1],xs2[1]],[y,y*1.04,y*1.04,y],color='black',lw=0.9)
ax.text(0.4, y*1.07, pfmt(p), ha='center', va='bottom', fontsize=8.5)
ax.set_ylim(0, m*1.36); ax.set_xlim(-0.6,1.4)
ax.set_xticks(xs2); ax.set_xticklabels(['WT','KO'])
ax.set_ylabel('MBP / loading control\n(normalized)')
ax.set_title('MBP (P6)', fontsize=10, fontweight='bold', pad=2)
ax.text(0.4, -0.155, 'P6', transform=ax.get_xaxis_transform(), ha='center',
        va='top', fontsize=9.5, fontweight='bold')
ax.plot([-0.35,1.15],[-0.115]*2, transform=ax.get_xaxis_transform(),
        color='black', lw=0.9, clip_on=False)

for ax in axes:
    ax.spines[['top','right']].set_visible(False)
    ax.tick_params(width=0.9, length=3)

h = [plt.Rectangle((0,0),1,1,fc=C[g],ec='none') for g in ['WT','KO']]
fig.legend(h, ['WT','Bcat1-KO'], loc='upper right', bbox_to_anchor=(0.995,1.02),
           frameon=False, fontsize=8.5, handlelength=1.1, ncol=2, columnspacing=0.9)

out=REPO_ROOT + '/manuscript bar plots/Bcat1KO_WB_BCAT2_MBP'
fig.savefig(out+'.png', dpi=400, bbox_inches='tight')
fig.savefig(out+'.svg', bbox_inches='tight')
fig.savefig(out+'.pdf', bbox_inches='tight')

st = pd.DataFrame(res, columns=['protein','age','mean_WT','mean_KO','t','p_two_tailed'])
st['test']="Welch's unpaired two-tailed t-test (n=3/group)"
st.to_csv(REPO_ROOT + '/manuscript bar plots/Bcat1KO_WB_stats.csv', index=False)
print(st.to_string(index=False))
