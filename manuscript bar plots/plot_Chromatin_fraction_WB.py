import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))  # publication: repo-relative root
import json, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

P=REPO_ROOT + '/'
pal = json.load(open(P+'BCAT1_figure_palette.json'))
C = {'shNT': pal['conditions']['KO']['hex'], 'shBCAT1': pal['genotypes']['KO']['hex']}

df = pd.read_excel(P+'manuscript bar plots/Source data.xlsx', sheet_name='Chromatin fraction WBs')
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['DejaVu Sans'],
                     'font.size':9,'axes.linewidth':0.9,'svg.fonttype':'none','pdf.fonttype':42,'ps.fonttype':42})
rng = np.random.default_rng(3)
def pfmt(p): return f"p = {p:.3f}" if p >= 0.001 else f"p = {p:.1e}"

targets = ['H3K9pr','ENL','BRD4']
labels  = {'H3K9pr':'H3K9pr','ENL':'ENL','BRD4':'BRD4'}
fig, axes = plt.subplots(1, 3, figsize=(5.6, 3.1), gridspec_kw={'wspace':0.55})
res = []

for ax, tgt in zip(axes, targets):
    nt = df[df.construct=='shNT'][tgt].values
    kd = df[df.construct=='shBCAT1'][tgt].values
    t, p = stats.ttest_1samp(kd, 1.0)     # ratio-paired one-sample t-test vs shNT = 1
    res.append((tgt, nt.mean(), kd.mean(), kd.std(ddof=1)/np.sqrt(len(kd)), t, p))

    for x, v, grp in [(0, nt, 'shNT'), (0.8, kd, 'shBCAT1')]:
        ax.bar(x, v.mean(), width=0.7, color=C[grp], edgecolor='none',
               yerr=v.std(ddof=1)/np.sqrt(len(v)), capsize=3,
               error_kw=dict(elinewidth=0.9, capthick=0.9), zorder=2)
        ax.scatter(x+rng.uniform(-0.13,0.13,len(v)), v, s=14, facecolor='white',
                   edgecolor='black', linewidth=0.7, zorder=3)

    hi = max(nt.max(), kd.max(), 1.0)
    y = hi*1.12
    ax.plot([0,0,0.8,0.8],[y,y*1.04,y*1.04,y], color='black', lw=0.9)
    ax.text(0.4, y*1.07, pfmt(p), ha='center', va='bottom', fontsize=8.5)
    ax.set_ylim(0, hi*1.34); ax.set_xlim(-0.6, 1.4)
    ax.set_xticks([0,0.8]); ax.set_xticklabels(['shNT','shBCAT1'], fontsize=8)
    ax.set_ylabel('%s (chromatin)\n(normalized to shNT)' % labels[tgt])
    ax.set_title(labels[tgt], fontsize=10, fontweight='bold', pad=2)
    ax.spines[['top','right']].set_visible(False)
    ax.tick_params(width=0.9, length=3)

out=P+'manuscript bar plots/Chromatin_fraction_WB_barplots'
fig.savefig(out+'.png', dpi=400, bbox_inches='tight')
fig.savefig(out+'.svg', bbox_inches='tight')
fig.savefig(out+'.pdf', bbox_inches='tight')

st = pd.DataFrame(res, columns=['target','mean_shNT','mean_shBCAT1','sem_shBCAT1','t','p_two_tailed'])
st['test']="ratio-paired one-sample t-test of shBCAT1/shNT vs 1 (n=3 replicates)"
st.to_csv(P+'manuscript bar plots/Chromatin_fraction_WB_stats.csv', index=False)
print(st.to_string(index=False))
