import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))  # publication: repo-relative root
import json, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

P=REPO_ROOT + '/'
pal = json.load(open(P+'BCAT1_figure_palette.json'))
C = {'shNT': pal['conditions']['KO']['hex'], 'shBCAT1': pal['genotypes']['KO']['hex']}  # #2166ac / #f0a020

df = pd.read_excel(P+'manuscript bar plots/Source data.xlsx', sheet_name='S6phos DMG WBs').dropna(subset=['construct'])
df['cell line'] = df['cell line'].str.replace('DIPXIII','DIPGXIII')
V='pS6/S6 signal'

plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['DejaVu Sans'],
                     'font.size':9,'axes.linewidth':0.9,'svg.fonttype':'none','pdf.fonttype':42,'ps.fonttype':42})

nt = df[df.construct=='shNT'].set_index('cell line')[V]
kd = df[df.construct=='shBCAT1'].set_index('cell line')[V]
lines = list(kd.index)
ratios = kd.values                      # each already normalized to its own shNT = 1
t, p = stats.ttest_1samp(ratios, 1.0)   # ratio-paired one-sample t-test vs 1

fig, ax = plt.subplots(figsize=(2.3, 3.1))
for x, grp, vals in [(0,'shNT',nt.reindex(lines).values), (0.8,'shBCAT1',ratios)]:
    ax.bar(x, vals.mean(), width=0.7, color=C[grp], edgecolor='none',
           yerr=vals.std(ddof=1)/np.sqrt(len(vals)), capsize=3,
           error_kw=dict(elinewidth=0.9, capthick=0.9), zorder=2)
    for cl, v, off in zip(lines, vals, [-0.13, 0.13]):
        ax.scatter(x+off, v, s=14, facecolor='white',
                   edgecolor='black', linewidth=0.7, zorder=3)

hi = max(ratios.max(), 1.0, ratios.mean()+ratios.std(ddof=1)/np.sqrt(len(ratios)))
y = hi*1.12
ax.plot([0,0,0.8,0.8],[y,y*1.04,y*1.04,y], color='black', lw=0.9)
ax.text(0.4, y*1.07, f"p = {p:.3f}", ha='center', va='bottom', fontsize=8.5)
ax.set_ylim(0, hi*1.34); ax.set_xlim(-0.6, 1.4)
ax.set_xticks([0,0.8]); ax.set_xticklabels(['shNT','shBCAT1'])
ax.set_ylabel('pS6 / total S6\n(normalized to shNT)')
ax.set_title('pS6 (DMG lines)', fontsize=10, fontweight='bold', pad=2)
ax.spines[['top','right']].set_visible(False)
ax.tick_params(width=0.9, length=3)

out=P+'manuscript bar plots/S6phos_DMG_WB_barplot'
fig.savefig(out+'.png', dpi=400, bbox_inches='tight')
fig.savefig(out+'.svg', bbox_inches='tight')
fig.savefig(out+'.pdf', bbox_inches='tight')

pd.DataFrame([{'comparison':'shBCAT1 vs shNT (pooled BT245 + DIPGXIII)',
               'n_cell_lines':len(ratios),'mean_shNT':1.0,'mean_shBCAT1':ratios.mean(),
               'BT245':kd.get('BT245'),'DIPGXIII':kd.get('DIPGXIII'),'t':t,'p_two_tailed':p,
               'test':'ratio-paired one-sample t-test of shBCAT1/shNT vs 1 (n=2 lines, 1 blot each)'}]
            ).to_csv(P+'manuscript bar plots/S6phos_DMG_WB_stats.csv', index=False)
print(dict(zip(lines,ratios)), 'mean=%.3f t=%.3f p=%.3f'%(ratios.mean(),t,p))
