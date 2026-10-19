import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))  # publication: repo-relative root
import json, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

P=REPO_ROOT + '/'
pal = json.load(open(P+'BCAT1_figure_palette.json'))
C = {k: pal['treatments'][k]['hex'] for k in ['DMSO','TDI-11055']}

df = pd.read_excel(P+'manuscript bar plots/Source data.xlsx', sheet_name='TDI WB')
V='POL2-S2phos/POL2'
if 'replicate' not in df.columns: df['replicate']=1
df['pair'] = df['Cell line'] + ' r' + df['replicate'].astype(str)

plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['DejaVu Sans'],
                     'font.size':9,'axes.linewidth':0.9,'svg.fonttype':'none','pdf.fonttype':42,'ps.fonttype':42})

pairs = df['pair'].unique().tolist()
veh = df[df.Treatment=='DMSO'].set_index('pair')[V].reindex(pairs).values
trt = df[df.Treatment=='TDI-11055'].set_index('pair')[V].reindex(pairs).values
t, p = stats.ttest_rel(veh, trt)                     # paired, n = 4 blots

# secondary: collapse replicates to per-cell-line means (n = 2 independent lines)
pl = df.groupby(['Cell line','Treatment'])[V].mean().unstack()
t2, p2 = stats.ttest_rel(pl['DMSO'].values, pl['TDI-11055'].values)

fig, ax = plt.subplots(figsize=(2.3, 3.1))
offs = np.linspace(-0.19, 0.19, len(pairs))
for x, v, grp in [(0, veh, 'DMSO'), (0.8, trt, 'TDI-11055')]:
    ax.bar(x, v.mean(), width=0.7, color=C[grp], edgecolor='none',
           yerr=v.std(ddof=1)/np.sqrt(len(v)), capsize=3,
           error_kw=dict(elinewidth=0.9, capthick=0.9), zorder=2)
    ax.scatter(x+offs, v, s=14, facecolor='white', edgecolor='black',
               linewidth=0.7, zorder=4)
for i in range(len(pairs)):
    ax.plot([0+offs[i], 0.8+offs[i]], [veh[i], trt[i]], color='0.45', lw=0.7, zorder=3)

hi = max(veh.max(), trt.max()); y = hi*1.12
ax.plot([0,0,0.8,0.8],[y,y*1.04,y*1.04,y], color='black', lw=0.9)
ax.text(0.4, y*1.07, f"p = {p:.3f}", ha='center', va='bottom', fontsize=8.5)
ax.set_ylim(0, hi*1.34); ax.set_xlim(-0.6, 1.4)
ax.set_xticks([0,0.8]); ax.set_xticklabels(['DMSO','TDI-11055'], fontsize=8)
ax.set_ylabel('POL2 S2phos / total POL2\n(normalized to DMSO)')
ax.spines[['top','right']].set_visible(False)
ax.tick_params(width=0.9, length=3)

out=P+'manuscript bar plots/TDI_WB_barplot'
fig.savefig(out+'.png', dpi=400, bbox_inches='tight')
fig.savefig(out+'.svg', bbox_inches='tight')
fig.savefig(out+'.pdf', bbox_inches='tight')

pd.DataFrame([
 {'comparison':'TDI-11055 vs DMSO','unit':'blot (2 lines x 2 replicates)','n_pairs':len(pairs),
  'mean_DMSO':veh.mean(),'mean_TDI':trt.mean(),'t':t,'p_two_tailed':p,
  'test':'paired two-tailed t-test, n=4 blots (plotted)'},
 {'comparison':'TDI-11055 vs DMSO','unit':'cell line (replicates averaged)','n_pairs':len(pl),
  'mean_DMSO':pl['DMSO'].mean(),'mean_TDI':pl['TDI-11055'].mean(),'t':t2,'p_two_tailed':p2,
  'test':'paired two-tailed t-test on per-line means, n=2 lines (conservative)'},
]).to_csv(P+'manuscript bar plots/TDI_WB_stats.csv', index=False)
print(dict(zip(pairs,trt)))
print('n=4 blots  t=%.3f p=%.4f'%(t,p))
print('n=2 lines  t=%.3f p=%.4f'%(t2,p2))
