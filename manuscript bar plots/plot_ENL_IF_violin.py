import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))  # publication: repo-relative root
import json, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

P=REPO_ROOT + '/'
pal = json.load(open(P+'BCAT1_figure_palette.json'))
BLUE, YELLOW = pal['conditions']['KO']['hex'], pal['genotypes']['KO']['hex']

df = pd.read_excel(P+'manuscript bar plots/ENL_IF/ENL_quant.xlsx', sheet_name='grouped')
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['DejaVu Sans'],
                     'font.size':9,'axes.linewidth':0.9,'svg.fonttype':'none','pdf.fonttype':42,'ps.fonttype':42})
rng = np.random.default_rng(5)
def pfmt(p): return f"p = {p:.3f}" if p >= 0.001 else f"p = {p:.1e}"

panels = [('shNT','shBCAT1'), ('SCM(+)','SCM(-)')]
fig, axes = plt.subplots(1, 2, figsize=(4.4, 3.1), gridspec_kw={'wspace':0.55})

for ax, (g1, g2) in zip(axes, panels):
    a = df[g1].dropna().values
    b = df[g2].dropna().values
    t, p = stats.ttest_ind(a, b, equal_var=False)

    parts = ax.violinplot([a, b], positions=[0, 0.8], widths=0.66,
                          showmeans=False, showmedians=False, showextrema=False)
    for body, col in zip(parts['bodies'], [BLUE, YELLOW]):
        body.set_facecolor(col); body.set_edgecolor('none')
        body.set_linewidth(0); body.set_alpha(1.0)

    for x, v in [(0, a), (0.8, b)]:
        q1, med, q3 = np.percentile(v, [25, 50, 75])
        ax.vlines(x, q1, q3, color='black', lw=3.2, zorder=4)
        ax.scatter(x+rng.uniform(-0.16,0.16,len(v)), v, s=7, facecolor='white',
                   edgecolor='0.25', linewidth=0.45, alpha=0.9, zorder=3)
        ax.scatter(x, med, s=16, facecolor='white', edgecolor='black',
                   linewidth=0.8, zorder=5)

    lo = min(a.min(), b.min()); hi = max(a.max(), b.max()); span = hi-lo
    y = hi + span*0.10
    ax.plot([0,0,0.8,0.8],[y, y+span*0.035, y+span*0.035, y], color='black', lw=0.9)
    ax.text(0.4, y+span*0.055, pfmt(p), ha='center', va='bottom', fontsize=8.5)
    ax.set_ylim(max(0, lo-span*0.12), hi+span*0.30); ax.set_xlim(-0.6, 1.4)
    ax.set_xticks([0,0.8]); ax.set_xticklabels([g1, g2], fontsize=8.5)
    ax.set_ylabel('ENL IF signal\n(normalized)')
    ax.spines[['top','right']].set_visible(False)
    ax.tick_params(width=0.9, length=3)

out=P+'manuscript bar plots/ENL_IF_violinplots'
fig.savefig(out+'.png', dpi=400, bbox_inches='tight')
fig.savefig(out+'.svg', bbox_inches='tight')
fig.savefig(out+'.pdf', bbox_inches='tight')
print('done')
