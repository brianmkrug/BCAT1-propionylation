import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))  # publication: repo-relative root
import json, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

P=REPO_ROOT + '/'
pal = json.load(open(P+'BCAT1_figure_palette.json'))
C = {g: pal['cell_line_groups'][g]['hex'] for g in ['DMG','non-DMG']}

df = pd.read_excel(P+'manuscript bar plots/Source data.xlsx', sheet_name='Cell line panel WBs')
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['DejaVu Sans'],
                     'font.size':9,'axes.linewidth':0.9,'svg.fonttype':'none','pdf.fonttype':42,'ps.fonttype':42})
rng = np.random.default_rng(4)
def pfmt(p): return f"p = {p:.3f}" if p >= 0.001 else f"p = {p:.1e}"

targets = ['BCAT1','BCAT2','H3K9pr']
fig, axes = plt.subplots(1, 3, figsize=(4.1, 1.85), gridspec_kw={'wspace':0.85})
res = []

for ax, tgt in zip(axes, targets):
    a = df[df.Group=='DMG'][tgt].values
    b = df[df.Group=='non-DMG'][tgt].values
    t, p = stats.ttest_ind(a, b, equal_var=False)
    res.append((tgt, a.mean(), b.mean(), t, p))

    for x, v, grp in [(0, a, 'DMG'), (0.62, b, 'non-DMG')]:
        ax.bar(x, v.mean(), width=0.46, color=C[grp], edgecolor='none',
               yerr=v.std(ddof=1)/np.sqrt(len(v)), capsize=3,
               error_kw=dict(elinewidth=0.9, capthick=0.9), zorder=2)
        ax.scatter(x+rng.uniform(-0.10,0.10,len(v)), v, s=14, facecolor='white',
                   edgecolor='black', linewidth=0.7, zorder=3)

    hi = max(df[tgt].max(),
             max(v.mean()+v.std(ddof=1)/np.sqrt(len(v)) for v in (a,b)))
    y = hi*1.12
    ax.plot([0,0,0.62,0.62],[y,y*1.04,y*1.04,y], color='black', lw=0.9)
    ax.text(0.31, y*1.07, pfmt(p), ha='center', va='bottom', fontsize=8.5)
    ax.set_ylim(0, hi*1.34); ax.set_xlim(-0.48, 1.10)
    ax.set_xticks([0,0.62])
    ax.set_xticklabels(['DMG','non-DMG'], fontsize=8.5,
                       rotation=45, ha='right', rotation_mode='anchor')
    ax.set_ylabel('%s (normalized)' % tgt)
    ax.set_title(tgt, fontsize=10, fontweight='bold', pad=2)
    ax.spines[['top','right']].set_visible(False)
    ax.tick_params(width=0.9, length=3)

out=P+'manuscript bar plots/Cell_line_panel_WB_barplots'
fig.savefig(out+'.png', dpi=400, bbox_inches='tight')
fig.savefig(out+'.svg', bbox_inches='tight')
fig.savefig(out+'.pdf', bbox_inches='tight')

st = pd.DataFrame(res, columns=['target','mean_DMG','mean_nonDMG','t','p_two_tailed'])
st['test']="Welch's unpaired two-tailed t-test (n=4 cell lines/group)"
st.to_csv(P+'manuscript bar plots/Cell_line_panel_WB_stats.csv', index=False)
print(st.to_string(index=False))
