import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))  # publication: repo-relative root
import json, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

pal = json.load(open(REPO_ROOT + '/BCAT1_figure_palette.json'))
RED    = pal['conditions']['K27M']['hex']   # #c0392b
BLUE   = pal['conditions']['KO']['hex']     # #2166ac
YELLOW = pal['genotypes']['KO']['hex']      # #f0a020

df = pd.read_excel(REPO_ROOT + '/manuscript bar plots/Source data.xlsx',
                   sheet_name='aKG assay')

plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['DejaVu Sans'],
                     'font.size':9,'axes.linewidth':0.9,'svg.fonttype':'none','pdf.fonttype':42,'ps.fonttype':42})
rng = np.random.default_rng(1)

def pfmt(p): return f"p = {p:.3f}" if p >= 0.001 else f"p = {p:.1e}"

panels = [(('K27M','KO'), (RED, BLUE)),
          (('shNT','shBCAT1'), (BLUE, YELLOW))]

fig, axes = plt.subplots(1, 2, figsize=(4.4, 3.1), gridspec_kw={'wspace':0.55})
res = []

for ax, ((ctrl, test), (c_ctrl, c_test)) in zip(axes, panels):
    c = df[df.Group==ctrl]['Normalized abundance'].values
    t_ = df[df.Group==test]['Normalized abundance'].values
    norm = c.mean()
    c, t_ = c/norm, t_/norm                    # control mean = 1
    tt, p = stats.ttest_ind(c, t_, equal_var=False)
    res.append((ctrl, test, c.mean(), t_.mean(), tt, p))

    for x, v, col in [(0, c, c_ctrl), (0.8, t_, c_test)]:
        ax.bar(x, v.mean(), width=0.7, color=col, edgecolor='none',
               yerr=v.std(ddof=1)/np.sqrt(len(v)), capsize=3,
               error_kw=dict(elinewidth=0.9, capthick=0.9), zorder=2)
        ax.scatter(x+rng.uniform(-0.13,0.13,len(v)), v, s=14, facecolor='white',
                   edgecolor='black', linewidth=0.7, zorder=3)

    hi = max((c.mean()+c.std(ddof=1)/np.sqrt(len(c))), (t_.mean()+t_.std(ddof=1)/np.sqrt(len(t_))),
             c.max(), t_.max())
    y = hi*1.12
    ax.plot([0,0,0.8,0.8],[y,y*1.04,y*1.04,y], color='black', lw=0.9)
    ax.text(0.4, y*1.07, pfmt(p), ha='center', va='bottom', fontsize=8.5)
    ax.set_ylim(0, hi*1.34); ax.set_xlim(-0.6, 1.4)
    ax.set_xticks([0,0.8]); ax.set_xticklabels([ctrl, test])
    ax.set_ylabel('α-ketoglutarate\n(normalized to %s)' % ctrl)
    ax.spines[['top','right']].set_visible(False)
    ax.tick_params(width=0.9, length=3)

out=REPO_ROOT + '/manuscript bar plots/aKG_assay_barplots'
fig.savefig(out+'.png', dpi=400, bbox_inches='tight')
fig.savefig(out+'.svg', bbox_inches='tight')
fig.savefig(out+'.pdf', bbox_inches='tight')

st = pd.DataFrame(res, columns=['control','test','mean_control','mean_test','t','p_two_tailed'])
st['test_used']="Welch's unpaired two-tailed t-test (n=3/group)"
st.to_csv(REPO_ROOT + '/manuscript bar plots/aKG_assay_stats.csv', index=False)
print(st.to_string(index=False))
