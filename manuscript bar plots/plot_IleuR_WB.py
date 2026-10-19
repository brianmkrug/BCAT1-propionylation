import os as _os; REPO_ROOT = _os.environ.get("BCAT1_REPO") or _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))  # publication: repo-relative root
import json, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

P=REPO_ROOT + '/'
pal = json.load(open(P+'BCAT1_figure_palette.json'))
G = pal['isoleucine_dose']
doses = [416, 70, 7]
C = {416: G['Ile_416uM']['hex'], 70: G['Ile_70uM']['hex'], 7: G['Ile_7uM']['hex']}

df = pd.read_excel(P+'manuscript bar plots/Source data.xlsx', sheet_name='IleuR WB').dropna(subset=['Cell line'])
df['pair'] = df['Cell line'] + ' r' + df['replicate'].astype(int).astype(str)
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['DejaVu Sans'],
                     'font.size':9,'axes.linewidth':0.9,'svg.fonttype':'none','pdf.fonttype':42,'ps.fonttype':42})
def pfmt(p): return f"p = {p:.3f}" if p >= 0.001 else f"p = {p:.1e}"

readouts = [('H3K9pr','H3K9pr\n(normalized to 416 µM)'),
            ('POL2-S2phos/POL2','POL2 S2phos / total POL2\n(normalized to 416 µM)')]
pairs = df['pair'].unique().tolist()
xs = [0, 0.8, 1.6]
fig, axes = plt.subplots(1, 2, figsize=(4.8, 3.1), gridspec_kw={'wspace':0.6})
res = []

for ax, (V, ylab) in zip(axes, readouts):
    vals = {d: df[df['Isoleucine (uM)']==d].set_index('pair')[V].reindex(pairs).values for d in doses}
    offs = np.linspace(-0.19, 0.19, len(pairs))
    for x, d in zip(xs, doses):
        v = vals[d]
        ax.bar(x, v.mean(), width=0.7, color=C[d], edgecolor='none',
               yerr=v.std(ddof=1)/np.sqrt(len(v)), capsize=3,
               error_kw=dict(elinewidth=0.9, capthick=0.9), zorder=2)
        ax.scatter(x+offs, v, s=14, facecolor='white', edgecolor='black',
                   linewidth=0.7, zorder=4)

    t, p = stats.ttest_rel(vals[416], vals[7])          # 416 vs 7 only
    res.append((V, vals[416].mean(), vals[70].mean(), vals[7].mean(), len(pairs), t, p))

    hi = max(np.concatenate([vals[d] for d in doses]))
    y = hi*1.12
    ax.plot([xs[0],xs[0],xs[2],xs[2]],[y,y*1.04,y*1.04,y], color='black', lw=0.9)
    ax.text((xs[0]+xs[2])/2, y*1.07, pfmt(p), ha='center', va='bottom', fontsize=8.5)
    ax.set_ylim(0, hi*1.34); ax.set_xlim(-0.6, 2.2)
    ax.set_xticks(xs); ax.set_xticklabels(['416','70','7'], fontsize=8.5)
    ax.set_xlabel('Isoleucine (µM)', fontsize=8.5)
    ax.set_ylabel(ylab)
    ax.spines[['top','right']].set_visible(False)
    ax.tick_params(width=0.9, length=3)

out=P+'manuscript bar plots/IleuR_WB_barplots'
fig.savefig(out+'.png', dpi=400, bbox_inches='tight')
fig.savefig(out+'.svg', bbox_inches='tight')
fig.savefig(out+'.pdf', bbox_inches='tight')

st = pd.DataFrame(res, columns=['readout','mean_416uM','mean_70uM','mean_7uM','n_pairs','t','p_416_vs_7'])
st['test']='paired two-tailed t-test, 416 uM vs 7 uM only; n=4 blots (2 lines x 2 replicates)'
st.to_csv(P+'manuscript bar plots/IleuR_WB_stats.csv', index=False)
print(st.to_string(index=False))
