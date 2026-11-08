"""H2452 initial guide comparison restricted to the verified 3 December assay."""
from pathlib import Path
import sys, json, hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[3]
FIG=ROOT/'figures/chapter3'
SRC=FIG/'source_data/h2452_ddpcr_inventory'
sys.path.insert(0,str(ROOT/'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
apply_style(8)
p=SRC/'h2452_03Dec2024_verified_concentration_ratios.csv'
d=pd.read_csv(p)
assert np.allclose(d.AMELY_AMELX_ratio,d.AMELY_copies_per_uL/d.AMELX_copies_per_uL)
assert np.allclose(d.AMELY_AMELX_ratio,d.instrument_ratio)
fig,ax=plt.subplots(figsize=(120/25.4,100/25.4))
fig.subplots_adjust(left=.20,right=.95,bottom=.22,top=.83)
summary=[]
for i,g in enumerate(['CTRL2','CP','CN','Q53']):
    x=d[d.guide.eq(g)];values=x.AMELY_AMELX_ratio.to_numpy();med=np.median(values)
    color='#3C5488' if i==0 else '#E64B35'
    for j,row in enumerate(x.itertuples()):
        off=0 if len(x)==1 else (-.10 if j==0 else .10)
        ax.scatter(i+off,row.AMELY_AMELX_ratio,s=35,color=color,
                   marker='o' if row.DNA_date=='2024-11-21' else '^',zorder=3)
    ax.plot([i-.23,i+.23],[med,med],color=color,lw=2)
    ax.text(i,max(values)+.09,f'{med:.3f}\nn={len(x)}',ha='center',fontsize=7)
    summary.append(dict(condition='sgCTRL' if i==0 else g,n_measurements=len(x),median=med))
ax.set(xticks=range(4),xticklabels=['sgCTRL','CP','CN','Q53'],ylim=(0,1.65),
       ylabel='AMELY / AMELX concentration ratio',title='H2452 | initial individual guides')
ax.grid(axis='y',color='#EEEEEE',lw=.5);ax.set_axisbelow(True)
fig.text(.20,.10,'Points: DNA measurements; lines: medians.',fontsize=7)
fig.text(.20,.05,'Verified 3 December 2024 assay; no inferential test.',fontsize=6.5,color='#555555')
stem='ch03_supp_h2452_ddpcr_verified_initial_guides'
for ext in ['pdf','svg','png']:fig.savefig(FIG/f'{stem}.{ext}',dpi=600,facecolor='white')
fig.savefig(FIG/'previews'/f'{stem}.png',dpi=170,facecolor='white')
plt.close(fig)
pd.DataFrame(summary).to_csv(SRC/'verified_figure_summary.csv',index=False)
(SRC/'verified_figure_provenance.json').write_text(json.dumps(dict(input=str(p),
    sha256=hashlib.sha256(p.read_bytes()).hexdigest(),n_measurements=len(d),
    reason='Exclude unresolved cross-batch channel orientation and legacy droplet-ratio conventions',
    denominator='AMELX concentration',normalization='No reference copy-number scaling',
    independent_nucleofection_count='not assigned from measurement counts'),indent=2),encoding='utf-8')
print(pd.DataFrame(summary).to_string(index=False))
