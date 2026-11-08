from pathlib import Path
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
root=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(root.parents[1]/'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
apply_style(8)
source=root/'source_data/a549_hallmark_ssgsea'
d=pd.read_csv(source/'A549_ssGSEA_scores.csv')
g=d.groupby('pathway').ssGSEA
d['z']=(d.ssGSEA-g.transform('mean'))/g.transform('std')
assert np.isfinite(d.z).all()
meta=d.drop_duplicates('sample_key').sort_values(['treatment','sample_index'])
t=d.pivot(index='pathway',columns='sample_key',values='z').reindex(columns=meta.sample_key).sort_index()
assert t.shape==(50,8)
assert np.allclose(t.mean(axis=1),0,atol=1e-10) and np.allclose(t.std(axis=1),1)
t.to_csv(source/'A549_Hallmark_zscores.csv')
fig=plt.figure(figsize=(7.1,10))
ax=fig.add_axes([.40,.08,.46,.85]);cax=fig.add_axes([.89,.30,.018,.40])
cmap=LinearSegmentedColormap.from_list('diverging',['#3C5488','white','#E64B35'])
limit=np.ceil(np.abs(t.to_numpy()).max()*10)/10
im=ax.imshow(t,aspect='auto',cmap=cmap,norm=TwoSlopeNorm(vmin=-limit,vcenter=0,vmax=limit))
ax.set_yticks(range(50),[x.removeprefix('HALLMARK_').replace('_',' ') for x in t.index],fontsize=6.5)
ax.set_xticks(range(8),[str(int(x)) for x in meta.sample_index],fontsize=8)
ax.set_xlabel('Sample')
ax.axvline(3.5,color='white',lw=2)
ax.text(.25,1.015,'sgCTRL',color='#3C5488',transform=ax.transAxes,ha='center')
ax.text(.75,1.015,'Y-depleted',color='#E64B35',transform=ax.transAxes,ha='center')
ax.tick_params(axis='y',length=0)
fig.colorbar(im,cax=cax,label='ssGSEA row z-score')
for ext in ['pdf','svg','png']:
 fig.savefig(root/f'ch03_supp_a549_rna_hallmark_ssgsea.{ext}',dpi=300)
plt.close(fig)
print('Exported 50 Hallmark sets across eight samples; z-scores verified.')
