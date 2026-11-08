"""Current-baseline add-back overlap; no tests or transformations beyond display."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

figroot = Path(__file__).resolve().parents[1]
root = figroot.parents[1]
font_manager.fontManager.ttflist = [f for f in font_manager.fontManager.ttflist if f.name != 'Open Sans']
for p in (root/'skills/thesis-figure-style/assets/fonts').rglob('*.ttf'):
    if 'variable' in p.name:
        continue
    font_manager.fontManager.addfont(str(p))
plt.rcParams.update({'font.family':'Open Sans','font.size':8,'pdf.fonttype':42,'svg.fonttype':'none'})
d = pd.read_csv(figroot/'source_data/addback_atac_overlap/gene_table.csv')
d = d[(d['set']=='joint_FDR005') & d.tested_in_addback].sort_values('baseline_logFC')
keys = ['UTY_v_LXHA','KDM5D_v_LXHA','DUAL_v_LXHA','LOY_v_LXHA']
cols = ['baseline_logFC']+[k+'_logFC' for k in keys]
a = d[cols].to_numpy(); lim = np.abs(a).max()
fig, ax = plt.subplots(figsize=(6.4,9.2))
fig.subplots_adjust(left=.17,right=.80,bottom=.03,top=.77)
im = ax.imshow(a,aspect='auto',cmap=LinearSegmentedColormap.from_list('effects',['#2166ac','white','#d6604d']),norm=TwoSlopeNorm(0,-lim,lim),interpolation='none')
ax.set_yticks(range(len(d)),d.gene_name,fontstyle='italic',fontsize=7)
ax.set_xticks(range(5),['Y-depleted vs sgCTRL','UTY vs vector','KDM5D vs vector','UTY + KDM5D vs vector','Untransduced Y-depleted vs vector'],rotation=45,ha='left',rotation_mode='anchor',fontsize=7)
ax.xaxis.tick_top();ax.tick_params(length=0,pad=3)
for j,k in enumerate(keys,1):
    rows=np.flatnonzero(d[k+'_FDR'].to_numpy()<.05)
    ax.scatter(np.full(len(rows),j),rows,s=6,c='black',linewidths=0)
for spine in ax.spines.values():spine.set_visible(False)
cax=fig.add_axes([.85,.42,.022,.19]);cb=fig.colorbar(im,cax=cax);cb.set_label('log2 fold change',fontsize=8);cb.outline.set_visible(False)
for ext in ['pdf','svg','png']:
    fig.savefig(figroot/f'ch03_supp_ht1080_addback_atac_overlap.{ext}',dpi=300,facecolor='white')
plt.close(fig)
