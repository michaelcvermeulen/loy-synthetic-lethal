from pathlib import Path
import hashlib
import json
import shutil
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle
from matplotlib.colors import LinearSegmentedColormap
R=Path(__file__).resolve().parents[3];F=R/'figures/chapter3';D=F/'source_data/du145_arm_recurrence'
sys.path.insert(0,str(R/'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
apply_style(7)
t=pd.read_csv(D/'arm_DNA_vs_RNA.csv')
samples=['LNCaP-LOY-2X','LNCaP-DMSO-LOY-1X','LNCaP-30','LNCaP-70','LNCaP-KU-30','LNCaP-KU-70','LNCaP-KU']
labels=['Guides 2×','Guides 1× + DMSO','Guides 1× + Pano 30 nM','Guides 1× + Pano 75 nM','Guides 1× + KU + Pano 30 nM','Guides 1× + KU + Pano 75 nM','Guides 1× + KU (low depth)','RNA: guides 2× vs sgCTRL']
arms=['chr'+str(c)+a for c in list(range(1,23))+['X','Y'] for a in ['p','q']]
mat=t.pivot(index='sample',columns='arm',values='DNA_median').reindex(index=samples,columns=arms)
rna=t[t['sample']==samples[0]].set_index('arm').RNA_mean_difference.reindex(arms)
values=np.vstack([mat.to_numpy(),rna.to_numpy()])
# Display-only exclusion: retain the source measurements and recurrence tables.
excluded_arms=['chr13p','chr14p','chr15p','chr21p','chr22p']
excluded_indices=[arms.index(arm) for arm in excluded_arms]
values[:,excluded_indices]=np.nan
pd.DataFrame(values,index=labels,columns=arms).to_csv(D/'displayed_matrix.csv')
(D/'display_mask.json').write_text(json.dumps({
    'excluded_arms':excluded_arms,
    'reason':'Acrocentric short arms (excluded / poorly resolved)',
    'scope':'Display only; all eight rows; source measurements and recurrence tables unchanged.',
    'source_table':'arm_DNA_vs_RNA.csv',
    'source_sha256':hashlib.sha256((D/'arm_DNA_vs_RNA.csv').read_bytes()).hexdigest(),
    'colour_limits':[-0.8,0.8],
    'style':{'facecolour':'#EEEEEE','hatch':'///','hatch_colour':'#BBBBBB'}
},indent=2)+'\n',encoding='utf-8')
z=t[(t.arm=='chr20p')&t['sample'].ne('LNCaP-KU')]
assert len(z)==6 and (z.DNA_median<-.2).all()
print(z[['sample','DNA_median','loss_fraction']].to_string(index=False))
fig,ax=plt.subplots(figsize=(210/25.4,90/25.4));fig.subplots_adjust(left=.265,right=.89,bottom=.16,top=.98)
cm=LinearSegmentedColormap.from_list('cn',['#3C5488','white','#E64B35']);cm.set_bad('#EEEEEE')
im=ax.imshow(values,aspect='auto',cmap=cm,vmin=-.8,vmax=.8,interpolation='nearest')
for i in excluded_indices:
    ax.add_patch(Rectangle((i-.5,-.5),1,len(labels),facecolor='#EEEEEE',
                           edgecolor='#BBBBBB',hatch='///',lw=0))
ax.set_yticks(range(8),labels,fontsize=6.5);ax.set_xticks(range(48),[x[3:] for x in arms],rotation=90,fontsize=5.8)
ax.tick_params(length=0,pad=3);ax.axhline(6.5,color='.3',lw=.9)
for x in np.arange(1.5,48,2):ax.axvline(x,color='white',lw=.4)
i=arms.index('chr20p');ax.add_patch(Rectangle((i-.5,-.5),1,8,fill=False,edgecolor='black',lw=.9))
ax.get_xticklabels()[i].set_fontweight('semibold')
for sp in ax.spines.values():sp.set_visible(False)
cb=fig.colorbar(im,cax=fig.add_axes([.905,.29,.013,.56]),extend='both');cb.set_label('Difference from sgCTRL (log₂)',fontsize=6.5);cb.ax.tick_params(labelsize=6);cb.outline.set_visible(False)
fig.legend(handles=[
    Patch(facecolor='#EEEEEE',edgecolor='#BBBBBB',hatch='///',lw=0,
          label='Acrocentric short arms (excluded / poorly resolved)'),
    Patch(facecolor='#EEEEEE',edgecolor='none',label='Unavailable RNA')
],loc='lower left',bbox_to_anchor=(.265,.015),ncol=2,fontsize=6,
   frameon=False,borderaxespad=0,handlelength=1.5,columnspacing=1.5)
stem='ch03_supp_du145_arm_recurrence'
selected_stem='CH3_S3_4_5_5_S03_du145_arm_recurrence'
selected=F/'selected'/selected_stem
selected.mkdir(parents=True,exist_ok=True)
for ext in ['pdf','svg','png']:
    working=F/(stem+'.'+ext)
    fig.savefig(working,dpi=300)
    shutil.copy2(working,selected/(selected_stem+'.'+ext))
plt.close(fig)
