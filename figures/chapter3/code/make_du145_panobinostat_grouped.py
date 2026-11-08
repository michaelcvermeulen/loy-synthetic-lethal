"""Independent nucleofection replicates grouped by panobinostat exposure."""
from pathlib import Path
import sys, json, hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[3]
FIG = ROOT / 'figures/chapter3'
SRC = FIG / 'source_data/du145_drug_optimisation'
sys.path.insert(0, str(ROOT/'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
apply_style(8)

conditions = ['DMSO vehicle', 'KU-55933', 'Panobinostat 30 nM',
              'Panobinostat 75 nM', 'KU-55933 + panobinostat 30 nM',
              'KU-55933 + panobinostat 75 nM']
labels = ['DMSO', 'KU', 'Pano 30 nM', 'Pano 75 nM', 'KU + pano 30 nM', 'KU + pano 75 nM']
colors = ['#777777', '#4DBBD5', '#E99B57', '#E64B35', '#00A087', '#8491B4']
markers = ['o', 'v', 's', 'D', 'P', 'X']
meta = {c:dict(label=l, color=k, marker=m, group=int(i>=2))
        for i,(c,l,k,m) in enumerate(zip(conditions,labels,colors,markers))}
dd = pd.read_csv(SRC/'all_target_measurements.csv')
wgs = pd.read_csv(SRC/'lpwgs/raw_depth_ratios.csv')
rows=[]
for target in ['AMELY','DDX3Y','KDM5D']:
    subset=dd[dd.target.eq(target)]
    control=float(subset.loc[subset.condition.eq('sgCTRL'),'concentration_ratio'].iloc[0])
    for c in conditions:
        value=float(subset.loc[subset.condition.eq(c),'concentration_ratio'].iloc[0])
        rows.append(dict(measure=target, condition=c, pano=bool(meta[c]['group']),
                         ratio=value, sgCTRL_ratio=control, residual_percent=100*value/control))
control=float(wgs.loc[wgs.condition.eq('sgCTRL'),'chrY_autosome_ratio'].iloc[0])
for c in conditions:
    value=float(wgs.loc[wgs.condition.eq(c),'chrY_autosome_ratio'].iloc[0])
    rows.append(dict(measure='lpWGS',condition=c,pano=bool(meta[c]['group']),
                     ratio=value,sgCTRL_ratio=control,residual_percent=100*value/control))
data=pd.DataFrame(rows)
assert len(data)==24 and not data.duplicated(['measure','condition']).any()
data.to_csv(SRC/'panobinostat_grouped_measurements.csv',index=False)
summary=data.groupby(['measure','pano']).agg(n_conditions=('ratio','size'),
    median_ratio=('ratio','median'),median_residual_percent=('residual_percent','median'),
    min_residual_percent=('residual_percent','min'),max_residual_percent=('residual_percent','max')).reset_index()
summary.to_csv(SRC/'panobinostat_grouped_summary.csv',index=False)

fig,axs=plt.subplots(2,2,figsize=(183/25.4,145/25.4),sharey=True)
fig.subplots_adjust(left=.12,right=.96,bottom=.23,top=.87,hspace=.64,wspace=.27)
for letter,measure,title,ax in zip('abcd',['AMELY','DDX3Y','KDM5D','lpWGS'],
                                  ['AMELY / AMELX','DDX3Y / RNaseP','KDM5D / RNaseP','lpWGS: whole-Y / autosomes'],axs.flat):
    ax.text(-.18,1.13,letter,transform=ax.transAxes,fontweight='semibold',fontsize=12)
    ax.set_title(title,fontsize=9,fontweight='semibold',pad=8)
    for group in [0,1]:
        cs=[c for c in conditions if meta[c]['group']==group]
        for c,offset in zip(cs,np.linspace(-.15,.15,len(cs))):
            row=data[data.measure.eq(measure)&data.condition.eq(c)].iloc[0]
            ax.scatter(group+offset,row.residual_percent,s=42,marker=meta[c]['marker'],
                       color=meta[c]['color'],edgecolors='white',linewidths=.5,zorder=3)
        med=float(summary.loc[summary.measure.eq(measure)&summary.pano.eq(bool(group)),
                              'median_residual_percent'].iloc[0])
        ax.plot([group-.22,group+.22],[med,med],color='#333333',lw=1.3,zorder=2)
        ax.text(group+.26,med,f'{med:.1f}%',va='center',ha='left',fontsize=7)
    ax.set(xlim=(-.4,1.65),ylim=(0,25),yticks=[0,5,10,15,20,25],xticks=[0,1],
           xticklabels=['No panobinostat\nn = 2','Panobinostat\nn = 4'])
    ax.tick_params(labelsize=7.5)
    ax.grid(axis='y',color='#EEEEEE',lw=.5,zorder=0)
    if measure in ['AMELY','KDM5D']: ax.set_ylabel('Residual signal (% of sgCTRL)',fontsize=8)
fig.text(.12,.965,'DU145 | one pooled CP + Q53 + P47 nucleofection',fontsize=10,fontweight='semibold')
fig.text(.12,.924,'Each point is an independent nucleofection replicate; horizontal lines show group medians.',fontsize=7.5)
handles=[Line2D([],[],marker=meta[c]['marker'],color=meta[c]['color'],ls='',markersize=6,label=meta[c]['label']) for c in conditions]
fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.53,.075),ncol=3,frameon=False,fontsize=7.5,columnspacing=2)
fig.text(.12,.043,'Biological replicates grouped by panobinostat exposure; colour and shape identify treatment.',fontsize=7,color='#555555')
fig.text(.12,.018,'ddPCR and lpWGS are linked by treatment condition, not confirmed matched DNA samples.',fontsize=6.7,color='#555555')
stem='ch03_supp_du145_panobinostat_grouped'
for ext in ['pdf','svg','png']: fig.savefig(FIG/f'{stem}.{ext}',dpi=600,facecolor='white')
fig.savefig(FIG/'previews'/f'{stem}.png',dpi=160,facecolor='white')
plt.close(fig)
record=dict(grouping={'no_panobinostat':conditions[:2],'panobinostat':conditions[2:]},
    normalization='100 * condition ratio / sgCTRL ratio, separately by assay; same-run ddPCR controls',
    unit='Independent biological nucleofection replicate; one per treatment condition, author-confirmed 23 September 2026',
    excluded=['Two-round population','LOY pass1 source label'],inferential_tests=False,
    inputs={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in
            [SRC/'all_target_measurements.csv',SRC/'lpwgs/raw_depth_ratios.csv']},
    outputs={ext:hashlib.sha256((FIG/f'{stem}.{ext}').read_bytes()).hexdigest() for ext in ['pdf','svg','png']})
(SRC/'panobinostat_grouped_provenance.json').write_text(json.dumps(record,indent=2),encoding='utf-8')
print(summary.to_string(index=False))
