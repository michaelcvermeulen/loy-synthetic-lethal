"""Second introduction variant: all six published Qi genomic Y-status classes."""
from pathlib import Path
import json
import hashlib
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

ROOT=Path(__file__).resolve().parents[3]
OUT=ROOT/'figures/chapter3'
DATA=OUT/'source_data/qi_tcga_loy_prevalence'
sys.path.insert(0,str(ROOT/'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style, condition_colours
font=apply_style(9)
core=condition_colours()
statuses=['LOY','rLOY','pLOY','Gain','Gain_Loss','WT']
labels={'LOY':'Complete LOY','rLOY':'Relative LOY (rLOY)',
        'pLOY':'Partial LOY (pLOY)','Gain':'Gain','Gain_Loss':'Gain + loss','WT':'WT'}
colours={'LOY':core['Y-depleted'],'rLOY':'#F39B7F','pLOY':'#EFC000',
         'Gain':'#00A087','Gain_Loss':'#8491B4','WT':core['Y-retaining']}
# Hatched mixed gain/loss remains distinguishable in greyscale.
hatches={'Gain_Loss':'////'}
workbook=DATA/'Qi_S1_Table.xlsx'
source=pd.read_excel(workbook,sheet_name='Table S1A')
called=source[source.Y_status.notna()].copy()
assert called.case_id.nunique()==len(called)==5014
assert set(called.Y_status)==set(statuses)
reference=pd.read_csv(DATA/'prevalence_by_tumour_type.csv')
order=reference.Tumor_type.tolist()
assert len(order)==28
counts=pd.crosstab(called.Cohort,called.Y_status).reindex(columns=statuses,fill_value=0)
counts['n_called']=counts[statuses].sum(axis=1)
assert np.array_equal(counts.loc[order,'n_called'],reference.n_called)
assert np.array_equal(counts.loc[order,'LOY'],reference.n_complete_loy)
overall=called.Y_status.value_counts().reindex(statuses)
assert overall.to_dict()=={'LOY':1390,'rLOY':114,'pLOY':206,'Gain':239,'Gain_Loss':77,'WT':2988}
counts.to_csv(DATA/'y_status_counts_all_cohorts.csv')
long=[]
for cohort in list(counts.index)+['All tumours']:
    n=int(counts.loc[cohort,'n_called']) if cohort!='All tumours' else len(called)
    for status in statuses:
        k=int(counts.loc[cohort,status]) if cohort!='All tumours' else int(overall[status])
        long.append(dict(cohort=cohort,status=status,n=k,denominator=n,percentage=100*k/n,
                         displayed=(cohort in order or cohort=='All tumours')))
pd.DataFrame(long).to_csv(DATA/'y_status_distribution_long.csv',index=False)
display=counts.loc[order,statuses].copy()
denominators=counts.loc[order,'n_called'].to_numpy()
percent=100*display.to_numpy()/denominators[:,None]
assert np.allclose(percent.sum(axis=1),100)
matrix=np.vstack([percent,100*overall.to_numpy()/len(called)])
assert np.allclose(matrix.sum(axis=1),100)
positions=np.r_[np.arange(28),29.2]
ticklabels=[f'{row.tumour_name} ({row.Tumor_type})' for row in reference.itertuples()]+['All tumours']
fig,ax=plt.subplots(figsize=(183/25.4,215/25.4))
fig.subplots_adjust(left=.49,right=.875,bottom=.10,top=.84)
left=np.zeros(len(matrix))
for j,status in enumerate(statuses):
    values=matrix[:,j]
    ax.barh(positions,values,left=left,height=.68,color=colours[status],
            edgecolor='white',linewidth=.35,hatch=hatches.get(status),zorder=2)
    for i,value in enumerate(values):
        if value>=8:
            ax.text(left[i]+value/2,positions[i],f'{value:.0f}',ha='center',va='center',
                    fontsize=7,color='white' if status in ['LOY','Gain','WT'] else '#222222',zorder=3)
    left+=values
assert np.allclose(left,100)
ax.axhline(28.15,xmin=0,xmax=1,color='#AAAAAA',lw=.6)
ax.set_xlim(0,100)
ax.set_ylim(30,-.7)
ax.set_yticks(positions,ticklabels,fontsize=8)
ax.get_yticklabels()[-1].set_weight('semibold')
ax.set_xticks([0,20,40,60,80,100])
ax.set_xlabel('Tumours by genomic Y-chromosome status (%)',fontsize=9,labelpad=8)
ax.tick_params(axis='y',length=0,pad=8)
ax.tick_params(axis='x',length=3,labelsize=8)
ax.spines['left'].set_visible(False)
for i,n in enumerate(list(denominators)+[len(called)]):
    ax.text(1.26,positions[i],f'{n:,}',transform=ax.get_yaxis_transform(),ha='right',va='center',
            fontsize=8,color='#444444',weight='semibold' if i==28 else 'normal',clip_on=False)
ax.text(1.26,1.025,'n',transform=ax.transAxes,ha='right',va='bottom',fontsize=8.5,weight='semibold')
fig.text(.035,.975,'a',fontsize=13,weight='semibold',va='top')
fig.text(.085,.975,'Y-chromosome status across TCGA cancer types',fontsize=11,weight='semibold',va='top')
fig.text(.085,.935,'5,014 male primary tumours | published genomic Y-status calls',fontsize=9.5,va='top')
handles=[Patch(facecolor=colours[s],edgecolor='white',hatch=hatches.get(s),label=labels[s]) for s in statuses]
fig.legend(handles=handles,ncol=3,loc='upper center',bbox_to_anchor=(.53,.908),
           fontsize=8.5,columnspacing=1.8,handlelength=1.7,handleheight=1.0,labelspacing=.7)
fig.text(.035,.025,'Source: Qi et al. (2023), Table S1A.',
         fontsize=7.3,color='#555555')
stem='ch03_main_intro_tcga_y_status_distribution'
for ext in ['pdf','svg','png']:fig.savefig(OUT/f'{stem}.{ext}',dpi=600)
fig.savefig(OUT/'previews'/f'{stem}.png',dpi=180)
# Greyscale review only; not a second scientific export.
from PIL import Image
Image.open(OUT/'previews'/f'{stem}.png').convert('L').save(OUT/'previews'/f'{stem}_greyscale.png')
plt.close(fig)
meta=dict(figure=stem,font=font,size_mm=[183,215],dpi=600,status_order=statuses,
          labels=labels,colours=colours,hatches=hatches,overall_counts=overall.to_dict(),
          source_doi='10.1016/j.cell.2023.06.006',source_sheet='Table S1A',
          source_sha256=hashlib.sha256(workbook.read_bytes()).hexdigest(),
          displayed_tumour_types=order,displayed_denominator=int(denominators.sum()),overall_denominator=5014,
          uncalled_excluded=27,segment_label_threshold_percent=8,
          interpretation='Six mutually exclusive source categories, not a reconstructed chronology. WT retains the published label.',
          script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
(DATA/'y_status_distribution_provenance.json').write_text(json.dumps(meta,indent=2),encoding='utf-8')
print('Rendered six-class stacked distributions for 28 cancer types and the complete 5,014-tumour cohort.')
