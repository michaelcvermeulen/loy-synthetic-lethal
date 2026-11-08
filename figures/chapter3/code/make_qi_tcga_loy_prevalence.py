"""Introduction panel: published complete genomic LOY by TCGA tumour type.

Use Qi et al. Table S1C, independently verify against individual S1A calls.
This is a descriptive redraw of published data, not a new LOY caller.
"""
from pathlib import Path
import sys
import json
import hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[3]
OUT=ROOT/'figures/chapter3'
DATA=OUT/'source_data/qi_tcga_loy_prevalence'
sys.path.insert(0,str(ROOT/'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style, condition_colours
font=apply_style(9)
colour=condition_colours()['Y-depleted']
workbook=DATA/'Qi_S1_Table.xlsx'
cases=pd.read_excel(workbook,sheet_name='Table S1A')
published=pd.read_excel(workbook,sheet_name='Table S1C')
assert len(cases)==5041 and cases.case_id.nunique()==len(cases)
assert cases.Gender.eq('male').all()
called=cases[cases.Y_status.notna()].copy()
assert len(called)==5014
assert called.Y_status.eq('LOY').sum()==1390
counts=called.groupby('Cohort').agg(
    n_called=('case_id','size'), n_complete_loy=('Y_status',lambda x:x.eq('LOY').sum()),
    n_relative_loy=('Y_status',lambda x:x.eq('rLOY').sum()),
    n_partial_loy=('Y_status',lambda x:x.eq('pLOY').sum()))
table=published[published.Male_sample_count.notna() & published.Fraction_LOY.notna()].copy()
table=table.merge(counts,left_on='Tumor_type',right_index=True,validate='one_to_one')
assert len(table)==28 and table.n_called.sum()==5003
assert np.array_equal(table.Male_sample_count.to_numpy(),table.n_called.to_numpy())
assert np.allclose(table.Fraction_LOY,table.n_complete_loy/table.n_called,atol=1e-12,rtol=0)
assert table.n_complete_loy.sum()==1390
omitted=counts.loc[~counts.index.isin(table.Tumor_type)]
assert omitted.index.tolist()==['BRCA']
assert omitted.loc['BRCA','n_called']==11 and omitted.loc['BRCA','n_complete_loy']==0
names={
 'ACC':'Adrenocortical carcinoma','BLCA':'Bladder urothelial carcinoma',
 'CHOL':'Cholangiocarcinoma','COAD':'Colon adenocarcinoma',
 'DLBC':'Diffuse large B-cell lymphoma','ESCA':'Oesophageal carcinoma',
 'GBM':'Glioblastoma','HNSC':'Head and neck squamous cell carcinoma',
 'KICH':'Chromophobe renal cell carcinoma','KIRC':'Clear cell renal cell carcinoma',
 'KIRP':'Papillary renal cell carcinoma','LAML':'Acute myeloid leukaemia',
 'LGG':'Lower-grade glioma','LIHC':'Hepatocellular carcinoma',
 'LUAD':'Lung adenocarcinoma','LUSC':'Lung squamous cell carcinoma',
 'MESO':'Mesothelioma','PAAD':'Pancreatic adenocarcinoma',
 'PCPG':'Phaeochromocytoma / paraganglioma','PRAD':'Prostate adenocarcinoma',
 'READ':'Rectal adenocarcinoma','SARC':'Sarcoma','SKCM':'Cutaneous melanoma',
 'STAD':'Gastric adenocarcinoma','TGCT':'Testicular germ cell tumours',
 'THCA':'Thyroid carcinoma','THYM':'Thymoma','UVM':'Uveal melanoma'}
table['tumour_name']=table.Tumor_type.map(names)
assert table.tumour_name.notna().all()
table['percentage_complete_loy']=100*table.n_complete_loy/table.n_called
table=table.sort_values(['percentage_complete_loy','Tumor_type'],ascending=[False,True]).reset_index(drop=True)
table.to_csv(DATA/'prevalence_by_tumour_type.csv',index=False)
counts.to_csv(DATA/'s1a_called_counts_all_cohorts.csv')
published.to_csv(DATA/'table_s1c_original_values.csv',index=False)
cases.loc[cases.Y_status.isna(),['case_id','Cohort','Y_status']].to_csv(DATA/'s1a_uncalled_cases.csv',index=False)

fig,ax=plt.subplots(figsize=(183/25.4,205/25.4))
fig.subplots_adjust(left=.49,right=.855,bottom=.095,top=.855)
y=np.arange(len(table))
ax.barh(y,table.percentage_complete_loy,height=.64,color=colour,linewidth=0,zorder=2)
ax.set_yticks(y,[f'{row.tumour_name} ({row.Tumor_type})' for row in table.itertuples()],fontsize=8)
ax.invert_yaxis()
ax.set_ylim(len(table)-.35,-.65)
ax.set_xlim(0,100)
ax.set_xticks([0,20,40,60,80,100])
ax.set_xlabel('Tumours with complete LOY (%)',fontsize=9,labelpad=8)
ax.tick_params(axis='y',length=0,pad=8)
ax.tick_params(axis='x',length=3,labelsize=8)
ax.spines['left'].set_visible(False)
for x in [20,40,60,80,100]:ax.axvline(x,color='#EEEEEE',lw=.6,zorder=0)
for i,row in enumerate(table.itertuples()):
    ax.text(row.percentage_complete_loy+1.8,i,f'{row.percentage_complete_loy:.1f}%',
            ha='left',va='center',fontsize=8,color='#333333')
    ax.text(1.34,i,f'{row.n_complete_loy:,} / {row.n_called:,}',transform=ax.get_yaxis_transform(),
            ha='right',va='center',fontsize=8,color='#444444',clip_on=False)
ax.text(1.34,1.027,'LOY / n',transform=ax.transAxes,ha='right',va='bottom',fontsize=8.5,weight='semibold')
fig.text(.035,.97,'a',fontsize=13,weight='semibold',va='top')
fig.text(.085,.97,'Complete LOY across TCGA cancer types',fontsize=11,weight='semibold',va='top')
fig.text(.085,.923,'Overall: 27.7% (1,390 / 5,014 male primary tumours)',fontsize=10,va='top')
fig.text(.085,.890,'Published cancer-type frequencies from Qi et al. (2023)',fontsize=8.5,color='#555555',va='top')
fig.text(.035,.018,'Complete genomic LOY; excludes relative and partial LOY. Source: Table S1C, checked against Table S1A.',
         fontsize=7.2,color='#555555')
stem='ch03_main_intro_tcga_loy_prevalence'
for ext in ['pdf','svg','png']:fig.savefig(OUT/f'{stem}.{ext}',dpi=600)
fig.savefig(OUT/'previews'/f'{stem}.png',dpi=180)
plt.close(fig)
meta=dict(figure=stem,font=font,size_mm=[183,205],dpi=600,colour=colour,
          source_doi='10.1016/j.cell.2023.06.006',source_sheets=['Table S1C','Table S1A'],
          workbook_sha256=hashlib.sha256(workbook.read_bytes()).hexdigest(),
          overall=dict(complete_loy=1390,called_male_tumours=5014,percentage=1390/5014*100),
          plotted=dict(cancer_types=28,complete_loy=1390,called_male_tumours=5003),
          omitted_from_breakdown='Table S1C omits a male BRCA prevalence; S1A has 11 called male BRCA tumours, none complete LOY. These remain in the published overall denominator.',
          uncalled=27,definition='Y_status exactly LOY; not rLOY, pLOY, Gain_Loss or expression-defined fLOY.',
          statistics='Descriptive published-cohort prevalence; no inferential tests or confidence intervals added.',
          script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
(DATA/'render_provenance.json').write_text(json.dumps(meta,indent=2),encoding='utf-8')
print('Rendered 28 cancer-type prevalences; all counts and fractions match S1A and S1C.')
