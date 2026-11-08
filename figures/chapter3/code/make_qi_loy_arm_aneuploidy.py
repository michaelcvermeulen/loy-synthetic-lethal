"""Reprocess Qi Fig. 2E and restyle published Fig. 2F GISTIC summaries."""
from pathlib import Path
import json
import hashlib
import sys
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from adjustText import adjust_text

ROOT=Path(__file__).resolve().parents[3]
OUT=ROOT/'figures/chapter3'
DATA=OUT/'source_data/qi_loy_arm_aneuploidy'
sys.path.insert(0,str(ROOT/'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style, condition_colours
font=apply_style(9)
highlight=condition_colours()['Y-depleted']
focus={'KIRP','KIRC','UVM'}
ref=pd.read_csv(OUT/'source_data/qi_tcga_loy_prevalence/prevalence_by_tumour_type.csv')
clinical=pd.read_excel(DATA/'TCGA-CDR-SupplementalTableS1.xlsx',sheet_name='TCGA-CDR',
                       usecols=['bcr_patient_barcode','type','gender'])
assert clinical.bcr_patient_barcode.is_unique
arms=pd.read_excel(DATA/'Taylor_2018_Table_S2.xlsx',skiprows=1)
arms['case_id']=arms.Sample.str[:12]
# Match the authors' notebook: all male Taylor sample records are used in
# each cancer-type mean, irrespective of availability of a Qi Y-status call.
# Multiple sample records per case are retained, as in its .loc lookup.
arms=arms.merge(clinical,left_on='case_id',right_on='bcr_patient_barcode',
                how='left',validate='many_to_one',indicator=True)
unmatched=arms[arms._merge!='both'][['Sample','case_id','Type','AS_del']]
unmatched.to_csv(DATA/'arm_records_without_cdr_match.csv',index=False)
assert arms.loc[arms._merge=='both','Type'].eq(arms.loc[arms._merge=='both','type']).all()
male=arms[arms.gender.eq('MALE')].copy()
assert len(male)==5039 and male.case_id.nunique()==5005
assert male.AS_del.notna().all()
assert np.allclose(male['AneuploidyScore(AS)'],male.AS_del+male.AS_amp)
male[['Sample','case_id','type','gender','AS_del','AS_amp','AneuploidyScore(AS)']].to_csv(
    DATA/'male_autosomal_arm_scores.csv',index=False)
male[male.case_id.duplicated(False)][['Sample','case_id','type','AS_del']].to_csv(
    DATA/'multiple_sample_records_per_case.csv',index=False)
summary=male.groupby('type').agg(mean_arm_losses=('AS_del','mean'),
    n_arm_loss=('AS_del','count'),n_arm_cases=('case_id','nunique'))
plot_a=ref.merge(summary,left_on='Tumor_type',right_index=True,how='left',validate='one_to_one')
assert len(plot_a)==28 and plot_a.mean_arm_losses.notna().all()
result=pearsonr(plot_a.percentage_complete_loy,plot_a.mean_arm_losses)
assert round(float(result.statistic),2)==.33 and round(float(result.pvalue),4)==.0833
plot_a.to_csv(DATA/'prevalence_vs_mean_arm_loss.csv',index=False)
gistic=pd.read_excel(DATA/'Qi_Table_S2.xlsx',sheet_name='Table S2A')
plot_b=gistic[gistic.Arm.eq('Yp') & gistic.Tumor_type.isin(ref.Tumor_type)].copy()
assert len(plot_b)==28 and plot_b.Tumor_type.is_unique
plot_b=plot_b.merge(ref[['Tumor_type','tumour_name']],on='Tumor_type',validate='one_to_one')
plot_b['yp_deletion_percentage']=100*plot_b['Del frequency']
plot_b.to_csv(DATA/'yp_gistic_published.csv',index=False)

fig,axes=plt.subplots(1,2,figsize=(240/25.4,145/25.4))
fig.subplots_adjust(left=.085,right=.97,bottom=.19,top=.79,wspace=.30)
specs=[(plot_a,'percentage_complete_loy','mean_arm_losses'),
       (plot_b,'yp_deletion_percentage','Del z-score')]
for ax,(data,x,y) in zip(axes,specs):
    ax.set_xlim(-3,91)
    ax.set_xticks([0,20,40,60,80])
    if ax is axes[0]:
        ax.set_ylim(-.6,15.3);ax.set_yticks([0,2,4,6,8,10,12,14])
    else:
        ax.set_ylim(-4,30);ax.set_yticks([0,5,10,15,20,25])
        ax.axhline(0,color='#DDDDDD',lw=.6,zorder=0)
    labels=[]
    for row in data.to_dict('records'):
        code=row['Tumor_type'];marked=code in focus
        colour=highlight if marked else '#747B82'
        ax.scatter(row[x],row[y],s=30 if marked else 23,color=colour,
                   edgecolor='white',linewidth=.55,zorder=4 if marked else 3)
        labels.append(ax.text(row[x]+1.0,row[y]+(.15 if ax is axes[0] else .3),code,
            fontsize=8,weight='semibold' if marked else 'normal',
            color=highlight if marked else '#444444',zorder=5))
    np.random.seed(37)
    adjust_text(labels,x=data[x].to_numpy(),y=data[y].to_numpy(),ax=ax,
        expand=(1.12,1.2),force_text=(.3,.4),force_static=(.25,.3),
        ensure_inside_axes=True,iter_lim=800,
        arrowprops=dict(arrowstyle='-',color='#AAAAAA',lw=.45))
    ax.tick_params(labelsize=8)
axes[0].set_xlabel('Male tumours with complete LOY (%)',labelpad=8)
axes[0].set_ylabel('Mean number of autosomal arms lost',labelpad=8)
axes[1].set_xlabel('Male tumours with Yp deletion (%)',labelpad=8)
axes[1].set_ylabel('GISTIC Yp deletion z-score',labelpad=8)
axes[0].text(.02,.975,f'Pearson r = {result.statistic:.3f}; P = {result.pvalue:.4f}\nn = 28 cancer types',
             transform=axes[0].transAxes,va='top',fontsize=8.5)
fig.text(.03,.96,'a',fontsize=13,weight='semibold',va='top')
fig.text(.085,.96,'LOY and average arm loss',fontsize=11,weight='semibold',va='top')
fig.text(.53,.96,'b',fontsize=13,weight='semibold',va='top')
fig.text(.585,.96,'Yp loss under the background model',fontsize=11,weight='semibold',va='top')
fig.text(.085,.913,'Reprocessed from Qi Fig. 2E source data',fontsize=9,va='top')
fig.text(.585,.913,'Published GISTIC scores from Qi Fig. 2F / Table S2A',fontsize=8.5,va='top')
fig.legend(handles=[Line2D([],[],marker='o',ls='',color=highlight,
    label='KIRP, KIRC and UVM (highlighted in the text)')],
    loc='upper center',bbox_to_anchor=(.51,.872),frameon=False,fontsize=8.5)
fig.text(.03,.045,'Source: Qi et al. (2023), Tables S1A/S1C and S2A; Taylor et al. (2018), Table S2; TCGA-CDR.',
         fontsize=7.2,color='#555555')
stem='ch03_main_intro_loy_arm_aneuploidy'
for ext in ['pdf','svg','png']:fig.savefig(OUT/f'{stem}.{ext}',dpi=600)
fig.savefig(OUT/'previews'/f'{stem}.png',dpi=180)
plt.close(fig)
urls={
 'Qi_Table_S2.xlsx':'https://ars.els-cdn.com/content/image/1-s2.0-S0092867423006463-mmc2.xlsx',
 'Taylor_2018_Table_S2.xlsx':'https://ars.els-cdn.com/content/image/1-s2.0-S1535610818301119-mmc2.xlsx',
 'TCGA-CDR-SupplementalTableS1.xlsx':'https://ars.els-cdn.com/content/image/1-s2.0-S0092867418302290-mmc1.xlsx',
 'LOY GISTIC FINAL.ipynb':'https://raw.githubusercontent.com/rheinbaylab/LOY_in_primary_tumors/master/Notebooks/LOY%20GISTIC%20FINAL.ipynb',
 'LOY_LOX WEX Final analysis.ipynb':'https://raw.githubusercontent.com/rheinbaylab/LOY_in_primary_tumors/master/Notebooks/LOY_LOX%20WEX%20Final%20analysis.ipynb'}
meta=dict(figure=stem,font=font,size_mm=[240,145],dpi=600,
  pearson_r=float(result.statistic),pearson_p_two_sided=float(result.pvalue),n_cancer_types=28,
  mean_arm_loss_population='All male Taylor sample records linked to TCGA-CDR, as in Qi notebook; not restricted to Qi-called cases',
  male_arm_sample_records=len(male),male_arm_unique_cases=male.case_id.nunique(),
  arm_records_without_cdr_match=len(unmatched),
  male_arm_sample_records_plotted=int(plot_a.n_arm_loss.sum()),
  endpoint='AS_del: losses only, not combined gain-plus-loss aneuploidy score',
  highlights=sorted(focus),gistic='Published Yp arm scores and frequencies retained; not recomputed or equated with complete LOY',
  regression='Original robust regression curves omitted; no new trend fit or background model',
  missing_fig2g='Autosomal broad_significance_results.txt files not found locally or in public GitHub; per-arm gene-count plots not reconstructed',
  sources={k:dict(url=u,sha256=hashlib.sha256((DATA/k).read_bytes()).hexdigest()) for k,u in urls.items()},
  prevalence_source_sha256=hashlib.sha256((OUT/'source_data/qi_tcga_loy_prevalence/prevalence_by_tumour_type.csv').read_bytes()).hexdigest(),
  script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
(DATA/'render_provenance.json').write_text(json.dumps(meta,indent=2),encoding='utf-8')
print(json.dumps({k:meta[k] for k in ['pearson_r','pearson_p_two_sided','n_cancer_types','male_arm_sample_records_plotted']},indent=2))
