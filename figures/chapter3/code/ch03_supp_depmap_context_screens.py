from pathlib import Path
import hashlib,json,shutil
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Patch
root=Path(__file__).resolve().parents[3]
out=root/'figures/chapter3'
source=root/'research/depmap_loy_reanalysis_v1/loy_reanalysis_v1_export_2026-10-07/results/loy_reanalysis_v1/figure_source_data'
data=out/'source_data/depmap_context_screens';data.mkdir(parents=True,exist_ok=True)
name='fig5c_context_mde.csv'
m=pd.read_csv(source/name)
assert len(m)==116
# Author decisions: omit genetic cancer-type-versus-rest tests and all drug interactions.
# This supplement now contains only the retained genetic context displays.
excluded_classes=['lineage','disease','subtype']
m=m.loc[~m['class'].isin(excluded_classes)].drop(columns=['interaction_fdr_global_hits']).copy()
assert len(m)==86 and m.tests.sum()==523654
m.to_csv(data/name,index=False)
for font in ['Regular','SemiBold']:font_manager.fontManager.addfont(str(root/f'skills/thesis-figure-style/assets/fonts/OpenSans-{font}.ttf'))
plt.rcParams.update({'font.family':'Open Sans','font.size':8.5,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.65,'pdf.fonttype':42,'svg.fonttype':'none'})
fig=plt.figure(figsize=(9.4,5.4))
black='#222222'
# a (previously c): retained genetic driver/genomic contexts, unchanged estimates.
ax=fig.add_axes([.255,.64,.695,.29])
order=m.groupby('class_label').median_interaction_mde.median().sort_values(ascending=False).index.tolist()
labels={'primary disease':'Primary disease','driver gene':'Driver gene','subtype':'Subtype','amplification':'Amplification','lineage':'Cancer lineage','molecular subtype':'Molecular subtype','genome doubling, MSI':'Genome doubling / MSI','deep deletion':'Deep deletion','CIN, ploidy, aneuploidy, LOH':'CIN / ploidy / aneuploidy / LOH'}
rng=np.random.default_rng(20261005)
for i,label in enumerate(order):
    vals=m[m.class_label==label].median_interaction_mde.to_numpy()
    ax.boxplot([vals],positions=[i],vert=False,widths=.58,showfliers=False,patch_artist=True,boxprops={'facecolor':'#F0F0F0','edgecolor':'#555555','lw':.7},medianprops={'color':black,'lw':1.1},whiskerprops={'color':'#555555','lw':.7},capprops={'color':'#555555','lw':.7},manage_ticks=False)
    ax.scatter(vals,i+rng.uniform(-.17,.17,len(vals)),s=8,c='#555555',alpha=.65,zorder=3)
ax.set(yticks=range(len(order)),yticklabels=[labels[x] for x in order],ylim=(len(order)-.5,-.5),xlabel='Median detectable interaction (Chronos units; 80% power, nominal P = 0.05)')
ax.tick_params(axis='y',length=0,pad=7)
fig.text(.02,.96,'a',fontsize=12,weight='semibold')
fig.text(.255,.96,'86 contexts; median detectable interaction = '+f'{m.median_interaction_mde.median():.2f}'+' Chronos units',fontsize=8)
# b (previously d): the approved CIN-stratified ARHGEF2 boxplot.
box_source=out/'source_data/arhgef2_cin_boxplot/plotted_points.csv'
shutil.copy2(box_source,data/'fig5d_arhgef2_cin_points.csv')
p=pd.read_csv(data/'fig5d_arhgef2_cin_points.csv')
assert len(p)==473 and set(p.assay)=={'CRISPR'}
cuts=p.cin.quantile([1/3,2/3]).to_numpy()
assert p.groupby('y_status').size().to_dict()=={'Y-negative':183,'Y-positive':290}
hits_source=root/'research/depmap_loy_reanalysis_v1/loy_reanalysis_v1_export_2026-10-07/derived/loy_reanalysis_v1/context_gene_hits.csv'
hit=pd.read_csv(hits_source).query("symbol == 'ARHGEF2' and context_id == 'genomic_continuous__CIN'")
assert len(hit)==1
hit.to_csv(data/'fig5d_arhgef2_discovery.csv',index=False)
ax=fig.add_axes([.105,.13,.50,.34])
colours={'Y-positive':'#3C5488','Y-negative':'#E64B35'}
rng=np.random.default_rng(20261007)
for i,context in enumerate(['Low CIN','Middle CIN','High CIN']):
    for j,(group,colour) in enumerate(colours.items()):
        vals=p.loc[(p.cin_group==context)&(p.y_status==group),'adjusted_dependency'].to_numpy()
        x=i+(-.19 if j==0 else .19);width=.29
        ax.scatter(x+rng.uniform(-width*.39,width*.39,len(vals)),vals,s=8,c=colour,alpha=.36,edgecolors='none',zorder=2)
        ax.boxplot([vals],positions=[x],widths=width,patch_artist=True,showfliers=False,manage_ticks=False,
            boxprops={'facecolor':colour,'edgecolor':colour,'alpha':.22,'linewidth':1},
            medianprops={'color':colour,'linewidth':1.5},whiskerprops={'color':colour,'linewidth':.8},
            capprops={'color':colour,'linewidth':.8},zorder=3)
        ax.text(x,3.07,f'n = {len(vals)}',ha='center',fontsize=7.5,color=colour)
ax.set(xlim=(-.55,2.55),ylim=(-3.2,3.5),xticks=[0,1,2],yticks=[-3,-2,-1,0,1,2,3],
    xticklabels=[f'Low CIN\n≤{cuts[0]:.3f}',f'Middle CIN\n>{cuts[0]:.3f}–{cuts[1]:.3f}',f'High CIN\n>{cuts[1]:.3f}'],
    ylabel='Adjusted ARHGEF2 gene effect\n(rank-normalised)')
ax.legend(handles=[Patch(facecolor=c,edgecolor=c,label=g,alpha=.65) for g,c in colours.items()],
    loc='lower right',bbox_to_anchor=(1,1.045),ncol=2,frameon=False,fontsize=8,handlelength=1.2,columnspacing=1.3)
fig.text(.02,.525,'b',fontsize=12,weight='semibold')
fig.text(.105,.525,f'LOY × CIN: class FDR = {hit.interaction_fdr_class.iloc[0]:.3f}',fontsize=8)
fig.text(.105,.022,'CIN thirds are for display only; lower values indicate stronger dependency.',fontsize=8,color='#555555')
# c (previously e): CIB1 within-cancer association, with raw Chronos measurements.
cib_source=out/'source_data/cib1_bladder_boxplot'
for original_name,new_name in [('plotted_points.csv','fig5e_cib1_points.csv'),('saved_within_cancer_result.csv','fig5e_cib1_discovery.csv')]:
    shutil.copy2(cib_source/original_name,data/new_name)
cib=pd.read_csv(data/'fig5e_cib1_points.csv')
cib_hit=pd.read_csv(data/'fig5e_cib1_discovery.csv').iloc[0]
assert cib.groupby('y_status').size().to_dict()=={'LOY':8,'Y-positive':10}
assert np.isclose(cib.loc[cib.y_status.eq('LOY'),'chronos'].mean()-cib.loc[cib.y_status.eq('Y-positive'),'chronos'].mean(),cib_hit.inside,atol=1e-9)
ax=fig.add_axes([.745,.13,.225,.34])
rng=np.random.default_rng(20261007)
for i,(group,colour) in enumerate({'Y-positive':'#3C5488','LOY':'#E64B35'}.items()):
    vals=cib.loc[cib.y_status.eq(group),'chronos'].to_numpy()
    ax.boxplot([vals],positions=[i],widths=.49,showfliers=False,manage_ticks=False,patch_artist=True,
        boxprops={'facecolor':colour,'edgecolor':colour,'alpha':.20,'linewidth':1},
        medianprops={'color':colour,'linewidth':1.6},whiskerprops={'color':colour,'linewidth':.9},
        capprops={'color':colour,'linewidth':.9})
    ax.scatter(i+rng.uniform(-.14,.14,len(vals)),vals,s=20,c=colour,alpha=.8,edgecolors='white',linewidth=.4,zorder=3)
ax.axhline(0,color='#999999',linestyle=(0,(3,3)),linewidth=.65,zorder=0)
ax.set(xlim=(-.58,1.58),ylim=(-.50,1.52),xticks=[0,1],xticklabels=['Y-positive\nn = 10','LOY\nn = 8'],
    yticks=[-.4,0,.4,.8,1.2],ylabel='CIB1 gene effect\n(Chronos)',xlabel='Bladder/urinary tract')
ax.plot([0,0,1,1],[1.28,1.33,1.33,1.28],color=black,linewidth=.75)
ax.text(.5,1.36,f'Within-cancer FDR = {cib_hit.own_fdr:.3f}',ha='center',va='bottom',fontsize=8)
fig.text(.66,.525,'c',fontsize=12,weight='semibold')
stem='ch03_supp_depmap_context_screens';code='CH3_S3_4_2_S02_depmap_context_screens'
selected=out/'selected'/code
for ext in ['pdf','png','svg']:
    f=out/f'{stem}.{ext}';fig.savefig(f,dpi=300,facecolor='white');shutil.copy2(f,selected/f'{code}.{ext}')
plt.close(fig)
(data/'provenance.json').write_text(json.dumps({
    'source_files':{'fig5c_context_mde.csv':hashlib.sha256((source/'fig5c_context_mde.csv').read_bytes()).hexdigest()},
    'panel_b_sources':{str(box_source):hashlib.sha256(box_source.read_bytes()).hexdigest(),str(hits_source):hashlib.sha256(hits_source.read_bytes()).hexdigest()},
    'panel_c_sources':{str(cib_source/name):hashlib.sha256((cib_source/name).read_bytes()).hexdigest() for name in ['plotted_points.csv','saved_within_cancer_result.csv']},
    'changes':'Author removed all drug interactions after the earlier removal of genetic cancer-type-versus-rest interactions. Former drug calibration panels a-b are removed; remaining genetic panels c-e are relabelled a-c. All genetic measurements, estimates and correction families are unchanged. Panel a: 86 driver/genomic contexts and 523654 tests. Panel b: 473 adjusted ARHGEF2 values and continuous-genomic class FDR. Panel c: 18 raw CIB1 Chronos scores and within-cancer FDR. No new fits or hypothesis tests.',
    'source_filename_mapping':{'a':'fig5c_context_mde.csv','b':'fig5d_arhgef2_cin_points.csv and fig5d_arhgef2_discovery.csv','c':'fig5e_cib1_points.csv and fig5e_cib1_discovery.csv'},
    'retired_sources':'fig5a_context_yield.csv and fig5b_context_tail.csv are historical drug-interaction tables, no longer loaded or displayed.',
    'excluded_genetic_context_classes':excluded_classes,
    'retained_genetic_contexts':int(len(m)),
    'retained_genetic_tests':int(m.tests.sum()),
    'panel_b_cin_cutoffs':cuts.tolist(),
    'panel_b_fdr_class':float(hit.interaction_fdr_class.iloc[0]),
    'panel_c_within_cancer_fdr':float(cib_hit.own_fdr),
    'panel_c_inference':'FDR from the original globally rank-normalised limma screen, BH across 6089 genes within bladder/urinary tract. Library covariates are constant in these 18 models. Scores shown are untransformed.',
    'original_screen_global_fdr_note':'The imported discovery row includes the historical global FDR across all 116 contexts. It is retained only as source provenance and is not used in the revised figure or chapter.',
    'font':'bundled Open Sans','code':str(Path(__file__).resolve())},indent=2),encoding='utf-8')
(selected/'selection_notes.md').write_text('Updated after the author removed drug interactions. The former drug-calibration panels a-b are removed. The retained genetic panels are now a (previously c): detection limits across 86 driver/genomic contexts; b (previously d): ARHGEF2 in 473 models across descriptive CIN thirds, class FDR = 0.038; c (previously e): raw CIB1 scores in 10 Y-positive and 8 LOY bladder/urinary tract models, within-cancer FDR = 0.030. All genetic values and correction families are unchanged. The source CSV filenames retain their historical figure/panel codes; the mapping is recorded in provenance.json. Original HPC outputs and superseded figures remain archived. Code: ../../code/ch03_supp_depmap_context_screens.py. Plotted data/provenance: ../../source_data/depmap_context_screens/. Caption: [central captions](../../../../Chapter3_figure_captions.md). No new fits or hypothesis tests.\n',encoding='utf-8')
print(out/f'{stem}.png')
