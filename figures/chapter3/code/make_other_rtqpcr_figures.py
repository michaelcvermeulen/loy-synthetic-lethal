"""Process available additional LOY RT-qPCR tables without imputing nondetects."""
from pathlib import Path
import sys
import json
import math
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.patches import Rectangle, Polygon

ROOT=Path(__file__).resolve().parents[3]
OUT=ROOT/'figures/chapter3'
DATA=OUT/'source_data/other_rtqpcr'
sys.path.insert(0,str(ROOT/'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
apply_style(9.5)
CONFIG=[
 ('h2452',1,'RT_QPCR_H2452_ATTEMPT1.csv','H2452 WT'),
 ('ht1080_cp2',1,'RT_QPCR_HT1080_TRANSFECTION2_ATTEMPT1_.csv','HT1080 WT'),
 ('ht1080_cp2',2,'RT_QPCR_HT1080_TRANSFECTION2_ATTEMPT1_2PM.csv','HT1080 WT'),
 ('ht1080_clones',1,'RT_QPCR_HT1080_CLONES_20_DEC_2024.csv','C7')]
rows=[]
wells=[]
qc=[]
for dataset,run,name,reference in CONFIG:
    raw=pd.read_csv(DATA/name)
    raw['source_row']=np.arange(len(raw))+2
    raw['ct_numeric']=pd.to_numeric(raw.CT,errors='coerce')
    raw['mean_numeric']=pd.to_numeric(raw['Ct Mean'],errors='coerce')
    raw=raw[raw['Sample Name'].notna() & raw['Target Name'].notna()].copy()
    if 'Omit' in raw:
        assert not raw.Omit.astype(str).str.lower().eq('true').any()
    # Harmonise only the obvious spelling variant; retain original label.
    raw['gene']=raw['Target Name'].replace({'E1F1AY':'EIF1AY'})
    means=raw.groupby(['Sample Name','gene']).mean_numeric.first()
    for (sample,gene),group in raw.groupby(['Sample Name','gene']):
        assert sample=='NTC' or group.mean_numeric.dropna().nunique()<=1
        mean=means[(sample,gene)]
        numeric=group.ct_numeric.dropna()
        n=len(group)
        valid=len(numeric)
        literal40=int(group.ct_numeric.eq(40).sum())
        mismatch=mean-numeric.mean() if valid and np.isfinite(mean) else np.nan
        qc.append(dict(dataset=dataset,run=run,sample=sample,gene=gene,n_wells=n,
                       n_detected=valid,reported_mean=mean,recomputed_mean=numeric.mean(),
                       mean_difference=mismatch,literal_ct40_count=literal40,
                       n_distinct_reported_means=group.mean_numeric.dropna().nunique()))
        if sample=='NTC':
            continue
        gapdh=means.get((sample,'GAPDH'),np.nan)
        ref_target=means.get((reference,gene),np.nan)
        ref_gapdh=means.get((reference,'GAPDH'),np.nan)
        ddct=(mean-gapdh)-(ref_target-ref_gapdh)
        ratio=2.0**(-ddct)
        status='quantified'
        if gene=='GAPDH': status='reference'
        elif valid==0: status='ND'
        elif not np.isfinite(ratio): status='NE'
        elif literal40 or (np.isfinite(mismatch) and abs(mismatch)>.01): status='QC'
        elif valid<n: status='partial'
        rows.append(dict(dataset=dataset,run=run,source_file=name,sample=sample,gene=gene,
                         source_target=group['Target Name'].iloc[0],reference_sample=reference,
                         n_wells=n,n_detected=valid,reported_mean_ct=mean,gapdh_mean_ct=gapdh,
                         reference_target_ct=ref_target,reference_gapdh_ct=ref_gapdh,
                         expression_ratio=ratio,reduction_percent=100*(1-ratio),status=status))
        for item in group.itertuples():
            ratio_well=2.0**(-((item.ct_numeric-gapdh)-(ref_target-ref_gapdh)))
            wells.append(dict(dataset=dataset,run=run,source_file=name,source_row=item.source_row,
                              sample=sample,gene=gene,ct=item.ct_numeric,
                              expression_ratio=ratio_well,summary_status=status))
RUNS=pd.DataFrame(rows)
RUNS.to_csv(DATA/'normalised_run_summaries.csv',index=False)
pd.DataFrame(wells).to_csv(DATA/'normalised_wells.csv',index=False)
pd.DataFrame(qc).to_csv(DATA/'assay_qc_inventory.csv',index=False)
assert np.allclose(RUNS.loc[RUNS['sample'].eq(RUNS.reference_sample),'expression_ratio'].dropna(),1)
summaries=[]
for (dataset,sample,gene),g in RUNS.groupby(['dataset','sample','gene'],sort=False):
    states=set(g.status)
    status=next((x for x in ['reference','QC','ND','NE','partial'] if x in states),'quantified')
    # No average between detected and undetermined runs: require full valid set.
    ratio=g.expression_ratio.mean() if g.expression_ratio.notna().all() else np.nan
    summaries.append(dict(dataset=dataset,sample=sample,gene=gene,
                          n_runs=len(g),n_wells=int(g.n_wells.sum()),n_detected=int(g.n_detected.sum()),
                          expression_ratio=ratio,reduction_percent=100*(1-ratio),status=status))
SUMMARY=pd.DataFrame(summaries)
SUMMARY.to_csv(DATA/'expression_reduction_summary.csv',index=False)
coords=pd.read_csv(OUT/'source_data/ht1080_initial_rtqpcr/y_gene_coordinates_grch38.csv').set_index('gene')
bands=pd.read_csv(OUT/'source_data/ht1080_initial_rtqpcr/chrY_cytobands_grch38.csv')
length=bands.end.max()
colours=LinearSegmentedColormap.from_list('reduction',['#3B6FB6','#FAFAFA','#C75A2B'])
colours.set_bad('#E6E6E6')
shades={'gneg':'#FFFFFF','gpos50':'#B9B9B9','gvar':'#EAEAEA'}
figure_records=[]

def draw(dataset, samples, genes, stem, title, footnote, full=False):
    subset=SUMMARY[SUMMARY.dataset.eq(dataset)].set_index(['sample','gene'])
    values=np.full((len(samples),len(genes)),np.nan)
    labels={}
    for i,(sample,label) in enumerate(samples):
        for j,gene in enumerate(genes):
            row=subset.loc[(sample,gene)]
            v=row.reduction_percent
            status=row.status
            if status in ['reference','ND','NE','QC']:
                text={'reference':'ref.','ND':'ND','NE':'NE','QC':'QC'}[status]
            else:
                values[i,j]=v
                text=f'{v:.0f}%'
                if v>=99.9: text='>99.9%'
                elif abs(v)>=99.5: text=f'{v:.1f}%'
                if status=='partial': text+='*'
            labels[i,j]=text
    finite=values[np.isfinite(values)]
    assert finite.size and finite.min()>=-100 and finite.max()<=100
    fig=plt.figure(figsize=(183/25.4,(124 if len(samples)>3 else 112)/25.4))
    left,width=.25,.57
    m=fig.add_axes([left,.70,width,.20])
    m.set_xlim(0,1); m.set_ylim(0,1); m.set_axis_off()
    for band in bands.itertuples():
        x0,x1=band.start/length,band.end/length
        if band.stain=='acen':
            verts=[(x0,.56),(x0,.72),(x1,.64)] if band.band.startswith('p') else [(x0,.64),(x1,.72),(x1,.56)]
            m.add_patch(Polygon(verts,fc='#777777',ec='#777777',lw=.5))
        else: m.add_patch(Rectangle((x0,.56),x1-x0,.16,fc=shades[band.stain],ec='#777777',lw=.5))
    for j,gene in enumerate(genes):
        position_gene={'UTY2':'UTY'}.get(gene,gene)
        if position_gene not in coords.index: continue
        row=coords.loc[position_gene]
        gx=((row.start-1+row.end)/2)/length
        dest=(j+.5)/len(genes)
        m.plot([gx,gx],[.55,.74],color='#353535',lw=1.2)
        m.plot([gx,gx,dest,dest],[.55,.49,.31,.24],color='#999999',lw=.7)
    m.text(0,.9,'Y chromosome',fontsize=8.5,ha='left')
    m.text(1,.9,'GRCh38 | 57.2 Mb',fontsize=8,ha='right',color='#555555')
    fig.text(.5,.955,title,ha='center',fontsize=11,weight='semibold')
    ax=fig.add_axes([left,.29,width,.39])
    im=ax.imshow(np.ma.masked_invalid(values),cmap=colours,norm=TwoSlopeNorm(vmin=-100,vcenter=0,vmax=100),aspect='auto')
    ax.set_xticks(range(len(genes)),[g if g!='CD99' else 'CD99' for g in genes],fontsize=9 if full else 10,weight='semibold')
    ax.xaxis.tick_top()
    ax.set_yticks(range(len(samples)),[s[1] for s in samples],fontsize=9)
    ax.tick_params(axis='both',length=0,pad=8)
    for spine in ax.spines.values():spine.set_visible(False)
    for x in np.arange(.5,len(genes),1):ax.axvline(x,color='white',lw=2)
    for y in np.arange(.5,len(samples),1):ax.axhline(y,color='white',lw=2)
    for (i,j),label in labels.items():
        value=values[i,j]
        ax.text(j,i,label,ha='center',va='center',fontsize=9 if full else 11,
                color='white' if np.isfinite(value) and abs(value)>65 else '#222222',weight='semibold')
    cax=fig.add_axes([.87,.33,.022,.30])
    cb=fig.colorbar(im,cax=cax,ticks=[-100,0,100])
    cb.outline.set_visible(False)
    cb.ax.tick_params(labelsize=8,length=2)
    cb.set_label('Expression reduction (%)',labelpad=8,fontsize=9)
    fig.text(.5,.225,'Reduction = 100 \u00d7 (1 \u2212 relative expression)',ha='center',fontsize=10)
    fig.text(.5,.165,footnote,ha='center',va='top',fontsize=8.5,color='#555555',linespacing=1.5)
    if dataset=='ht1080_clones':
        extra='ND: all target wells undetermined; no numerical reduction assigned\n* Detected-well estimate only (1/4 wells); not a fully quantified sample'
    elif full and dataset=='ht1080_cp2':
        extra='*: fewer than 3 wells detected; ND: all target wells undetermined\nref.: GAPDH normalisation; CD99 is an X/Y target'
    elif full:
        extra='ref.: GAPDH normalisation; QC: source Ct=40 / mean inconsistency\nNE: ratio not estimable; CD99 is an X/Y target'
    else:extra='Positive: lower expression than reference     Negative: higher expression'
    fig.text(.5,.065,extra,ha='center',va='top',fontsize=8,color='#555555',linespacing=1.5)
    for ext in ['pdf','svg','png']:fig.savefig(OUT/f'{stem}.{ext}',dpi=600)
    fig.savefig(OUT/'previews'/f'{stem}.png',dpi=180)
    plt.close(fig)
    figure_records.append(dict(stem=stem,dataset=dataset,samples=[s[0] for s in samples],genes=genes,
                               title=title,colour_limits=[-100,100],size_mm=[183,124 if len(samples)>3 else 112]))

core=['RPS4Y1','KDM5D','EIF1AY']
hguides=[('H2452 '+g,g) for g in ['CP','Q53','P47','P51']]
draw('h2452',hguides,core,'ch03_main_h2452_rtqpcr_three_genes_ideogram',
     'H2452 | individual-guide testing','Relative to H2452 WT; one assay; three technical wells per target')
draw('h2452',hguides+[('KP4','KP4 (context)'),('A375','A375 (context)')],
     ['RPS4Y1','UTY','KDM5D','EIF1AY','CD99','GAPDH'],
     'ch03_supp_h2452_rtqpcr_all_genes_ideogram','H2452 | all recorded targets',
     'H2452 WT reference; KP4 and A375 are separate cell-line comparators',True)
draw('ht1080_cp2',[('HT1080 CP2','CP2')],core,
     'ch03_main_ht1080_cp2_rtqpcr_three_genes_ideogram','HT1080 | CP2',
     'Relative to run-matched HT1080 WT; mean of two repeat assays\nWorking interpretation: second CP transfection')
draw('ht1080_cp2',[('HT1080 CP2','CP2 (2 runs)'),('A375','A375 (Run 1)'),('KP4','KP4 (Run 2)')],
     ['RPS4Y1','UTY2','KDM5D','EIF1AY','CD99','GAPDH'],
     'ch03_supp_ht1080_cp2_rtqpcr_all_genes_ideogram','HT1080 CP2 | all recorded targets',
     'Run-matched HT1080 WT reference; three technical wells per target\nEIF1AY recorded as E1F1AY; UTY2 assay label retained',True)
draw('ht1080_clones',[('C24','C24'),('KP4','KP4 (context)')],core,
     'ch03_main_ht1080_clones_rtqpcr_three_genes_ideogram','HT1080 clones | C7 Y-positive reference',
     'C24 and KP4 relative to Y-positive C7; four technical wells per target\nKP4 is a separate cell-line comparator')
draw('ht1080_clones',[('C24','C24'),('KP4','KP4 (context)')],core+['GAPDH'],
     'ch03_supp_ht1080_clones_rtqpcr_all_genes_ideogram','HT1080 clones | all recorded targets',
     'Y-positive C7 reference; four technical wells per target\nref.: GAPDH normalising reference',True)
(DATA/'render_provenance.json').write_text(json.dumps(dict(
    figures=figure_records, formula='2^-[(target mean Ct - sample GAPDH mean Ct) - (reference target mean Ct - reference GAPDH mean Ct)]',
    reduction='100*(1-expression_ratio); CP2 uses equal-weight arithmetic mean of normalised run ratios.',
    missing='No Ct=40 imputation; recorded literal 40 values flagged. All-undetermined samples ND; partial detection flagged.',
    clone_reference='Author confirmed C7 is Y-positive on 17 September 2026; C7 is used as reference.',
    cp2='Author response: it must; working interpretation of repeat assays after a second CP transfection.',
    statistics='Descriptive only; no inference from technical wells or controls from different cell lines.'),indent=2),encoding='utf-8')
print('Processed four tables; rendered six figures with matching source tables and QC flags.')
