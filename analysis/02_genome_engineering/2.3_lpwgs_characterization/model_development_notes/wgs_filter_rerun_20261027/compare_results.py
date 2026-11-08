"""Compare historical, mask-only and fully rerun WGS; never update selected assets."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.backends.backend_pdf import PdfPages

OUT = Path(__file__).resolve().parent
THESIS = OUT.parents[2]
FIG = THESIS / 'figures/chapter3/drafts/wgs_filter_rerun_20261027'
FIG.mkdir(parents=True, exist_ok=True)
for font in (THESIS/'skills/thesis-figure-style/assets/fonts').rglob('*.ttf'):
    font_manager.fontManager.addfont(str(font))
plt.rcParams.update({'font.family':'Open Sans','font.size':9,'axes.spines.top':False,
    'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
SAMPLES = {'HT1080':('HT1080-sgCTRL-D0','HT1080-LOY-D0'),
           'A549':('A549-CTRL2','A549-LOY'), 'DU145':('LNCaP-sgCTRL','LNCaP-LOY-2X')}
annotation = pd.read_csv(OUT/'bin_annotation.csv').set_index('bin_id')
regions = [('PAR1',0,2458320),('MSY Yp',2458320,10565750),
           ('centromere',10565750,10883085),('euchromatic Yq',10883085,27449937),
           ('Yq12',27449937,62122809),('PAR2',62122809,62460029)]
def load_pair(cell, mode):
    d = pd.read_csv(OUT/cell/f'{mode}_bins.csv',dtype={'chromosome':str})
    c,l = SAMPLES[cell]
    keys = ['bin_id','chromosome','start','end']
    z = d[d['sample'].eq(c)].merge(d[d['sample'].eq(l)],on=keys,suffixes=('_control','_depleted'),validate='one_to_one')
    z = z.set_index('bin_id').join(annotation[['use']].rename(columns={'use':'original_mask'}))
    z['midpoint'] = (z.start+z.end)/2e6
    z['region']='boundary-crossing'
    ranges = regions.copy()
    if cell=='DU145':
        ranges[1:2]=[('distal MSY Yp',2458320,4800000),('proximal MSY Yp',4800000,10565750)]
    for label, left,right in ranges:
        z.loc[z.chromosome.eq('Y') & z.start.gt(left) & z.end.le(right),'region']=label
    z['valid'] = z.use_control & z.use_depleted & np.isfinite(z[['copynumber_control','copynumber_depleted','segmented_control','segmented_depleted']]).all(axis=1) & z.copynumber_control.gt(.001)
    z['included'] = z.valid & ~z.region.isin(['boundary-crossing','centromere','Yq12']) & z.chromosome.eq('Y')
    z['percent'] = 100*z.copynumber_depleted.clip(lower=.001)/z.copynumber_control.clip(lower=.001)
    z['percent_without_floor'] = 100*z.copynumber_depleted.clip(lower=0)/z.copynumber_control
    z['depleted_at_floor'] = z.copynumber_depleted.le(.001)
    z['delta_cn'] = z.copynumber_depleted-z.copynumber_control
    z['log2_ratio'] = np.log2(z.copynumber_depleted.clip(lower=.001)/z.copynumber_control.clip(lower=.001))
    z['segment_log2_ratio']=np.log2(z.segmented_depleted.clip(lower=.001)/z.segmented_control.clip(lower=.001))
    return z

all_pairs={}; summary=[]; genome=[]; chromosomes=[]; validations=[]
for cell in SAMPLES:
    old=load_pair(cell,'historical'); new=load_pair(cell,'corrected')
    baseline=load_pair(cell,'baseline_reproduction') if (OUT/cell/'baseline_reproduction_bins.csv').exists() else None
    all_pairs[cell]=(old,new)
    # Reproduce the published source median, including the legacy numerical floor.
    published=pd.read_csv(THESIS/f'figures/chapter3/source_data/regional_y_atlases_normalised/{cell.lower()}/regional_median_paired_ratios.csv').set_index('region')
    for r,row in published.iterrows():
        v=old.loc[old.included & old.region.eq(r),'percent']
        assert len(v)==int(row.n_bins)
        assert np.isclose(v.median(),row.median_percent_of_sgCTRL,atol=1e-8)
    for mode,z in [('current',old),('mask_only',old),('full_rerun',new),('baseline_reproduction',baseline)]:
        if z is None: continue
        eligible=z.included & (z.original_mask if mode=='mask_only' else True)
        for r in [r for r in old.region.unique() if r not in ('boundary-crossing','centromere','Yq12')]:
            v=z.loc[eligible & z.region.eq(r)]
            summary.append({'cell':cell,'region':r,'analysis':mode,'n_bins':len(v),
                'median_percent':v.percent.median(),'median_percent_without_floor':v.percent_without_floor.median(),
                'n_at_floor':int(v.depleted_at_floor.sum()),'bin_ids':';'.join(v.index)})
    for mode,z in [('current',old),('full_rerun',new)]:
        z[z.chromosome.eq('Y')].to_csv(OUT/cell/f'{mode}_Y_bin_comparison.csv')
    common=old.valid & new.valid & new.original_mask
    autosomes=common & old.chromosome.isin([str(i) for i in range(1,23)])
    for label,keep in [('autosomes',autosomes),('non-Y',common & old.chromosome.ne('Y'))]:
        x=old.loc[keep,'delta_cn']; y=new.loc[keep,'delta_cn']
        genome.append({'cell':cell,'scope':label,'common_bins':int(keep.sum()),
            'old_new_bin_delta_pearson':x.corr(y),
            'median_abs_change_in_delta_cn':(y-x).abs().median(),
            'p95_abs_change_in_delta_cn':(y-x).abs().quantile(.95),
            'old_median_delta_cn':x.median(),'new_median_delta_cn':y.median(),
            'old_new_segment_log2_ratio_pearson':old.loc[keep,'segment_log2_ratio'].corr(new.loc[keep,'segment_log2_ratio']),
            'median_abs_segment_log2_ratio_change':(new.loc[keep,'segment_log2_ratio']-old.loc[keep,'segment_log2_ratio']).abs().median(),
            'p95_abs_segment_log2_ratio_change':(new.loc[keep,'segment_log2_ratio']-old.loc[keep,'segment_log2_ratio']).abs().quantile(.95),
            'fraction_segment_log2_ratio_changes_gt_0_1':(new.loc[keep,'segment_log2_ratio']-old.loc[keep,'segment_log2_ratio']).abs().gt(.1).mean()})
    for chrom in [str(i) for i in range(1,23)]+['X','Y']:
        keep=common & old.chromosome.eq(chrom)
        chromosomes.append({'cell':cell,'chromosome':chrom,'n_common_bins':int(keep.sum()),
            'old_median_log2_ratio':old.loc[keep,'log2_ratio'].median(),
            'new_median_log2_ratio':new.loc[keep,'log2_ratio'].median(),
            'old_median_segment_log2_ratio':old.loc[keep,'segment_log2_ratio'].median(),
            'new_median_segment_log2_ratio':new.loc[keep,'segment_log2_ratio'].median()})
    # CBS is stochastic and its post-segmentation normalization can rescale bins.
    err=None if baseline is None else float((old[['copynumber_control','copynumber_depleted']]-baseline[['copynumber_control','copynumber_depleted']]).abs().max().max())
    scale_checks=[]
    if baseline is not None:
        for condition in ('control','depleted'):
            col=f'copynumber_{condition}'
            positive=old[col].gt(.05) & baseline[col].gt(.05)
            scale=(baseline.loc[positive,col]/old.loc[positive,col]).median()
            shape_error=(baseline[col]-old[col]*scale).abs().max()
            scale_checks.append({'condition':condition,'baseline_to_historical_scale':float(scale),'max_error_after_scaling':float(shape_error)})
            assert shape_error < 1e-7, (cell,condition,shape_error)
    validations.append({'cell':cell,'published_regional_values_reproduced':True,
        'max_baseline_CN_error':err,'original_mask_Y_bins':int((new.original_mask & new.chromosome.eq('Y')).sum()),
        'eligible_new_Y_bins':int(new.included.sum()),'PAR_bins':int((new.included & new.region.isin(['PAR1','PAR2'])).sum()),
        'baseline_sample_scale_checks':scale_checks})
    assert not (new.included & new.region.isin(['PAR1','PAR2'])).any()

summary=pd.DataFrame(summary);summary.to_csv(OUT/'regional_comparison.csv',index=False)
pd.DataFrame(genome).to_csv(OUT/'genomewide_comparison.csv',index=False)
pd.DataFrame(chromosomes).to_csv(OUT/'chromosome_comparison.csv',index=False)
(OUT/'comparison_validation.json').write_text(json.dumps(validations,indent=2))

# A review figure, not a replacement for any draft-2 figure.
fig,axes=plt.subplots(3,2,figsize=(11.2,8.6),gridspec_kw={'width_ratios':[1.5,1]},layout='constrained')
track_top=max(150,50*np.ceil(max(z.loc[z.included & z.midpoint.lt(27.5),'percent'].max()
    for pair in all_pairs.values() for z in pair)/50))
for i,(cell,(old,new)) in enumerate(all_pairs.items()):
    ax=axes[i,0]
    a=old[old.included & old.midpoint.lt(27.5)]; b=new[new.included & new.midpoint.lt(27.5)]
    ax.scatter(a.midpoint,a.percent,s=9,c='#b0b0b0',alpha=.7,label='Current analysis')
    ax.scatter(b.midpoint,b.percent,s=23,c='#7E3F8F',marker='D',label='Corrected rerun',zorder=4)
    ax.axhline(100,color='#555555',ls='--',lw=.7)
    ax.axvspan(0,2.45832,color='#eeeeee',zorder=0);ax.axvspan(10.56575,10.883085,color='#dddddd',zorder=0)
    ax.text(1.2,track_top*.95,'PAR1\nexcluded',ha='center',va='top',fontsize=8)
    ax.set(xlim=(0,27.5),ylim=(-5,track_top),ylabel=f'{cell}\nDNA (% of sgCTRL)',xlabel='Chromosome Y position (Mb)')
    ax.text(.98,.95,f'{len(b)} usable bins',transform=ax.transAxes,ha='right',va='top',fontsize=8)
    if i==0:ax.legend(loc='upper center',bbox_to_anchor=(.65,.87),frameon=False,fontsize=8)
    ax=axes[i,1]
    dd=summary[(summary.cell==cell)&(~summary.region.isin(['PAR1','PAR2']))]
    rr=dd.region.unique()
    labels=[]
    for j,r in enumerate(rr):
        vals=dd[dd.region.eq(r)].set_index('analysis')
        before=vals.loc['current','median_percent'];after=vals.loc['full_rerun','median_percent']
        ax.plot([before,after],[j,j],color='#999999',lw=1.5)
        ax.scatter(before,j,c='#aaaaaa',s=28)
        ax.scatter(after,j,c='#7E3F8F',marker='D',s=35,zorder=4)
        ax.annotate(f'{before:.2f} to {after:.2f}%',(max(before,after)+2,j),va='center',fontsize=9)
        labels.append(r.replace('euchromatic ','euchr. ')+f"\n(n = {int(vals.loc['full_rerun','n_bins'])})")
    ax.set(yticks=np.arange(len(rr)),yticklabels=labels,xlim=(-2,88),ylim=(len(rr)-.4,-.6),xlabel='Regional median (% of sgCTRL)')
    ax.spines['left'].set_visible(False);ax.tick_params(axis='y',length=0)
fig.savefig(FIG/'ch03_supp_wgs_filter_before_after.png',dpi=300)
fig.savefig(FIG/'ch03_supp_wgs_filter_before_after.svg')
fig.savefig(FIG/'ch03_supp_wgs_filter_before_after.pdf')
plt.close(fig)

# Genome-wide old/new tracks, equal axes within each model; rejected bins omitted in new tracks.
fig,axes=plt.subplots(3,1,figsize=(11.2,8.2),layout='constrained')
offset=0; offsets={};ticks=[];labels=[]
for chrom in [str(i) for i in range(1,23)]+['X','Y']:
    offsets[chrom]=offset
    length=annotation.loc[annotation.chromosome.astype(str).eq(chrom),'end'].max()/1e6
    ticks.append(offset+length/2);labels.append(chrom);offset+=length+8
for ax,(cell,(old,new)) in zip(axes,all_pairs.items()):
    for z,c,label,s in [(old,'#bbbbbb','Current',2),(new,'#7E3F8F','Corrected',3)]:
        d=z[z.valid & z.chromosome.isin(offsets)].copy()
        x=d.midpoint+d.chromosome.map(offsets)
        ax.scatter(x,d.delta_cn,s=s,c=c,alpha=.6,label=label,rasterized=True)
    ax.axhline(0,color='#333333',lw=.7)
    ax.set(xticks=ticks,xticklabels=labels,xlim=(-10,offset),ylabel=f'{cell}\nCN difference\n(depleted − control)')
    limit=max(1.5,max(z.loc[z.valid & z.chromosome.isin(offsets),'delta_cn'].abs().max() for z in (old,new))*1.03)
    ax.set_ylim(-limit,limit)
    ax.legend(loc='upper right',frameon=False,markerscale=3,ncol=2)
axes[-1].set_xlabel('Chromosome')
fig.savefig(FIG/'ch03_supp_wgs_filter_genomewide.png',dpi=300)
fig.savefig(FIG/'ch03_supp_wgs_filter_genomewide.pdf')
plt.close(fig)
print(summary[['cell','region','analysis','n_bins','median_percent']].to_string(index=False))
print(pd.DataFrame(genome).to_string(index=False))
print('Figures:',FIG)
