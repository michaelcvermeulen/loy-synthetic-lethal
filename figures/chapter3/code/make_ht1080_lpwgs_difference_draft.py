"""Draft F03 alternative; leaves reviewed exports and original measurements intact."""
from pathlib import Path
import sys, io, json, hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

FIG = Path(__file__).resolve().parents[1]
THESIS = FIG.parents[1]
SOURCE = FIG/'source_data/ht1080_lpwgs_confirmation/qdnaseq_bins.csv'
CYTO = THESIS.parent/'WGS/chm13v2.0_cytobands_allchrs.bed'
OUT = FIG/'source_data/ht1080_lpwgs_difference_draft'
OUT.mkdir(exist_ok=True)
sys.path.insert(0, str(THESIS/'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
apply_style(8)
STEM = 'ch03_main_ht1080_lpwgs_difference_draft'
chroms = ['chr'+str(i) for i in range(1, 23)]+['chrX', 'chrY']
cb = pd.read_csv(CYTO, sep='\t', header=None, names=['chrom','start','end','band','stain'])
lengths = cb.groupby('chrom').end.max()
yq = cb[cb.chrom.eq('chrY') & cb.band.eq('q12')].iloc[0]
d = pd.read_csv(SOURCE)
d = d[d.use & np.isfinite(d.log2)]
paired = d.pivot(index=['chromosome','start','end'], columns='label', values='log2').dropna().reset_index()
assert len(paired) == 31185
excluded = paired.chromosome.eq('chrY') & paired.start.le(yq.end) & paired.end.gt(yq.start)
paired.loc[excluded].to_csv(OUT/'Yq12_excluded_source_bins.csv', index=False)
paired = paired.loc[~excluded].copy()
paired['mid'] = (paired.start+paired.end)/2
paired['difference'] = paired.LOY-paired.sgCTRL

def aggregate(data, size):
    data = data.copy()
    data['window_start'] = (data.mid//size).astype(int)*size
    b = data.groupby(['chromosome','window_start']).agg(
        sgCTRL=('sgCTRL','median'), LOY=('LOY','median'), n_bins=('start','size'),
        source_start=('start','min'), source_end=('end','max')).reset_index()
    b['window_end'] = np.minimum(b.window_start+size, b.source_end)
    b['mid'] = (b.window_start+b.window_end)/2
    b['difference'] = b.LOY-b.sgCTRL
    return b

bins = aggregate(paired, 5_000_000)
bins.to_csv(OUT/'genome_5Mb.csv', index=False)
arms = []
for ch in chroms[:-1]:
    cen = cb[cb.chrom.eq(ch) & cb.stain.eq('acen')]
    assert len(cen) == 2
    for arm, left, right in [('p',0,cen.start.min()),('q',cen.end.max(),lengths[ch])]:
        v = paired[paired.chromosome.eq(ch) & paired.mid.gt(left) & paired.mid.le(right)]
        arms.append(dict(chromosome=ch, arm=ch[3:]+arm, n_bins=len(v),
                         difference=v.difference.median(),
                         low_mappability_arm=ch in ['chr13','chr14','chr15','chr21','chr22'] and arm=='p'))
arms = pd.DataFrame(arms)
arms.to_csv(OUT/'arm_medians.csv', index=False)
# Retain the original close-up scope; arm boundaries and the old segment boundary
# are annotations only. No new segment fitting or absolute copy-number calls.
detail = aggregate(paired[paired.chromosome.eq('chr13') & paired.start.ge(18_100_001)], 1_000_000)
detail.to_csv(OUT/'chr13q_1Mb.csv', index=False)
segment_delta = 0.479599094060819-0.111904199476798
offsets = {}; end = 0
for ch in chroms:
    offsets[ch] = end
    end += lengths[ch]/1e6+20

fig = plt.figure(figsize=(8.4,6.6))
gs = fig.add_gridspec(3,1, height_ratios=[1.15,1,1], left=.11,right=.98,top=.97,bottom=.085,hspace=.52)
ax = fig.add_subplot(gs[0])
for i,ch in enumerate(chroms):
    z = bins[bins.chromosome.eq(ch)]
    ax.scatter(z.mid/1e6+offsets[ch],z.difference,s=8,
               c='#E64B35' if ch=='chrY' else ['#555555','#AAAAAA'][i%2],lw=0)
ax.axhline(0,c='#333333',lw=.5,zorder=0)
ax.axvspan(offsets['chrY']+yq.start/1e6,offsets['chrY']+yq.end/1e6,color='#EEEEEE',zorder=-1)
ax.annotate('Yq12 excluded',xy=(offsets['chrY']+(yq.start+yq.end)/2e6,0),xytext=(.985,.86),
            textcoords='axes fraction',ha='right',fontsize=7,
            arrowprops=dict(arrowstyle='-',lw=.5,color='#888888'),color='#666666')
ax.set_yscale('symlog',linthresh=.3)
extent=max(10,bins.difference.abs().max()*1.12)
ax.set(ylim=(-extent,extent),xlim=(-10,end),ylabel='$\u0394$ log$_2$(relative CN)')
ax.set_yticks([-10,-1,-.3,0,.3,1,10],['−10','−1','−0.3','0','0.3','1','10'])
ax.tick_params(axis='y',labelsize=7)
ax.set_xticks([offsets[ch]+lengths[ch]/2e6 for ch in chroms],[ch[3:] for ch in chroms],fontsize=7)
ax.set_xlabel('Chromosome',labelpad=3)
ax.text(.015,.96,'5-Mb windows',transform=ax.transAxes,va='top',fontsize=7.5)
ax.text(.015,.06,'Y-depleted − sgCTRL',transform=ax.transAxes,va='bottom',fontsize=7.5)

ax = fig.add_subplot(gs[1]); x=np.arange(len(arms))
bars=ax.bar(x,arms.difference,width=.78,color=['#00A087' if a=='13q' else '#666666' for a in arms.arm])
for bar,low in zip(bars,arms.low_mappability_arm):
    if low: bar.set_facecolor('#DDDDDD');bar.set_edgecolor('#999999');bar.set_hatch('///')
for i in range(2,len(arms),2):ax.axvline(i-.5,c='#EEEEEE',lw=.5,zorder=0)
ax.axhline(0,c='#333333',lw=.5)
ax.set(xlim=(-.8,len(arms)-.2),ylim=(-.41,.43),ylabel='Median $\u0394$ log$_2$(relative CN)')
ax.set_xticks(x,arms.arm,rotation=90,fontsize=6.8)
for tick,low in zip(ax.get_xticklabels(),arms.low_mappability_arm):
    if low: tick.set_color('#AAAAAA')
ax.set_xlabel('Chromosome arm',labelpad=4)
ix=arms.index[arms.arm.eq('13q')][0]
ax.text(ix,arms.loc[ix,'difference']+.025,'13q',ha='center',color='#007C68',fontsize=8)
ax.text(.015,.95,'Autosomes and X',transform=ax.transAxes,va='top',fontsize=7.5)

ax = fig.add_subplot(gs[2])
ax.axhline(0,c='#888888',lw=.5)
ax.scatter(detail.mid/1e6,detail.difference,s=13,color='#555555',lw=0,zorder=3,label='1-Mb windows')
ax.axvline(34.7,c='#AAAAAA',lw=.6,ls=':')
ax.plot([34.7,113.566686],[segment_delta]*2,color='#00A087',lw=1.5,label='Original segment difference')
ax.text(77,.512,f'$\u0394$ log$_2$ = {segment_delta:.3f}  (1.29-fold)',
        color='#007C68',fontsize=8,ha='center')
ax.set(xlim=(18.1,113.567),ylim=(min(-.04,detail.difference.min()-.05),max(.55,detail.difference.max()+.07)),
       xlabel='Chromosome 13q position (Mb)',ylabel='$\u0394$ log$_2$(relative CN)')
ax.legend(loc='lower right',frameon=False,fontsize=7,ncol=2,handlelength=1.5,columnspacing=1)
for a,letter in zip(fig.axes,'abc'):
    a.text(-.105,1.035,letter,transform=a.transAxes,fontsize=12,weight='semibold')
    a.spines[['top','right']].set_visible(False)
for ext in ['pdf','svg','png']:
    buf=io.BytesIO();fig.savefig(buf,format=ext,dpi=300,facecolor='white')
    (FIG/f'{STEM}.{ext}').write_bytes(buf.getvalue())
plt.close(fig)
assert len(bins)==630 and excluded.sum()==351
assert len(arms)==46 and np.isfinite(arms.difference).all()
assert arms.difference.between(-.41,.43).all()
assert np.isclose(2**segment_delta,1.29029,rtol=1e-4)
report=dict(input_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [SOURCE,CYTO]},
  comparison='HT1080-LOY-D0 (CP + Q53 2x) minus HT1080-sgCTRL-D0; one WGS specimen per condition',
  source_bins_before_mask=31185,excluded_Yq12_bins=int(excluded.sum()),genome_windows=len(bins),
  Yq12_BED=[int(yq.start),int(yq.end)],coordinate_rule='exclude every overlapping one-based inclusive source bin',
  window_estimator='median corrected unsegmented log2 per sample, then subtract sgCTRL from LOY',
  arm_estimator='median paired 100-kb log2 difference; midpoints outside centromere; all 46 autosomal/X arms',
  arm_caution='Acrocentric short arms are retained, but their repetitive sequence and floor-level values limit interpretation.',
  genome_axis='symlog; linear within +/-0.3; all window values retained',
  chr13_original_segment_delta=segment_delta,chr13_relative_ratio=2**segment_delta,
  font='Open Sans',source_measurements_unchanged=True,selected_figure_unchanged=True)
(OUT/'validation.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
