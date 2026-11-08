"""Read-only diagnostic of the apparent GFP-side split; no figure/gate changes."""
from pathlib import Path
import csv
import json
import os
import sys

THESIS=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(THESIS.parent/'.figure_python_deps'))
sys.path.insert(0,str(THESIS.parent/'comp_assay'))
import numpy as np
from spectral.a8_pipeline.fcs_io import load_fcs
from spectral.a8_pipeline.competition_report import load_cef_gates
from spectral.a8_pipeline.gating import points_in_polygon
OUT=THESIS/'figures/chapter4/review/gfp_split_diagnostic_20261019'
OUT.mkdir(parents=True,exist_ok=True)
os.environ.setdefault('MPLCONFIGDIR',str(OUT/'mpl_cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0,str(THESIS/'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
apply_style(8)
REPO=THESIS.parent/'comp_assay'
EXPORT=REPO/'20261016_171604_Export/comp_growth_assay_D0_OCT_15'
saved=json.loads((REPO/'analysis/spectral_growth_assay_D0_20261015/results.json').read_text())
gates=load_cef_gates(next(EXPORT.glob('*.cef')))
weights=np.asarray(saved['classifier']['weights'])
intercept=saved['classifier']['intercept']
names=['A01-LOY_CHERRY_WT_1','A04-LOY_GFP_WT_1_001',
       'A10-CTRL_GFP_SOLO_WT_003','A08-LOY_GFP_SOLO_WT_001',
       'A07-LOY_CHERRY_SOLO_WT','A09-CTRL_CHERRY_SOLO_WT_002']
stored={}
summary=[]
for name in names:
    data=load_fcs(EXPORT/'Plate_004'/(name+'.fcs'))
    events=data.events
    keep=np.ones(len(events),dtype=bool)
    for gname in ['Scatter','SSC Singlets','FSC Singlets']:
        gate=gates[gname]
        keep &= points_in_polygon(events.loc[:,list(gate.parameters)].to_numpy(float),gate.vertices)
    ev=events.loc[keep]
    raw=ev[['eGFP-A','mCherry-A']].to_numpy(float)
    transformed=np.arcsinh(raw/1000)
    negative=(transformed[:,0]<6)&(transformed[:,1]<6.5)
    gfp=(transformed@weights+intercept>0)&~negative
    mcherry=(transformed@weights+intercept<=0)&~negative
    stored[name]=(raw,transformed,gfp,mcherry,negative)
    for cls,mask in [('GFP_upper',gfp&(raw[:,1]>=0)),('GFP_lower',gfp&(raw[:,1]<0)),('mCherry',mcherry)]:
        if not mask.any(): continue
        selected=ev.loc[mask]
        entry={'sample':name,'group':cls,'events':int(mask.sum()),'total_gfp_events':int(gfp.sum()),
               'gfp_unmixed_median':float(selected['eGFP-A'].median()),
               'mcherry_unmixed_median':float(selected['mCherry-A'].median()),
               'raw_YG4_median':float(selected['YG4 (625)-A'].median()),
               'raw_B3_median':float(selected['B3 (530)-A'].median()),
               'raw_YG4_above_historical_unstained_pct':float(100*(selected['YG4 (625)-A']>saved['unstained_reference']['mcherry_threshold']).mean())}
        summary.append(entry)
    if gfp.any():
        print(json.dumps({'sample':name,'GFP_events':int(gfp.sum()),'positive_mCherry_fraction':float((raw[gfp,1]>0).mean()),'GFP_mCherry_raw_quantiles':np.quantile(raw[gfp,1],[.01,.05,.25,.5,.75,.95,.99]).tolist()}))
with (OUT/'cluster_summary.csv').open('w',newline='',encoding='utf-8') as handle:
    writer=csv.DictWriter(handle,fieldnames=summary[0].keys());writer.writeheader();writer.writerows(summary)

fig,axes=plt.subplots(2,3,figsize=(183/25.4,125/25.4),gridspec_kw={'wspace':.45,'hspace':.65})
chosen=[('A01-LOY_CHERRY_WT_1','WT mixture A01'),('A10-CTRL_GFP_SOLO_WT_003','CTRL GFP-only A10')]
bins=np.linspace(-500000,500000,61)
plot_audit=[]
for row,(name,label) in enumerate(chosen):
    raw,transformed,gfp,mcherry,negative=stored[name]
    for column,cofactor in [(0,1000),(2,100000)]:
        ax=axes[row,column]
        display=np.arcsinh(raw/cofactor)
        for mask,color in [(negative,'#BBBBBB'),(mcherry,'#D34F4F'),(gfp,'#32965B')]:
            ax.scatter(display[mask,0],display[mask,1],s=1.6,c=color,alpha=.3,linewidths=0,rasterized=True)
        ax.axhline(0,color='#555555',lw=.55,ls=':')
        ax.set_xlabel('GFP',fontsize=7);ax.set_ylabel('mCherry',fontsize=7)
        ax.set_xlim((-8,12.5) if cofactor==1000 else (-3,8))
        ax.set_ylim((-8,12.5) if cofactor==1000 else (-5,8))
        ax.set_title(f'{label}\narcsinh cofactor {cofactor:,}',fontsize=7.5,pad=5)
        ax.tick_params(labelsize=6.5,length=2)
    hist=axes[row,1]
    hist.hist(raw[gfp,1],bins=bins,color='#32965B',alpha=.65)
    hist.axvline(0,color='#555555',lw=.7,ls=':')
    hist.set_title('Same GFP-class events\nLinear unmixed mCherry',fontsize=7.5,pad=5)
    hist.set_xlabel('Unmixed mCherry signal',fontsize=7)
    hist.set_ylabel('Events',fontsize=7)
    hist.set_xlim(-500000,500000)
    hist.set_xticks([-400000,0,400000],['−400k','0','400k'])
    hist.tick_params(labelsize=6.5,length=2)
    visible=((raw[gfp,1]>=-500000)&(raw[gfp,1]<=500000)).mean()
    hist.text(.98,.97,f'{100*visible:.1f}% in view',transform=hist.transAxes,
              ha='right',va='top',fontsize=6,color='#666666')
    plot_audit.append({'sample':name,'GFP_total':int(gfp.sum()),'histogram_range':[-500000,500000],
                       'GFP_outside_histogram_range':int(((raw[gfp,1]<-500000)|(raw[gfp,1]>500000)).sum())})
fig.subplots_adjust(left=.075,right=.985,bottom=.15,top=.84)
fig.text(.5,.975,'Display diagnostic: identical events, unchanged classification',ha='center',va='top',fontsize=9,fontweight='semibold')
fig.text(.5,.025,'Cofactor 100,000 is a display sensitivity check, not a newly selected analysis setting.',ha='center',fontsize=7,color='#666666')
fig.savefig(OUT/'gfp_split_display_check.png',dpi=210)
fig.savefig(OUT/'gfp_split_display_check.pdf',dpi=600)
plt.close(fig)
(OUT/'diagnostic_metadata.json').write_text(json.dumps({'samples':names,'histogram_display_audit':plot_audit,'gating_and_classifier':'Original October gates and saved classifier; no fitting or editing.','primary_figure_changed':False},indent=2)+'\n')
print(json.dumps(summary,indent=2))
