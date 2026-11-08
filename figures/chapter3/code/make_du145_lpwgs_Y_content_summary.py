"""Summarise whole-Y depth and regional corrected signal for eight DU145 samples."""
from pathlib import Path
import sys, json, hashlib, time
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

ROOT=Path(__file__).resolve().parents[3]
FIG=ROOT/'figures/chapter3'; SRC=FIG/'source_data/du145_drug_optimisation'
sys.path.insert(0,str(ROOT/'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
apply_style(8)
ORDER=['LNCaP-sgCTRL','LNCaP-DMSO-LOY-1X','LNCaP-30','LNCaP-70','LNCaP-KU','LNCaP-KU-30','LNCaP-KU-70','LNCaP-LOY-2X']
LABELS=['sgCTRL','1× vehicle','1× pano 30 nM','1× pano 75 nM','1× KU','1× KU + pano 30 nM','1× KU + pano 75 nM','2×, no drug']
COLORS=json.loads((SRC/'figure_render_provenance.json').read_text())['color_mapping']
depth=pd.read_csv(SRC/'lpwgs/raw_depth_ratios.csv').set_index('sample').loc[ORDER]
bins=pd.read_csv(SRC/'lpwgs/chrY_qdnaseq_bins.csv')
bins=bins[bins.use.eq(True)&bins.log2_seg.notna()].copy()
bins['relative_cn']=2**bins.log2
# Shared atlas boundaries; centromere excluded from the five region summaries.
cyto=pd.read_csv(FIG/'source_data/ht1080_y_expression_confirmation/T2T_chrY_cytobands.csv')
cyto=cyto[cyto.chrom.eq('chrY')]
cen=cyto[cyto.stain.eq('acen')]
p1=2.458320;p2=62.122809
c0=cen.start.min()/1e6;c1=cen.end.max()/1e6
q0=cyto.loc[cyto.band.eq('q12'),'start'].iloc[0]/1e6
length=cyto.end.max()/1e6
REGIONS=[('PAR1',0,p1),('MSY Yp',p1,c0),('MSY Yq',c1,q0),('Yq12',q0,p2),('PAR2',p2,length)]
rows=[]
for s in ORDER:
    for name,a,b in REGIONS:
        v=bins.loc[bins['sample'].eq(s)&bins.midpoint.ge(a)&bins.midpoint.lt(b),'relative_cn']
        rows.append(dict(sample=s,region=name,start_Mb=a,end_Mb=b,n_bins=len(v),median_relative_cn=v.median()))
regional=pd.DataFrame(rows)
# Confirm identical statistics for the two samples shared with the main atlas.
atlas=pd.read_csv(FIG/'source_data/du145_lpwgs_confirmation/regional_medians_verified.csv')
region_alias={'MSY Yq':'MSY Yq (euchromatic)','Yq12':'Yq12 heterochromatin'}
for sample,condition in [('LNCaP-sgCTRL','sgCTRL'),('LNCaP-LOY-2X','LOY')]:
 for row in regional[regional['sample'].eq(sample)].itertuples():
  ref=atlas[atlas.condition.eq(condition)&atlas.region.eq(region_alias.get(row.region,row.region))]
  assert len(ref)==1 and row.n_bins==ref.n_bins.iloc[0] and np.isclose(row.median_relative_cn,ref.median_relative_cn.iloc[0])
regional.to_csv(SRC/'lpwgs/Y_content_summary_regions.csv',index=False)
values=regional.pivot(index='sample',columns='region',values='median_relative_cn').loc[ORDER,[r[0] for r in REGIONS]].values
assert np.isfinite(values).all() and (values>=.001-1e-10).all()
fig=plt.figure(figsize=(183/25.4,126/25.4))
fig.text(.025,.954,'a',fontsize=12,fontweight='semibold',va='top')
fig.text(.29,.94,'Whole-Y coverage',fontsize=9,fontweight='semibold')
fig.text(.565,.954,'b',fontsize=12,fontweight='semibold',va='top')
fig.text(.61,.94,'Regional Y content',fontsize=9,fontweight='semibold')
ax=fig.add_axes([.29,.21,.24,.65]); y=np.arange(8)
for i,s in enumerate(ORDER):
    v=depth.loc[s,'chrY_autosome_ratio']
    ax.hlines(i,0,v,color=COLORS[s],lw=1,alpha=.8)
    ax.scatter(v,i,s=29,color=COLORS[s],zorder=3)
    ax.text(v+.008,i,f'{v:.4f}',va='center',fontsize=7.4,color='#333333')
ax.set(xlim=(0,.345),ylim=(7.6,-.6),yticks=y,yticklabels=LABELS,xticks=[0,.1,.2,.3],xlabel='Y/autosome mean depth')
ax.tick_params(axis='y',length=0,pad=8);ax.spines['left'].set_visible(False)
ax.axhline(.5,color='#DDDDDD',lw=.6);ax.axhline(6.5,color='#DDDDDD',lw=.6)
ax.grid(axis='x',color='#EEEEEE',lw=.5);ax.set_axisbelow(True)
for i,t in enumerate(ax.get_yticklabels()):t.set_color(COLORS[ORDER[i]])
heat=fig.add_axes([.585,.21,.315,.65]);norm=LogNorm(vmin=.001,vmax=1)
im=heat.imshow(values,cmap='magma',norm=norm,aspect='auto',interpolation='nearest')
heat.set(ylim=(7.6,-.6),yticks=[],xticks=range(len(REGIONS)),xticklabels=['PAR1','MSY\nYp','MSY\nYq','Yq12','PAR2'])
heat.tick_params(axis='x',length=0,pad=6)
heat.spines[['top','right','left','bottom']].set_visible(False)
for i in range(8):
    for j in range(len(REGIONS)):
        v=values[i,j]
        text='Floor' if v<=.0010000001 else f'{v:.3f}'
        heat.text(j,i,text,ha='center',va='center',fontsize=7.3,color='#222222' if norm(v)>.62 else 'white')
for yy in np.arange(.5,7.5):heat.axhline(yy,color='white',lw=.7)
for xx in np.arange(.5,len(REGIONS)-.5):heat.axvline(xx,color='white',lw=.7)
cbax=fig.add_axes([.92,.30,.012,.48]);cb=fig.colorbar(im,cax=cbax,ticks=[.001,.01,.1,1]);cb.ax.set_yticklabels(['0.001','0.01','0.1','1'])
cb.ax.tick_params(labelsize=7);cb.set_label('Median corrected relative CN',fontsize=7,labelpad=2)
fig.text(.29,.109,'1× / 2×: one / two nucleofections; pano: panobinostat; KU: KU-55933.',fontsize=7)
fig.text(.29,.075,'One WGS sample per condition. Floor: upstream numerical floor (0.001).',fontsize=7)
fig.text(.04,.025,'Source: DU145 treatment lpWGS; mosdepth and saved QDNAseq outputs.',fontsize=6.7,color='#666666')
stem='ch03_main_du145_lpwgs_Y_content_summary'
for dest,dpi in [(FIG/f'{stem}.{ext}',600) for ext in ['pdf','svg','png']]+[(FIG/'previews'/f'{stem}.png',160)]:
    tmp=dest.with_name(dest.stem+'.staging'+dest.suffix)
    fig.savefig(tmp,dpi=dpi,facecolor='white')
    for attempt in range(10):
        try:tmp.replace(dest);break
        except OSError:
            if attempt==9:raise
            time.sleep(.5)
plt.close(fig)
depth.reset_index().to_csv(SRC/'lpwgs/Y_content_summary_depth.csv',index=False)
provenance={'figure':stem,'size_mm':[183,126],'font':'Open Sans','png_dpi':600,
    'whole_Y_metric':'chrY mean depth / length-weighted chr1-22 mean depth, computed from bases/length',
    'regional_metric':'Median 2**saved log2 in retained midpoint-selected bins; floor 0.001 retained',
    'regions':REGIONS,'atlas_regions_match':True,'shared_sample_statistics_verified':True,'centromere_omitted':True,'samples':ORDER,'heatmap_scale':'magma; log norm 0.001 to 1','biological_n_per_condition':1,
    'tests':None,'outputs':{f.name:hashlib.sha256(f.read_bytes()).hexdigest() for ext in ['pdf','svg','png'] for f in [FIG/f'{stem}.{ext}']}}
(SRC/'Y_content_summary_provenance.json').write_text(json.dumps(provenance,indent=2),encoding='utf-8')
print('Rendered',stem)
