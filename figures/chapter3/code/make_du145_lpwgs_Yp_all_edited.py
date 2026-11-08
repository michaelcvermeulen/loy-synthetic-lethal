"""DU145 PAR1/Yp zoom from saved QDNAseq bins and segments; no new segmentation."""
from pathlib import Path
import sys,json,hashlib,time
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle,Polygon
from matplotlib.lines import Line2D
ROOT=Path(__file__).resolve().parents[3]; FIG=ROOT/'figures/chapter3'; PROJECTS=ROOT.parent
sys.path.insert(0,str(ROOT/'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
apply_style(8)
CELL='DU145'
assert CELL in ['A549','DU145']
SLUG=CELL.lower()
W=FIG/f'source_data/{SLUG}_lpwgs_confirmation'; R=FIG/f'source_data/{SLUG}_y_expression_confirmation'
SHARED=FIG/'source_data/ht1080_y_expression_confirmation'
COMBO='CP + Q53' if CELL=='A549' else 'CP + Q53 + P47'
COL={'sgCTRL':'#3C5488','LOY':'#E64B35'}; LABEL={'sgCTRL':'sgCTRL','LOY':'Y-depleted'}
WGS_LABEL={'sgCTRL':f'{CELL} | sgCTRL','LOY':f'{CELL} | Y-depleted ({COMBO} 2\u00d7)'}
cyto=pd.read_csv(SHARED/'T2T_chrY_cytobands.csv');cyto=cyto[cyto.chrom.eq('chrY')]
length=cyto.end.max()/1e6; cen=cyto[cyto.stain.eq('acen')];c0,c1=cen.start.min()/1e6,cen.end.max()/1e6
p1=2.458320;p2=62.122809;q0=float(cyto[cyto.band.eq('q12')].start.iloc[0])/1e6
regions=[('PAR1',0,p1),('MSY Yp',p1,c0),('centromere',c0,c1),('MSY Yq (euchromatic)',c1,q0),('Yq12 heterochromatin',q0,p2),('PAR2',p2,length)]
exported=[]
def save(fig,stem):
 for ext in ['pdf','svg','png']:
  tmp=FIG/(stem+'.staging.'+ext);fig.savefig(tmp,dpi=600,facecolor='white')
  for attempt in range(10):
   try:tmp.replace(FIG/(stem+'.'+ext));break
   except PermissionError:
    if attempt==9:raise
    time.sleep(.5)
 fig.savefig(FIG/'previews'/(stem+'.png'),dpi=160,facecolor='white');plt.close(fig);exported.append(stem)
def tag(ax,letter):ax.text(-.09,1.10,letter,transform=ax.transAxes,fontsize=12,fontweight='semibold')
def ideogram(ax):
 for b in cyto.itertuples():
  x0,x1=b.start/1e6,b.end/1e6
  if b.stain=='acen':
   pts=[(x0,.2),(x0,.8),(x1,.5)] if b.band.startswith('p') else [(x0,.5),(x1,.2),(x1,.8)]
   ax.add_patch(Polygon(pts,facecolor='#777777',edgecolor='#555555',lw=.4))
  else:ax.add_patch(Rectangle((x0,.2),x1-x0,.6,facecolor={'gneg':'white','gpos50':'#999999','gvar':'#AFC0CD'}.get(b.stain,'#BBBBBB'),edgecolor='#555555',lw=.4))
 for a,b in [(0,p1),(p2,length)]:ax.add_patch(Rectangle((a,.2),b-a,.6,facecolor='#41B6C4',edgecolor='#555555',lw=.4))
 ax.set(xlim=(0,length),ylim=(0,1));ax.axis('off')

d=pd.read_csv(FIG/'source_data/du145_drug_optimisation/lpwgs/plotted_chrY_bins.csv');d=d[d.use & d.log2_seg.notna() & d.chromosome.eq('chrY')].copy()
W=FIG/'source_data/du145_all_edited_Yp_overlay';W.mkdir(exist_ok=True)
d['condition']=d['sample'];d['relative_cn']=2**d.log2;d['segment_relative_cn']=2**d.log2_seg
xmax=c1;d=d[d.midpoint.lt(xmax)].copy()
segments=[]
for condition,g in d.groupby('condition'):
 g=g.sort_values('start');g['run']=g.log2_seg.ne(g.log2_seg.shift()).cumsum()
 for _,h in g.groupby('run'):
  segments.append(dict(condition=condition,start_bp=int(h.start.min()),end_bp=int(h.end.max()),n_bins=len(h),segment_relative_cn=float(h.segment_relative_cn.iloc[0]),median_bin_relative_cn=float(h.relative_cn.median())))
seg=pd.DataFrame(segments);seg.to_csv(W/'PAR1_Yp_saved_segments.csv',index=False);d.to_csv(W/'PAR1_Yp_bins.csv',index=False)
samples=['LNCaP-DMSO-LOY-1X','LNCaP-30','LNCaP-70','LNCaP-LOY-2X','LNCaP-KU','LNCaP-KU-30','LNCaP-KU-70']
labels=['1× DMSO','1× pano 30 nM','1× pano 75 nM','2×, no drug','1× KU (low depth)','1× KU + pano 30 nM','1× KU + pano 75 nM']
colors=['#777777','#E99B57','#E64B35','#705A8C','#4DBBD5','#00A087','#8491B4']
fig=plt.figure(figsize=(183/25.4,160/25.4));gs=fig.add_gridspec(3,1,height_ratios=[1,1.5,.66],left=.10,right=.965,top=.94,bottom=.12,hspace=.55)
for i,shown in enumerate([['LNCaP-sgCTRL'],samples]):
 ax=fig.add_subplot(gs[i]);ax.axvspan(0,p1,color='#EAF6F7');ax.axvspan(c0,c1,color='#EAEAEA');ax.axvline(p1,color='#007786',ls=':',lw=.7)
 for sample in shown:
  color='#3C5488' if i==0 else colors[samples.index(sample)]
  g=d[d['sample'].eq(sample)]
  ax.scatter(g.midpoint,g.relative_cn,color=color,s=6,alpha=.50 if i==0 else .25,lw=0,rasterized=False)
  for row in seg[seg.condition.eq(sample)].itertuples():
   ax.plot([(row.start_bp-1)/1e6,min(row.end_bp/1e6,xmax)],[row.segment_relative_cn]*2,color=color,lw=1.2,ls='--' if sample=='LNCaP-KU' else '-')
 ax.set(xlim=(0,xmax),ylim=(-.025,1.8),yticks=[0,.5,1,1.5],ylabel='Relative CN')
 ax.tick_params(axis='x',labelbottom=False);tag(ax,chr(97+i))
 ax.set_title('DU145 | sgCTRL' if i==0 else 'DU145 | all edited samples',loc='left',fontsize=9)
 if i==0:
  ax.text(p1/2,1.55,'PAR1',ha='center',color='#007786',fontsize=8);ax.text((p1+c0)/2,1.55,'Male-specific Yp',ha='center',fontsize=8)
 else:
  ax.legend(handles=[Line2D([],[],color=c,lw=1.4,ls='--' if 'low depth' in l else '-',label=l) for c,l in zip(colors,labels)],loc='upper right',frameon=False,fontsize=6,ncol=2,columnspacing=.9,handlelength=1.8)
ax=fig.add_subplot(gs[2]);ideogram(ax);ax.set_axis_on();ax.set(xlim=(0,xmax),ylim=(-1.15,2.7),yticks=[],xlabel='HG002 T2T-Y position (Mb)')
# Draw the cytobands without PAR colour overlay; show PAR as its own aligned track.
# White/grey cytobands overlay the generic ideogram to keep band and PAR boundaries distinct.
for row in cyto[cyto.start.lt(xmax*1e6)].itertuples():
 l,r=row.start/1e6,min(row.end/1e6,xmax)
 if row.stain!='acen':ax.add_patch(Rectangle((l,.2),r-l,.6,facecolor='white' if row.stain=='gneg' else '#999999',edgecolor='#555555',lw=.5))
 label_x={'p11.32':.22,'p11.31':1.35,'p11.2':5.5,'p11.1':10.45,'q11.1':10.84}.get(row.band,(l+r)/2)
 if row.band in ['p11.32','p11.31']:
  ax.plot([(l+r)/2,label_x],[.82,1.4],color='#888888',lw=.5)
  ax.text(label_x,1.45,row.band,fontsize=6.7,ha='center')
 elif row.band=='p11.2':ax.text(label_x,.5,row.band,fontsize=8,ha='center',va='center')
ax.add_patch(Rectangle((0,2.05),p1,.25,facecolor='#41B6C4',edgecolor='none'))
ax.text(p1/2,2.4,'PAR1',ha='center',fontsize=7,color='#007786')
ax.text((p1+c0)/2,2.4,'MSY Yp',ha='center',fontsize=7)
ax.plot([p1,c0],[2.175,2.175],color='#555555',lw=.8)
ax.text(c0-.1,1.5,'Centromere',ha='right',fontsize=6.7,color='#666666')
ax.plot([(c0+c1)/2,c0-.1],[.85,1.4],color='#888888',lw=.5)
mapping=pd.read_csv(R/'gene_mapping_plotted.csv');probes=pd.read_csv(FIG/'source_data/ht1080_lpwgs_confirmation/ddpcr_target_loci.csv')
gene_map=pd.concat([mapping,probes],ignore_index=True).drop_duplicates('gene_name');gene_map=gene_map[gene_map.t2t_start.lt(c0*1e6)].sort_values('t2t_start')
last=-1
for row in gene_map.itertuples():
 pos=(row.t2t_start+row.t2t_end)/2e6;label_x=max(pos,last+.8);last=label_x
 ax.plot([pos,pos,label_x],[.18,-.10,-.33],color='#888888',lw=.5)
 ax.text(label_x,-.4,row.gene_name,ha='center',va='top',fontsize=7,fontstyle='italic',color='#007786' if row.region in ['PAR1','PAR2'] else '#222222')
ax.spines[['left','top','right']].set_visible(False);ax.set_xticks(np.arange(0,11,1));tag(ax,'c')
save(fig,'ch03_supp_du145_lpwgs_PAR1_Yp_all_edited_draft')
(W/'README.md').write_text("Draft overlay of CH3-S3.4.5.3-S01. Source: ../du145_drug_optimisation/lpwgs/plotted_chrY_bins.csv. All seven edited samples plus sgCTRL; corrected 100-kb points and original saved segments, without new segmentation. Colours follow existing treatment figures. KU-only is marked low depth; one sequenced sample per condition. Raw LNCaP sample IDs denote DU145 as established in project provenance. Raw 70 labels are displayed as author-confirmed 75 nM. Y limits expanded to 1.8 to include all bins. The original selected S01 is unchanged.",encoding='utf-8')
print('Saved all-edited draft.')
