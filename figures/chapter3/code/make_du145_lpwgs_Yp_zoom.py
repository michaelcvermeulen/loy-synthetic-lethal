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

d=pd.read_csv(W/'qdnaseq_bins.csv');d=d[d.use & d.log2_seg.notna() & d.chromosome.eq('chrY')].copy()
d['condition']=d.label.replace({'LOY-2X':'LOY'});d['relative_cn']=2**d.log2;d['segment_relative_cn']=2**d.log2_seg
xmax=c1;d=d[d.midpoint.lt(xmax)].copy()
segments=[]
for condition,g in d.groupby('condition'):
 g=g.sort_values('start');g['run']=g.log2_seg.ne(g.log2_seg.shift()).cumsum()
 for _,h in g.groupby('run'):
  segments.append(dict(condition=condition,start_bp=int(h.start.min()),end_bp=int(h.end.max()),n_bins=len(h),segment_relative_cn=float(h.segment_relative_cn.iloc[0]),median_bin_relative_cn=float(h.relative_cn.median())))
seg=pd.DataFrame(segments);seg.to_csv(W/'PAR1_Yp_saved_segments.csv',index=False);d.to_csv(W/'PAR1_Yp_bins.csv',index=False)
fig=plt.figure(figsize=(183/25.4,145/25.4));gs=fig.add_gridspec(3,1,height_ratios=[1,1,.66],left=.10,right=.965,top=.92,bottom=.13,hspace=.44)
for i,condition in enumerate(COL):
 ax=fig.add_subplot(gs[i]);g=d[d.condition.eq(condition)]
 ax.axvspan(0,p1,color='#EAF6F7');ax.axvspan(c0,c1,color='#EAEAEA')
 ax.axvline(p1,color='#007786',ls=':',lw=.7)
 ax.scatter(g.midpoint,g.relative_cn,color='#888888',s=8,alpha=.7,lw=0)
 for row in seg[seg.condition.eq(condition)].itertuples():
  ax.plot([(row.start_bp-1)/1e6,min(row.end_bp/1e6,xmax)],[row.segment_relative_cn]*2,color=COL[condition],lw=1.7)
 ax.set(xlim=(0,xmax),ylim=(-.025,1.5),yticks=[0,.5,1,1.5],ylabel='Relative CN')
 ax.tick_params(axis='x',labelbottom=False);tag(ax,chr(97+i))
 ax.set_title(WGS_LABEL[condition],loc='left',fontsize=9,color=COL[condition])
 if i==0:
  ax.text(p1/2,1.32,'PAR1',ha='center',color='#007786',fontsize=8)
  ax.text((p1+c0)/2,1.32,'Male-specific Yp',ha='center',fontsize=8)
  ax.legend(handles=[Line2D([],[],marker='o',linestyle='none',color='#888888',markersize=3,label='100-kb bins'),Line2D([],[],color=COL[condition],lw=2,label='Saved segments')],loc='upper right',frameon=False,fontsize=6.5,bbox_to_anchor=(1,1.30),ncol=2)
 else:
  for left,right,level,text in [(0,4.8,1.14,'0-4.8 Mb'),(8.3,8.9,.83,'8.3-8.9 Mb')]:
   ax.plot([left,left,right,right],[level-.04,level,level,level-.04],color='#555555',lw=.6)
   ax.text((left+right)/2,level+.07,text,ha='center',fontsize=7)
  for boundary in [4.8,6,8.3,8.9]:ax.axvline(boundary,color='#BBBBBB',ls=':',lw=.5,zorder=0)
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
save(fig,'ch03_supp_du145_lpwgs_PAR1_Yp_zoom')
(W/'PAR1_Yp_zoom_validation.json').write_text(json.dumps({'source':'qdnaseq_bins.csv; original saved log2_seg values','new_segmentation':False,'display_range_Mb':[0,xmax],'samples':sorted(d['sample'].unique()),'max_plotted_relative_cn':float(d.relative_cn.max()),'retention_interpretation':'Stronger residual signal in source segments 0-4.8 and 8.3-8.9 Mb; not whole-cytoband retention or exact structural breakpoints','cytoband_p11_2_Mb':[.454123,10.565750],'PAR1_end_Mb':p1,'numerical_floor':.001,'files':{f.name:hashlib.sha256(f.read_bytes()).hexdigest() for ext in ['pdf','svg','png'] for f in [FIG/f'ch03_supp_du145_lpwgs_PAR1_Yp_zoom.{ext}']}},indent=2)+'\n')
print(seg.to_string(index=False))
