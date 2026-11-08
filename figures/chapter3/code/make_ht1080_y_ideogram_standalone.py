"""Standalone export of HT1080 regional Y atlas panel d; original plotting geometry and labels, without panel letter."""
"""Restyle source-backed HT1080 genome, regional Y and expression figures."""
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
W=FIG/'source_data/ht1080_lpwgs_confirmation'; R=FIG/'source_data/ht1080_y_expression_confirmation'
COL={'sgCTRL':'#3C5488','LOY':'#E64B35'}; LABEL={'sgCTRL':'sgCTRL','LOY':'Y-depleted'}
WGS_LABEL={'sgCTRL':'sgCTRL','LOY':'Y-depleted (CP + Q53 2\u00d7)'}
cyto=pd.read_csv(R/'T2T_chrY_cytobands.csv');cyto=cyto[cyto.chrom.eq('chrY')]
length=cyto.end.max()/1e6; cen=cyto[cyto.stain.eq('acen')];c0,c1=cen.start.min()/1e6,cen.end.max()/1e6
p1=2.458320;p2=62.122809;q0=float(cyto[cyto.band.eq('q12')].start.iloc[0])/1e6
regions=[('PAR1',0,p1),('MSY Yp',p1,c0),('centromere',c0,c1),('MSY Yq (euchromatic)',c1,q0),('Yq12 heterochromatin',q0,p2),('PAR2',p2,length)]
def ideogram(ax):
 for b in cyto.itertuples():
  x0,x1=b.start/1e6,b.end/1e6
  if b.stain=='acen':
   pts=[(x0,.2),(x0,.8),(x1,.5)] if b.band.startswith('p') else [(x0,.5),(x1,.2),(x1,.8)]
   ax.add_patch(Polygon(pts,facecolor='#777777',edgecolor='#555555',lw=.4))
  else:ax.add_patch(Rectangle((x0,.2),x1-x0,.6,facecolor={'gneg':'white','gpos50':'#999999','gvar':'#AFC0CD'}.get(b.stain,'#BBBBBB'),edgecolor='#555555',lw=.4))
 for a,b in [(0,p1),(p2,length)]:ax.add_patch(Rectangle((a,.2),b-a,.6,facecolor='#41B6C4',edgecolor='#555555',lw=.4))
 ax.set(xlim=(0,length),ylim=(0,1));ax.axis('off')
# Match panel d dimensions in the original five-row atlas.
layout=plt.figure(figsize=(183/25.4,205/25.4))
gs=layout.add_gridspec(5,1,height_ratios=[1,1,.60,.72,1.22],left=.105,right=.90,top=.96,bottom=.10,hspace=.49)
position=gs[3].get_position(layout)
width,height=layout.get_size_inches()*[position.width,position.height]
plt.close(layout)
fig=plt.figure(figsize=(width,height))
ax=fig.add_axes([0,0,1,1])
ideogram(ax)
ax.set_ylim(-1.0,2.55)
ax.set_axis_on()
ax.set_yticks([])
mapping=pd.read_csv(R/'gene_mapping_plotted.csv').sort_values('t2t_start').reset_index(drop=True)
probes=pd.read_csv(W/'ddpcr_target_loci.csv')
for row in probes.itertuples():
 left,right=row.t2t_start/1e6,row.t2t_end/1e6;mid=(left+right)/2
 ax.add_patch(Rectangle((left,.2),right-left,.6,facecolor='#222222',edgecolor='#222222',lw=.5))
 label_x={'RPS4Y1':2.8,'AMELY':7.6}.get(row.gene_name,mid)
 ax.plot([mid,mid,label_x],[.18,-.10,-.25],color='#222222',lw=.7)
 ax.text(label_x,-.30,row.gene_name,ha='center',va='top',fontsize=7.1,fontweight='semibold')
# Light expression-gene labels above the chromosome; target labels remain below.
expression_labels=mapping.copy()
positions=((expression_labels.t2t_start+expression_labels.t2t_end)/2e6).to_numpy()
label_positions=positions.copy()
# Spread crowded Yp/Yq labels while maintaining genomic order and explicit leaders.
for j in range(1,len(label_positions)):
 label_positions[j]=max(label_positions[j],label_positions[j-1]+2.35)
label_positions=np.minimum(label_positions,length-1.0)
for row,label_x,actual_x in zip(expression_labels.itertuples(),label_positions,positions):
 ax.plot([actual_x,actual_x,label_x],[.82,1.04,1.40],color='#B5B5B5',lw=.45)
 ax.text(label_x,1.43,row.gene_name,rotation=60,ha='left',va='bottom',fontsize=6.3,fontweight='normal',fontstyle='italic',color='#888888',clip_on=False)
ax.text(length,-.45,'ddPCR target loci',ha='right',fontsize=6.8,color='#555555')
ax.set_xticks([0,10,20,30,40,50,60]);ax.tick_params(axis='x',labelsize=7,pad=1);ax.set_xlabel('HG002 T2T-Y position (Mb)',fontsize=7.5,labelpad=3);ax.spines[['left','bottom','right','top']].set_visible(False)

stem='ch03_main_ht1080_y_ideogram_standalone'
for ext in ['pdf','svg','png']:
    fig.savefig(FIG/(stem+'.'+ext),dpi=600,facecolor='white',bbox_inches='tight',pad_inches=.08)
plt.close(fig)
print(stem)
