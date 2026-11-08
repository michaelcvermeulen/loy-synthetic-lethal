from pathlib import Path
import sys,gzip,json,hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
ROOT=Path(__file__).resolve().parents[3];F=ROOT/'figures/chapter3';W=ROOT.parent/'WGS'
sys.path.insert(0,str(ROOT/'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
apply_style(8)
CELL=sys.argv[1] if len(sys.argv)>1 else 'DU145'
CTRL={'HT1080':'HT1080-sgCTRL-D0','A549':'A549-CTRL2','DU145':'LNCaP-sgCTRL'}[CELL]
OUT=F/f'source_data/{CELL.lower()}_Yp_10kb_filters';OUT.mkdir(exist_ok=True)
ref=ROOT.parent.parent/'data/chm13v2.0_chrY.fa'
with ref.open() as h:seq=''.join(l.strip() for l in h if not l.startswith('>')).upper()
assert len(seq)==62460029,len(seq)
intervals=[]
with gzip.open(W/'difficulty_regions_T2T.gz','rt') as h:
 for l in h:
  if l.startswith('chrY\t'):
   _,a,b,*_=l.split();intervals.append((int(a),int(b)))
merged=[]
for a,b in sorted(intervals):
 if merged and a<=merged[-1][1]:merged[-1][1]=max(b,merged[-1][1])
 else:merged.append([a,b])
meta=pd.read_csv(F/'source_data/du145_drug_optimisation/lpwgs/Y_content_summary_depth.csv')
if CELL!='DU145':
 samples=[CTRL,{'HT1080':'HT1080-LOY-D0','A549':'A549-LOY'}[CELL]]
 metadata=[]
 for sample in samples:
  summary=pd.read_csv(W/f'coverage/{sample}.genome.100kb.mosdepth.summary.txt',sep='\t')
  autosomes=summary[summary.chrom.isin(['chr'+str(i) for i in range(1,23)])]
  assert len(autosomes)==22
  metadata.append(dict(sample=sample,condition='sgCTRL' if sample==CTRL else 'Y-depleted',autosome_depth=autosomes.bases.sum()/autosomes.length.sum()))
 meta=pd.DataFrame(metadata)
meta.to_csv(OUT/'normalisation.csv',index=False)

rows=[]
for m in meta.itertuples():
 with gzip.open(W/f'coverage/{m.sample}.chrY.regions.bed.gz','rt') as h:
  for l in h:
   ch,a,b,depth=l.split()[:4];a=int(a);b=int(b)
   if a>=10883085:break
   assert b-a==10000
   bases=seq[a:b];n=sum(bases.count(c) for c in 'ACGT');gc=100*(bases.count('G')+bases.count('C'))/n if n else np.nan
   bl=100*sum(max(0,min(b,y)-max(a,x)) for x,y in merged if x<b and y>a)/(b-a)
   rows.append(dict(sample=m.sample,condition=m.condition,start=a,end=b,mid=(a+b)/2e6,depth=float(depth),relative_depth=float(depth)/m.autosome_depth,gc=gc,blacklist_pct=bl,blacklist_flag=bl>=20))
d=pd.DataFrame(rows);d.to_csv(OUT/'all_bins.csv',index=False)
fig=plt.figure(figsize=(10,6.7));gs=fig.add_gridspec(3,2,height_ratios=[1,1.5,.3],left=.075,right=.90,top=.88,bottom=.10,hspace=.40,wspace=.15)
axes=[];ymax=3.0
assert d.relative_depth.max()<=ymax
for row in range(2):
 q=d[d['sample'].eq(CTRL) if row==0 else ~d['sample'].eq(CTRL)]
 for col in range(2):
  ax=fig.add_subplot(gs[row,col]);axes.append(ax)
  ax.axvspan(0,2.45832,color='#EAF6F7');ax.axvspan(10.56575,10.883085,color='#EEEEEE');ax.axvline(2.45832,color='#999999',ls=':',lw=.6)
  if col==0:
   for flag,color in [(False,'#555555'),(True,'#D55E00')]:
    z=q[q.blacklist_flag==flag];ax.scatter(z.mid,z.relative_depth,c=color,s=4,alpha=.5 if row==0 or CELL!='DU145' else .25,lw=0,rasterized=True)
  else:sc=ax.scatter(q.mid,q.relative_depth,c=q.gc,cmap='viridis',vmin=20,vmax=80,s=4,alpha=.75 if row==0 or CELL!='DU145' else .4,lw=0,rasterized=True)
  ax.set(xlim=(0,10.883085),ylim=(-.03,ymax),yticks=[0,1,2,3]);ax.tick_params(labelbottom=False)
  if col==0:ax.set_ylabel('Depth / mean autosomal depth')
  else:ax.tick_params(labelleft=False)
  ax.text(.02,.95,'sgCTRL' if row==0 else ('All seven edited samples' if CELL=='DU145' else 'Y-depleted'),transform=ax.transAxes,va='top',fontsize=8)
  if row==0:ax.set_title('Blacklist status' if col==0 else 'GC content',fontsize=10,pad=10)
fig.suptitle(CELL,fontsize=11,x=.90,y=.98,ha='right')
fig.legend(handles=[Line2D([],[],marker='o',ls='',color=c,markersize=4,label=l) for c,l in [('#555555','<20% blacklist overlap'),('#D55E00','≥20% blacklist overlap')]],loc='upper left',bbox_to_anchor=(.07,.985),frameon=False,ncol=2,fontsize=8)
cax=fig.add_axes([.925,.36,.013,.4]);cb=fig.colorbar(sc,cax=cax);cb.set_label('GC (% of A/C/G/T bases)');cb.set_ticks([20,40,60,80]);cb.set_ticklabels(['≤20','40','60','≥80'])
genes=pd.concat([pd.read_csv(F/'source_data/du145_y_expression_confirmation/gene_mapping_plotted.csv'),pd.read_csv(F/'source_data/ht1080_lpwgs_confirmation/ddpcr_target_loci.csv')]).drop_duplicates('gene_name')
for col in range(2):
 ax=fig.add_subplot(gs[2,col]);ax.axvspan(0,2.45832,ymin=.55,ymax=.85,color='#41B6C4');ax.plot([2.45832,10.56575],[.7,.7],c='#555555',lw=1)
 ax.text(1.2,.96,'PAR1',ha='center',fontsize=7);ax.text(6.5,.96,'MSY Yp',ha='center',fontsize=7)
 for gene in genes[genes.gene_name.isin(['RPS4Y1','ZFY','AMELY'])].itertuples():
  x=(gene.t2t_start+gene.t2t_end)/2e6;g=gene.gene_name
  ax.plot([x,x],[.62,.32],c='#777777',lw=.6);ax.text(x,.15,g,ha='right' if g=='RPS4Y1' else 'left',fontsize=6,fontstyle='italic')
 ax.set(xlim=(0,10.883085),ylim=(0,1.3),yticks=[],xlabel='HG002 T2T-Y position (Mb)');ax.spines[['left','top','right']].set_visible(False)
fig.text(.075,.02,'10-kb bins; all bins retained. GC shown continuously (no established GC pass/fail cutoff). No new GC correction or segmentation.',fontsize=7)
stem=f'ch03_supp_{CELL.lower()}_Yp_10kb_filter_draft'
for ext in ['pdf','svg','png']:fig.savefig(F/f'{stem}.{ext}',dpi=300)
fig.savefig(F/'previews'/f'{stem}.png',dpi=150)
(OUT/'README.md').write_text(f'{CELL} exploratory 10-kb coverage view. Primary controls and Y-depleted samples for HT1080/A549; all seven edited conditions for DU145. Original selected figures unchanged. Coverage is divided by length-weighted mean autosomal depth; see normalisation.csv. No QDNAseq correction or segmentation. All bins retained. Orange: >=20% overlap with union of T2T difficulty regions; this is the blacklist component only, not full QDNAseq filtering. GC is calculated from the local CHM13v2 chrY FASTA and displayed continuously; no GC failure cutoff is inferred. Coordinates zero-based, same Yp window as S01. all_bins.csv records every point. A549-BC14 is not part of the primary comparison.',encoding='utf-8')
print(d.groupby('sample').size().to_string());print('max relative depth',d.relative_depth.max());print('wrote',stem)
