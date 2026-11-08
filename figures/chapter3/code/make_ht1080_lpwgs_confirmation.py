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
# The adopted regional atlas is rendered without rewriting the source tables.
if '--atlas-only' in sys.argv:
 import subprocess
 subprocess.run([sys.executable,str(Path(__file__).with_name('make_regional_y_atlas.py')),'HT1080'],check=True)
 raise SystemExit(0)
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
exported=[]
validation_path=W/'figure_validation.json'
prior_validation=json.loads(validation_path.read_text()) if validation_path.exists() else {}
def save(fig,stem):
 if stem=='ch03_main_ht1080_lpwgs_chrY_atlas':
  plt.close(fig)
  import subprocess
  subprocess.run([sys.executable,str(Path(__file__).with_name('make_regional_y_atlas.py')),'HT1080'],check=True)
  exported.append(stem)
  return
 if stem=='ch03_main_ht1080_lpwgs_genomewide':
  plt.close(fig)
  if '--atlas-only' not in sys.argv:
   import subprocess
   subprocess.run([sys.executable,str(Path(__file__).with_name('make_ht1080_lpwgs_difference.py'))],check=True)
   exported.append(stem)
  return
 if '--atlas-only' in sys.argv and stem!='ch03_main_ht1080_lpwgs_chrY_atlas':
  plt.close(fig);return
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
d=pd.read_csv(W/'qdnaseq_bins.csv');d=d[d.use & d.log2_seg.notna()].copy();d['rel_raw']=2**d.log2;d['rel_seg']=2**d.log2_seg
assert set(d.label)=={'sgCTRL','LOY'}
chroms=['chr'+str(c) for c in list(range(1,23))+['X','Y']]
# Genome-wide profiles, retaining the original bin and segment values.
fig=plt.figure(figsize=(183/25.4,145/25.4));gs=fig.add_gridspec(3,1,height_ratios=[1,1,.8],left=.10,right=.975,top=.94,bottom=.09,hspace=.43)
axes=[fig.add_subplot(gs[i]) for i in range(3)]
for ax,condition,letter in zip(axes[:2],['sgCTRL','LOY'],['a','b']):
 g=d[d.label.eq(condition)]
 mids=[]
 for j,ch in enumerate(chroms):
  h=g[g.chromosome.eq(ch)].sort_values('genome_x');lo,hi=h.genome_x.min(),h.genome_x.max();mids.append((lo+hi)/2)
  if j%2:ax.axvspan(lo,hi,color='#F2F2F2',zorder=0)
  ax.scatter(h.genome_x.where(h.log2.between(-3,1.5)),h.log2.where(h.log2.between(-3,1.5)),s=.7,color='#A0A0A0',alpha=.45,lw=0,rasterized=True)
  ax.plot(h.genome_x,h.log2_seg.where(h.log2_seg.between(-3,1.5)),color=COL[condition],lw=.75)
  clipped=h[h.log2_seg.lt(-3)];ax.scatter(clipped.genome_x,np.full(len(clipped),-3),marker='_',s=1,color=COL[condition],lw=.5)
 ax.axhline(0,color='#888888',lw=.5);ax.set(xlim=(d.genome_x.min(),d.genome_x.max()),ylim=(-3.05,1.55),yticks=[-3,-2,-1,0,1])
 ax.set_xticks(mids,[c[3:] for c in chroms]);ax.tick_params(axis='x',labelsize=6.3,pad=2,length=2)
 ax.set_ylabel('log$_2$(relative CN)',fontsize=8);ax.set_title(WGS_LABEL[condition],loc='left',fontsize=9,color=COL[condition]);tag(ax,letter)
axes[1].set_xlabel('Chromosome')
ax=axes[2]
for condition in COL:
 h=d[d.label.eq(condition)&d.chromosome.eq('chr13')].sort_values('midpoint')
 ax.scatter(h.midpoint,np.clip(h.log2,-.5,1),s=1,color=COL[condition],alpha=.22,lw=0,rasterized=True)
 ax.plot(h.midpoint,np.clip(h.log2_seg,-.5,1),color=COL[condition],lw=1,label=LABEL[condition])
ax.set(xlim=(18.1,113.567),ylim=(-.45,.95),xlabel='Chromosome 13 position (Mb)',ylabel='log$_2$(relative CN)');ax.axhline(0,color='#888888',lw=.5);ax.legend(frameon=False,ncol=2,loc='upper right',fontsize=7)
ax.set_title('Chromosome 13q copy-number difference',loc='left',fontsize=9);tag(ax,'c')
ax.text(65,.52,'0.480',color=COL['LOY'],fontsize=7);ax.text(65,.15,'0.112',color=COL['sgCTRL'],fontsize=7)
save(fig,'ch03_main_ht1080_lpwgs_genomewide')
long=pd.read_csv(R/'HT1080_Y_PAR_per_sample.csv');mapping=pd.read_csv(R/'gene_mapping_plotted.csv').sort_values('t2t_start').reset_index(drop=True);meta=pd.read_csv(R/'sample_metadata.csv')
assert len(long)==132 and long.groupby('gene_name')['sample'].nunique().eq(12).all()
assert np.allclose(long.cpm,long.raw_count/long.library_size_all_genes*1e6)
assert np.allclose(long.log2_cpm_plus_1,np.log2(long.cpm+1))
summary=long.groupby(['condition','gene_name']).log2_cpm_plus_1.mean().unstack(0)
windows=long.groupby(['batch','gene_name','condition']).cpm.mean().unstack('condition');windows['edited_to_control_mean_cpm_ratio']=windows.LOY/windows.sgCTRL;windows.to_csv(R/'window_mean_cpm_ratios.csv')
genes=mapping.gene_name.tolist();n=len(genes)
# Normalize unlogged CPM to the gene's control mean within each sampling window.
reference=long[long.condition.eq('sgCTRL')].groupby(['batch','gene_name']).cpm.mean()
relative=long.join(reference.rename('control_mean_cpm'),on=['batch','gene_name'])
relative['relative_expression_percent']=100*relative.cpm/relative.control_mean_cpm.where(relative.control_mean_cpm.gt(0))
relative.to_csv(R/'relative_expression_per_sample.csv',index=False)
relative_summary=relative.groupby(['condition','gene_name']).relative_expression_percent.mean().unstack(0)
relative_summary.to_csv(R/'relative_expression_group_means.csv')
assert np.allclose(relative_summary.sgCTRL.dropna(),100)
assert np.nanmax(relative.relative_expression_percent)<160


# Regional Y atlas: region medians from unsegmented corrected values.
y=d[d.chromosome.eq('chrY')].copy();computed=[]
for name,a,b in regions:
 for condition in COL:
  v=y[y.label.eq(condition)&y.midpoint.ge(a)&y.midpoint.lt(b)].rel_raw
  computed.append(dict(region=name,condition=condition,n_bins=len(v),median_relative_cn=float(v.median())))
stats=pd.DataFrame(computed);stats.to_csv(W/'regional_medians_verified.csv',index=False)
old=pd.read_csv(W/'figure_chrY_region_stats.tsv',sep='\t')
for row in stats.itertuples():
 ref=old[old.region.eq(row.region)&old.label.str.contains('sgCTRL' if row.condition=='sgCTRL' else 'LOY')]
 assert len(ref)==1 and np.isclose(row.median_relative_cn,ref.median_rel_cn.iloc[0]) and row.n_bins==ref.n_bins.iloc[0]
fig=plt.figure(figsize=(183/25.4,205/25.4));gs=fig.add_gridspec(5,1,height_ratios=[1,1,.60,.72,1.22],left=.105,right=.90,top=.96,bottom=.10,hspace=.49)
for i,condition in enumerate(COL):
 ax=fig.add_subplot(gs[i]);g=y[y.label.eq(condition)]
 for name,a,b in regions:
  if name in ['centromere','Yq12 heterochromatin']:ax.axvspan(a,b,color='#EAEAEA')
  elif name.startswith('PAR'):ax.axvspan(a,b,color='#EAF6F7')
 ax.scatter(g.midpoint,np.minimum(g.rel_raw,1.3),s=4,color='#888888',alpha=.6,lw=0)
 high=g[g.rel_raw>1.3];ax.scatter(high.midpoint,np.full(len(high),1.3),s=8,marker='^',color='#555555',lw=0)
 for row in stats[stats.condition.eq(condition)].itertuples():
  if row.region in ['centromere','Yq12 heterochromatin']:continue
  _,a,b=next(r for r in regions if r[0]==row.region);ax.plot([a,b],[row.median_relative_cn]*2,color=COL[condition],lw=2)
 ax.set(xlim=(0,length),ylim=(-.015,1.35),yticks=[0,.5,1],ylabel='Relative CN');ax.tick_params(axis='x',labelbottom=False)
 ax.set_title(WGS_LABEL[condition],loc='left',fontsize=9,color=COL[condition]);tag(ax,chr(97+i))
 if i==0:ax.text((q0+p2)/2,1.12,'Yq12: limited mappability',ha='center',fontsize=8,color='#666666')
# Chromosome-aligned region strips preserve PAR and heterochromatin context.
ax=fig.add_subplot(gs[2]);norm=matplotlib.colors.Normalize(0,.70);cmap=plt.get_cmap('viridis')
short={'PAR1':'PAR1','MSY Yp':'MSY Yp','MSY Yq (euchromatic)':'MSY Yq','Yq12 heterochromatin':'Yq12','PAR2':'PAR2'}
for name,left,right in regions:
 if name=='centromere':
  ax.axvspan(left,right,color='#DDDDDD',zorder=0);continue
 mid=(left+right)/2
 label_x=0 if name=='PAR1' else length if name=='PAR2' else mid
 ax.text(label_x,1.27,short[name],ha='left' if name=='PAR1' else 'right' if name=='PAR2' else 'center',fontsize=7.4,va='bottom')
 for i,condition in enumerate(COL):
  value=float(stats[stats.region.eq(name)&stats.condition.eq(condition)].median_relative_cn.iloc[0]);yc=.79 if i==0 else .25
  ax.add_patch(Rectangle((left,yc-.22),right-left,.44,facecolor=cmap(norm(value)),edgecolor='white',lw=.6))
  if name in ['PAR1','PAR2']:
   # These regions are narrow at chromosome scale; keep values outside the strips.
   ytext=-.27 if i==0 else -.56
   ax.text(label_x,ytext,f'{value:.3f}',ha='left' if name=='PAR1' else 'right',va='center',fontsize=7,color=COL[condition],fontweight='semibold')
  else:
   txt='Floor' if value<=.0010000001 else f'{value:.3f}'
   ax.text(mid,yc,txt,ha='center',va='center',fontsize=7.4,color='white' if value<.40 else '#222222')
 if name in ['PAR1','PAR2']:ax.plot([mid,mid],[0,-.1],color='#777777',lw=.5)
ax.set(xlim=(0,length),ylim=(-.77,1.57));ax.set_yticks([.79,.25],['sgCTRL','Y-depleted']);ax.set_xticks([]);ax.tick_params(length=0,labelsize=7)
for sp in ax.spines.values():sp.set_visible(False)
tag(ax,'c')
# Compact shared scale, with no rescaling between specimens or regions.
pos=ax.get_position();cax=fig.add_axes([.925,pos.y0+.05*pos.height,.012,.70*pos.height]);cb=fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm,cmap=cmap),cax=cax,ticks=[0,.35,.7]);cb.ax.tick_params(labelsize=6);cb.set_label('Relative CN',fontsize=7);cb.outline.set_visible(False)
# Original coordinate scale with the Y-linked ddPCR assay target loci.
ax=fig.add_subplot(gs[3]);ideogram(ax);ax.set_ylim(-1.0,2.55);tag(ax,'d');ax.set_axis_on();ax.set_yticks([])
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
# RNA-seq confirms regional expression patterns; gene columns are evenly spaced.
ax=fig.add_subplot(gs[4]);tag(ax,'e')
for i,row in mapping.iterrows():
 if row.region!='MSY':ax.axvspan(i-.5,i+.5,color='#EAF6F7',zorder=0)
for condition,off in [('sgCTRL',-.19),('LOY',.19)]:
 ax.scatter(np.arange(n)+off,relative_summary.loc[genes,condition],marker='_',s=95,color=COL[condition],linewidths=2,zorder=4)
 for batch,mark,dx in [('August 2025','o',-.055),('May 2026','^',.055)]:
  for passage in [1,2,3]:
   values=relative[relative.condition.eq(condition)&relative.batch.eq(batch)&relative.passage_label.eq(passage)].set_index('gene_name').loc[genes].relative_expression_percent
   ax.scatter(np.arange(n)+off+dx+(passage-2)*.045,values,s=9,facecolor='white',edgecolor=COL[condition],marker=mark,lw=.7,zorder=3)
ax.axhline(100,color='#888888',lw=.7,ls='--',zorder=1)
ax.set(xlim=(-.5,n-.5),ylim=(0,185),yticks=[0,50,100,150],ylabel='Expression relative\nto sgCTRL (%)');ax.set_xticks(range(n),genes,rotation=60,ha='right',fontsize=7,fontstyle='italic');ax.tick_params(axis='x',length=0)
for tick,row in zip(ax.get_xticklabels(),mapping.itertuples()):
 if row.region!='MSY':tick.set_color('#007786')
handles=[Line2D([],[],color=COL[c],lw=4,label=LABEL[c]) for c in COL]+[Line2D([],[],marker=m,linestyle='none',color='#333333',markersize=3,label=b) for b,m in [('Aug 2025','o'),('May 2026','^')]]
ax.legend(handles=handles,ncol=4,frameon=False,fontsize=6.3,loc='upper right',bbox_to_anchor=(.99,.985),columnspacing=.9,handlelength=1.1,borderaxespad=0)
save(fig,'ch03_main_ht1080_lpwgs_chrY_atlas')
# Descriptive RNA values are the original Ben_Thesis values, not DE-normalized.
fig=plt.figure(figsize=(183/25.4,150/25.4));gs=fig.add_gridspec(3,1,height_ratios=[.85,1.05,1.3],left=.16,right=.94,top=.95,bottom=.14,hspace=.53)
ax=fig.add_subplot(gs[0]);ideogram(ax);ax.set_ylim(-2.4,1.4)
for i,row in mapping.iterrows():
 pos=(row.t2t_start+row.t2t_end)/2e6;target=(i+.5)*length/n;ax.plot([pos,pos,target],[.2,.02,-.8],color='#8191A8',lw=.4);ax.text(target,-.88,row.gene_name,rotation=55,ha='right',va='top',fontsize=6.4,fontstyle='italic',color='#007786' if row.region!='MSY' else '#222222')
ax.text(p1/2,1.02,'PAR1',ha='center',fontsize=7,color='#007786');ax.text(length,1.02,'PAR2',ha='right',fontsize=7,color='#007786');ax.text(6,1.02,'Yp',fontsize=8);ax.text(36,1.02,'Yq',fontsize=8)
ax.text((c0+c1)/2,1.02,'Centromere',fontsize=7,ha='left');tag(ax,'a')
ax=fig.add_subplot(gs[1])
for i,row in mapping.iterrows():
 if row.region!='MSY':ax.axvspan(i-.5,i+.5,color='#EAF6F7',zorder=0)
for condition,off in [('sgCTRL',-.19),('LOY',.19)]:
 ax.bar(np.arange(n)+off,summary.loc[genes,condition],width=.34,color=COL[condition],label=LABEL[condition],zorder=2)
 for batch,mark,dx in [('August 2025','o',-.055),('May 2026','^',.055)]:
  for passage in [1,2,3]:
   v=long[long.condition.eq(condition)&long.batch.eq(batch)&long.passage_label.eq(passage)].set_index('gene_name').loc[genes].log2_cpm_plus_1
   ax.scatter(np.arange(n)+off+dx+(passage-2)*.045,v,s=5,color='#333333',marker=mark,lw=0,zorder=3)
ax.set(xlim=(-.5,n-.5),ylim=(0,11.4),ylabel='log$_2$(CPM + 1)');ax.set_xticks(range(n),genes,rotation=60,ha='right',fontsize=7,fontstyle='italic');ax.tick_params(axis='x',length=0);tag(ax,'b')
handles=[Line2D([],[],color=COL[c],lw=5,label=LABEL[c]) for c in COL]+[Line2D([],[],marker=m,linestyle='none',color='#333333',markersize=3,label=b) for b,m in [('Aug 2025','o'),('May 2026','^')]]
fig.legend(handles=handles,ncol=4,frameon=False,fontsize=6.7,loc='upper center',bbox_to_anchor=(.55,.725),columnspacing=1,handlelength=1.1)
ax=fig.add_subplot(gs[2]);matrix=long.pivot(index='sample',columns='gene_name',values='log2_cpm_plus_1').loc[meta['sample'],genes]
im=ax.imshow(matrix,aspect='auto',cmap='magma',vmin=0,vmax=11)
labels=[f'{"Aug" if row.batch=="August 2025" else "May"} {"CTRL" if row.condition=="sgCTRL" else "Y-depl."} {row.passage_label}' for row in meta.itertuples()]
ax.set_yticks(range(12),labels,fontsize=6.5);ax.set_xticks(range(n),genes,rotation=60,ha='right',fontsize=7,fontstyle='italic');ax.tick_params(length=0)
for line in [2.5,5.5,8.5]:ax.axhline(line,color='white',lw=2 if line==5.5 else .6)
for i in range(n-1):ax.axvline(i+.5,color='white',lw=.35)
for spine in ax.spines.values():spine.set_visible(False)
tag(ax,'c');cbax=fig.add_axes([.957,.14,.012,.22]);cb=fig.colorbar(im,cax=cbax);cb.set_label('log$_2$(CPM + 1)',fontsize=7);cb.ax.tick_params(labelsize=6);cb.outline.set_visible(False)
save(fig,'ch03_main_ht1080_y_expression_confirmation')
# Exact source-derived 13q segment summary; no absolute CN inference.
d[d.chromosome.eq('chr13')].groupby(['label','log2_seg']).agg(start=('start','min'),end=('end','max'),n_bins=('bin_id','size')).reset_index().to_csv(W/'chr13_segment_summary.csv',index=False)
prior_validation.update({'bin_selection':'use == TRUE and nonmissing segmented log2','original_regional_statistics_match':True,'n_WGS_specimens':2,'RNA_values_unchanged':True,'RNA_original_pools_per_condition':1,'atlas_RNA_normalization':'100 * sample CPM / gene-specific mean sgCTRL CPM within the same RNA sampling window; no logarithmic transformation','atlas_RNA_summary':'mean of six sample percentages per condition; circles August 2025, triangles May 2026; three per condition per window','atlas_RNA_undefined_control_references':int(reference.eq(0).sum()),'atlas_RNA_y_limits':[0,185]})
prior_validation['figures']=list(dict.fromkeys(prior_validation.get('figures',[])+exported))
for stem in exported:
 entry=next((v for v in prior_validation.get('exports',[]) if v['stem']==stem),None)
 if entry is not None:entry['sha256']={ext:hashlib.sha256((FIG/(stem+'.'+ext)).read_bytes()).hexdigest() for ext in ['pdf','svg','png']}
validation_path.write_text(json.dumps(prior_validation,indent=2)+'\n')
print('Exported '+', '.join(exported))
