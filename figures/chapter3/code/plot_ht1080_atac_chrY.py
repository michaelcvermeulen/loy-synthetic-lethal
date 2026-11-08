from pathlib import Path
import sys,numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
F=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(F.parents[1]/'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
apply_style(8)
D=F/'source_data/ht1080_atac_chrY'
d=pd.read_csv(D/'coverage_10kb.csv');par=pd.read_csv(D/'chrX_PAR1_10kb.csv');dar=pd.read_csv(D/'DAR_expanded_Y_PAR.csv')
genes=pd.read_csv(D/'labelled_genes_hg38.csv')
par_genes=pd.DataFrame([{'gene_name':'P2RY8','start':1462572,'end':1537187},{'gene_name':'AKAP17A','start':1591604,'end':1602532}])
genes=pd.concat([genes,par_genes],ignore_index=True).sort_values('start')
fig=plt.figure(figsize=(9.3,6.7))
blue,red,gold='#3C5488','#E64B35','#AC7918'
samples=['CTRL_11','CTRL_12','CTRL_21','LOY_11','LOY_12','LOY_21']
for i,s in enumerate(samples):
 ax=fig.add_axes([.13,.83-i*.104,.43,.083])
 q=d[d['sample'].eq(s)];x=(q.start+q.end)/2e6
 p=par[par['sample'].eq(s)];px=(p.start+p.end)/2e6
 col=blue if i<3 else red
 ax.axvspan(.01,2.781479,color=gold,alpha=.13,zorder=0)
 ax.axvline(2.781479,color=gold,ls='--',lw=.6)
 ax.fill_between(x,q.mean_fragment_CPM,color=col,lw=.35)
 ax.fill_between(px,p.mean_fragment_CPM,color=gold,lw=.35)
 ax.set_xlim(0,27);ax.set_ylim(0,.75 if i<3 else .32)
 ax.set_yticks([0,.75] if i<3 else [0,.3]);ax.tick_params(labelsize=7,pad=2,length=2)
 ax.set_xticks([0,5,10,15,20,25]);ax.tick_params(labelbottom=i==5)
 ax.text(.98,.72,s.replace('CTRL','sgCTRL').replace('LOY','Y-depleted'),transform=ax.transAxes,ha='right',va='center',color=col,fontsize=7.5)
 if i==0:ax.text(1.4,.64,'PAR1',color=gold,ha='center',fontsize=7)
ga=fig.add_axes([.13,.135,.43,.12]);ga.set_xlim(0,27);ga.set_ylim(0,1)
ga.axhline(.08,color='#999999',lw=.5)
for k,(r,lx) in enumerate(zip(genes.itertuples(),np.linspace(.7,26,len(genes)))):
 mid=(r.start+r.end)/2e6
 gene_color=gold if r.gene_name in ['P2RY8','AKAP17A'] else '#333333'
 ga.plot([(r.start-1)/1e6,r.end/1e6],[.08,.08],color=gene_color,lw=2.5)
 ga.annotate(r.gene_name,xy=(mid,.08),xytext=(lx,.35 if k%2==0 else .72),ha='center',va='bottom',rotation=55,fontsize=6.5,color=gene_color,arrowprops=dict(arrowstyle='-',color='#777777',lw=.45))
ga.set_yticks([]);ga.spines['left'].set_visible(False);ga.set_xlabel('Chromosome Y position (Mb; hg38)',labelpad=3)
ga.set_xticks([0,5,10,15,20,25]);ga.tick_params(labelsize=7)
fig.text(.13,.945,'a',fontsize=12,weight='semibold')
fig.text(.065,.62,'Fragment coverage (CPM)',rotation=90,ha='center',va='center',fontsize=8)
va=fig.add_axes([.68,.31,.30,.60])
other=dar[dar.region.eq('Other')]
va.scatter(other.logFC,-np.log10(other.FDR.clip(lower=1e-300)),s=3,c='#B5B5B5',alpha=.4,rasterized=True)
for region,col,marker in [('chrY non-PAR','#7C3F98','o'),('chrX PAR',gold,'^')]:
 q=dar[dar.region.eq(region)]
 va.scatter(q.logFC,-np.log10(q.FDR.clip(lower=1e-300)),s=18,c=col,marker=marker,edgecolor='white',lw=.25,zorder=3)
va.axhline(-np.log10(.05),ls='--',color='#777777',lw=.6);va.axvline(0,color='#999999',lw=.5)
va.set_xlabel('log2 fold change\n(Y-depleted / sgCTRL)');va.set_ylabel('-log10 FDR');va.set_xticks([-10,-8,-6,-4,-2,0,2])
fig.text(.64,.945,'b',fontsize=12,weight='semibold')
# Compact volcano and four locus views; retain source statistics unchanged.
fig.set_size_inches(10,7.0)
for ax in list(fig.axes):
 if ax is va:continue
 pos=ax.get_position();ax.set_position([.08,pos.y0,.41,pos.height])
for text in fig.texts:
 x,y=text.get_position()
 if text.get_text()=='b':text.set_position((.53,.945))
 else:text.set_position((.08 if text.get_text()=='a' else .025,y))
va.set_position([.56,.655,.42,.27])
va.set_xlabel('log2 fold change (Y-depleted / sgCTRL)',fontsize=7)
va.tick_params(labelsize=7)
handles=[Line2D([],[],marker='o',color='#B5B5B5',ls='',label='Other'),Line2D([],[],marker='o',color='#7C3F98',ls='',label='chrY'),Line2D([],[],marker='^',color=gold,ls='',label='PAR')]
va.legend(handles=handles,loc='upper right',fontsize=6.5,handletextpad=.3)
fig.text(.56,.582,'FDR < 0.05: 88/89 chrY and 6/41 PAR peaks lower',fontsize=7)
loci=pd.read_csv(D/'locus_coverage_100bp.csv')
for i,gene in enumerate(['P2RY8','AKAP17A','RPS4Y1','DDX3Y']):
 x0=.56+(i%2)*.225;y0=.355-(i//2)*.215
 ax=fig.add_axes([x0,y0,.195,.135]);q=loci[loci.gene.eq(gene)]
 for prefix,col in [('CTRL',blue),('LOY',red)]:
  v=q[q['sample'].str.startswith(prefix)].groupby(['start','end']).CPM.agg(['mean','min','max']).reset_index()
  xx=(v.start+v.end)/2e6
  ax.fill_between(xx,v['min'],v['max'],color=col,alpha=.18,lw=0)
  ax.plot(xx,v['mean'],color=col,lw=.8)
 ymax=ax.get_ylim()[1]
 ax.set_ylim(-ymax*.12,ymax)
 ax.plot([(q.gene_start.iloc[0]-1)/1e6,q.gene_end.iloc[0]/1e6],[-ymax*.065]*2,color='#333333',lw=2)
 ax.set_xlim(q.start.min()/1e6,q.end.max()/1e6)
 ax.tick_params(labelsize=6.5,length=2,pad=2)
 from matplotlib.ticker import MaxNLocator,ScalarFormatter
 ax.xaxis.set_major_locator(MaxNLocator(2));ax.xaxis.set_major_formatter(ScalarFormatter(useOffset=False))
 ax.yaxis.set_major_locator(MaxNLocator(3,prune='lower'))
 ax.set_ylabel('CPM',fontsize=7,labelpad=2)
 ax.set_xlabel(q.chrom.iloc[0]+' position (Mb)',fontsize=7,labelpad=2)
 ax.text(0,1.08,gene+(' (PAR1)' if i<2 else ''),transform=ax.transAxes,fontsize=8,fontstyle='italic')
 ax.text(-.22,1.08,'cdef'[i],transform=ax.transAxes,fontsize=11,weight='semibold')
 if gene in ['P2RY8','AKAP17A','RPS4Y1','DDX3Y']:
  candidates=q[q['sample'].str.startswith('CTRL') & (q.start>=q.gene_start.iloc[0]) & (q.end<=q.gene_end.iloc[0])].groupby('start').CPM.mean()
  centre=int(candidates.idxmax())+50;lo=centre-1500;hi=centre+1500
  ins=ax.inset_axes([.48,.48,.49,.50])
  for prefix,col in [('CTRL',blue),('LOY',red)]:
   vv=q[q['sample'].str.startswith(prefix)&(q.start>=lo)&(q.end<=hi)].groupby(['start','end']).CPM.agg(['mean','min','max']).reset_index()
   xx=(vv.start+vv.end)/2e6
   ins.fill_between(xx,vv['min'],vv['max'],color=col,alpha=.18,lw=0)
   ins.plot(xx,vv['mean'],color=col,lw=.8)
  ins.set_xlim(lo/1e6,hi/1e6);ins.set_ylim(bottom=0)
  ins.tick_params(labelsize=5,length=1,pad=1)
  ins.xaxis.set_major_locator(MaxNLocator(2));ins.xaxis.set_major_formatter(ScalarFormatter(useOffset=False))
  ins.yaxis.set_major_locator(MaxNLocator(2))
  ax.axvspan(lo/1e6,hi/1e6,color='#999999',alpha=.15,zorder=0)
  ins.text(.5,1.06,'3-kb peak inset',transform=ins.transAxes,ha='center',fontsize=6)
  (D/f'{gene}_inset_coordinates.txt').write_text(f'hg38 {q.chrom.iloc[0]}:{lo}-{hi} (0-based); strongest sgCTRL mean 100-bp bin within annotated {gene} gene span.\n')

fig.legend(handles=[Line2D([],[],color=blue,label='sgCTRL'),Line2D([],[],color=red,label='Y-depleted')],bbox_to_anchor=(.56,.545),loc='upper left',ncol=2,fontsize=7,borderaxespad=0)
fig.text(.56,.075,'Lines: mean; shading: replicate range (n = 3/group).',fontsize=7)
stem='ch03_supp_ht1080_atac_Y_PAR'
for ext in ['pdf','svg','png']:fig.savefig(F/(stem+'.'+ext),dpi=400)
print(stem)
