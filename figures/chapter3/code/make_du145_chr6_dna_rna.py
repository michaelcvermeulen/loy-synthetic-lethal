"""Exploratory chromosome 6 DNA/RNA comparison; no model refitting."""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from du145_rna_common import FIG,T,LPWGS,BLUE,RED,export
from adjustText import adjust_text
stem='ch03_supp_du145_chr6_dna_rna'
g=pd.read_csv(T/'CN_gene_level.csv')
g=g[g.t2t_chrom.eq('chr6') & g.t2t_mid.notna()].copy()
g['Mb']=g.t2t_mid/1e6
g['histone']=g.gene_name.str.match(r'^(H[234](?:AC|BC|C)\d+|H1-\d+)$',na=False)
b=pd.read_csv(LPWGS/'qdnaseq_bins.csv')
b=b[b.chromosome.eq('chr6') & b['use'].eq(True)].copy()
w=b.pivot(index=['start','end','midpoint'],columns='sample',values=['log2','log2_seg']).reset_index()
w.columns=['_'.join(filter(None,x)) for x in w.columns]
w['DNA_bin_delta']=w['log2_LNCaP-LOY-2X']-w['log2_LNCaP-sgCTRL']
w['DNA_segment_delta']=w['log2_seg_LNCaP-LOY-2X']-w['log2_seg_LNCaP-sgCTRL']
h=g[g.histone & g.t2t_arm.eq('chr6p')]
lo,hi=h.Mb.min(),h.Mb.max()
zoom=(np.floor(lo)-.5,np.ceil(hi)+.5)
fig,axs=plt.subplots(3,2,figsize=(7.2,6.0),gridspec_kw={'width_ratios':[1.4,1], 'height_ratios':[1,1,1.5]})
fig.subplots_adjust(left=.12,right=.98,bottom=.10,top=.94,wspace=.24,hspace=.19)
for j,limits in enumerate([(0,b.end.max()/1e6),zoom]):
 ax=axs[0,j]
 for sample,color,label in [('LNCaP-sgCTRL',BLUE,'sgCTRL'),('LNCaP-LOY-2X',RED,'Y-depleted')]:
  z=b[b['sample'].eq(sample)].sort_values('start')
  ax.scatter(z.midpoint,z.log2,s=2,color=color,alpha=.17,rasterized=True,lw=0)
  # Separate segments across excluded bins; no interpolation over missing intervals.
  ax.hlines(z.log2_seg,z.start/1e6,z.end/1e6,color=color,lw=1.5,label=label)
 ax.set_ylim(-1.5,1.5)
 ax=axs[1,j]
 ax.scatter(w.midpoint,w.DNA_bin_delta,s=2,color='#B7B7B7',alpha=.6,lw=0,rasterized=True)
 ax.hlines(w.DNA_segment_delta,w.start/1e6,w.end/1e6,color='#333333',lw=1.5)
 for v in [-.2,.2]: ax.axhline(v,color='#888888',ls='--',lw=.7)
 ax.set_ylim(-.8,.8)
 if j==1: ax.text(.97,.90,f'Histone interval: {h.CN_delta_LOY2X.median():+.3f}',transform=ax.transAxes,ha='right',va='top',fontsize=7)
 ax=axs[2,j]
 other=g[~g.histone];hist=g[g.histone]
 ax.scatter(other.Mb,other.logFC,s=7,color='#AAAAAA',alpha=.6,lw=0,rasterized=True)
 ax.scatter(hist.Mb,hist.logFC,s=16,color='#00A087',lw=.3,edgecolor='white',zorder=4)
 # Draw complete gene range; significant histones use black outlines.
 sig=hist[hist.FDR<.05]
 ax.scatter(sig.Mb,sig.logFC,s=20,facecolors='none',edgecolors='#333333',lw=.6,zorder=5)
 ax.set_ylim(min(-1.3,g.logFC.min()-.15),max(1.3,g.logFC.max()+.2))
 if j==1:
  labels=h.nlargest(6,'logFC')
  texts=[ax.text(r.Mb,r.logFC,r.gene_name,fontsize=6.2) for r in labels.itertuples()]
  ax.set_xlim(*limits)
  adjust_text(texts,ax=ax,arrowprops=dict(arrowstyle='-',color='#888888',lw=.5),expand=(1.2,1.4))
 for i in range(3):
  ax=axs[i,j]
  ax.axhline(0,color='#CCCCCC',lw=.6,zorder=0)
  ax.axvspan(lo,hi,color='#00A087',alpha=.065,zorder=0)
  ax.set_xlim(*limits)
  ax.text(-.10,1.03,'abcdef'[i*2+j],transform=ax.transAxes,weight='semibold',fontsize=11)
  if i<2: ax.tick_params(labelbottom=False)
  if j==1: ax.tick_params(labelleft=False)
 axs[2,j].set_xlabel('Chromosome 6 position (Mb; T2T)')
axs[0,0].set_ylabel('DNA log$_2$ copy ratio')
axs[1,0].set_ylabel('DNA log$_2$ difference\n(Y-depleted \u2212 sgCTRL)')
axs[2,0].set_ylabel('RNA log$_2$ fold change\n(Y-depleted / sgCTRL)')
fig.legend(handles=[Line2D([],[],color=BLUE,label='sgCTRL'),Line2D([],[],color=RED,label='Y-depleted'),Line2D([],[],marker='o',ls='',color='#00A087',label='Histone genes',ms=4)],loc='upper center',ncol=3,frameon=False,fontsize=7,bbox_to_anchor=(.54,1))
export(fig,stem,{'chr6_genes':g,'chr6_DNA_bins':b,'chr6_DNA_differences':w},{'inputs':[T/'CN_gene_level.csv',LPWGS/'qdnaseq_bins.csv'],'notes':'Exploratory draft. Shared T2T coordinates. DNA: usable 100-kb QDNAseq bins and saved segments, sgCTRL and LOY2X (separately collected from RNA). Difference is LOY2X minus sgCTRL. RNA: original passage-blocked voom-limma logFC, 3 passages per population; no refit. Histone protein gene symbols H2/H3/H4 cluster names and H1-number symbols; pseudogenes and antisense excluded. Grey/teal dots show all primary chr6 genes with mapped coordinates. Black rings indicate histone genes with primary FDR <0.05. Full RNA y range, DNA bins outside displayed limits are omitted visually; all values retained in source data. Light teal marks span of mapped chr6p histone genes; right panels zoom around that interval. No segment lines across missing bins.'})
print('Histone interval Mb',lo,hi,'zoom',zoom,'histones',len(h),'higher',sum(h.logFC>0),'FDR05',sum(h.FDR<.05))
