"""A549 primary nuclear-gene RINe model, with condition-only sensitivity."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from adjustText import adjust_text
ROOT=Path(__file__).resolve().parents[3];OUT=ROOT/'figures/chapter3'
DATA=OUT/'source_data/a549_primary_nuclear_RINe'
SRC=OUT/'source_data/a549_rna_thesis'
sys.path.insert(0,str(ROOT/'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
apply_style(8)
BLUE,RED,TEAL,BROWN='#3C5488','#E64B35','#00A087','#7E6148'
de=pd.read_csv(DATA/'DE_nuclear_RINe_adjusted.csv');un=pd.read_csv(DATA/'DE_nuclear_unadjusted.csv')
assert len(de)==22843 and de.gene_id.is_unique and (de.FDR<.05).sum()==27
assert set(de.gene_id)==set(un.gene_id)
xlim=(min(-10,np.floor(min(de.logFC.min(),un.logFC.min()))-1),max(5,np.ceil(max(de.logFC.max(),un.logFC.max()))+1))
ylim=(-.15,max(-np.log10(de.PValue.min()),-np.log10(un.PValue.min()))+1.1)
def volcano(ax,d,genes):
 y=d.chrom.eq('chrY') & ~d.PAR;other=~y & ~d.PAR;sig=d.FDR<.05
 for mask,color,size in [(other & ~sig,'#C8C8C8',3),(other & sig & (d.logFC<0),BLUE,13),(other & sig & (d.logFC>0),RED,13),(d.PAR,TEAL,17),(y,BROWN,18)]:
  q=d[mask];ax.scatter(q.logFC,-np.log10(q.PValue),s=size,c=color,alpha=.6 if size==3 else 1,linewidths=.25,edgecolors='white' if size>3 else 'none',rasterized=size==3,zorder=2 if size==3 else 4)
 ax.set(xlim=xlim,ylim=ylim,xlabel='log₂FC (Y-depleted / sgCTRL)',ylabel='−log₁₀(nominal P)')
 for x in [-1,1]:ax.axvline(x,c='#AAAAAA',ls=':',lw=.5)
 cut=d.loc[sig,'PValue'].max()
 if np.isfinite(cut):ax.axhline(-np.log10(cut),c='#666666',ls='--',lw=.6)
 texts=[]
 for gene in genes:
  q=d[d.gene_name.eq(gene)].iloc[0];texts.append(ax.text(q.logFC,-np.log10(q.PValue),gene,fontsize=6))
 np.random.seed(42);adjust_text(texts,ax=ax,arrowprops=dict(arrowstyle='-',color='#888888',lw=.45))
 n=int((sig & other).sum());ax.text(.98,.98,f'{sig.sum()} significant\n{n} non-Y/non-PAR',transform=ax.transAxes,ha='right',va='top',fontsize=7)
def save(fig,stem):
 for ext in ['pdf','svg','png']:fig.savefig(OUT/f'{stem}.{ext}',dpi=300,facecolor='white')
 plt.close(fig)
from make_a549_global_pathways_combined import render as render_combined
render_combined()
# The standalone S02 comparison is retained in combined S01 panels d and e.
from make_a549_qc_sensitivity_combined import render as render_supplement
render_supplement()
