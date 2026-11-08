from pathlib import Path
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[3]
FIG=ROOT/'figures/chapter3'
sys.path.insert(0,str(ROOT/'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
apply_style(8)
D=FIG/'source_data/ht1080_atac_copy_number'
r=pd.read_csv(D/'region_CN_context.csv');r=r[r.chrom=='chr19'].copy();r['Mb']=(r.start+r.end)/2e6
b=pd.read_csv(D/'WGS_hg38_100kb_bins.csv');b=b[(b.chrom=='chr19')&b.use&b.delta_segment_log2.notna()]
summary=r.groupby(np.floor(r.Mb).astype(int)).agg(n=('logFC','size'),median_log2FC=('logFC','median'),mean_log2FC=('logFC','mean'),negative_fraction=('logFC',lambda v:(v<0).mean()))
summary.index.name='bin_start_Mb'
summary.to_csv(FIG/'previews/ht1080_chr19_regional_ATAC_1Mb.csv')
f,axes=plt.subplots(2,1,figsize=(183/25.4,145/25.4),gridspec_kw={'height_ratios':[1,1.2]})
for ax,limits,title in zip(axes,[(0,59.2),(50,59.2)],['Chromosome 19','Distal chromosome 19: 50-59 Mb']):
 d=r[r.Mb.between(*limits)]
 ax.scatter(d.Mb,d.logFC,s=5,color='#BBBBBB',alpha=.7,lw=0,label='All tested peaks',rasterized=True)
 sig=d[d.FDR<.05];ax.scatter(sig.Mb,sig.logFC,s=13,color='#3C5488',lw=0,label='ATAC FDR <0.05',zorder=3)
 ss=summary[(summary.index>=limits[0])&(summary.index<limits[1])]
 ax.plot(ss.index+.5,ss.median_log2FC,color='#222222',marker='o',ms=2,lw=1,label='1 Mb median')
 ax.scatter((b.start0+b.end0)/2e6,b.delta_segment_log2,color='#8C6BB1',s=3,label='Earlier D0 lpWGS',lw=0)
 ax.axhline(0,color='#AAAAAA',lw=.5);ax.axvspan(52,53,color='#F5E9D6',alpha=.6,zorder=-1)
 ax.set(xlim=limits,ylim=(-3.2,2.3),ylabel='ATAC / lpWGS difference (log2)',title=title)
 ax.set_xlabel('Position (Mb, hg38)')
axes[0].legend(ncol=2,loc='upper right',fontsize=6,frameon=False)
axes[1].set_xticks(range(50,60))
axes[1].text(52.5,1.9,'52-53 Mb: 26 peaks\nmedian ATAC log2FC = -0.53',ha='center',va='top',fontsize=6.3)
axes[1].text(56.5,1.9,'56-57 Mb: 22 peaks\nmedian ATAC log2FC = -0.07',ha='center',va='top',fontsize=6.3)
f.tight_layout()
for ext in ['pdf','png','svg']:
 f.savefig(FIG/f'previews/ht1080_chr19_regional_ATAC_review.{ext}',dpi=220,bbox_inches='tight')
print(summary.loc[52:56].to_string())
