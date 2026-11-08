from pathlib import Path
import sys,textwrap
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
F=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(F.parents[1]/'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
apply_style(7)
d=pd.read_csv(F/'source_data/ht1080_atac_gene_set_context/fgsea_wide.csv')
parts=[]
for coll in ['Hallmark','Reactome','GO:BP']:
 g=d[d.collection==coll]
 keep=(g.padj_ATAC_median<.05)|(g.padj_ATAC_trended<.05)|g.index.isin(g.nsmallest(5,'padj_ATAC_trended').index)
 parts.append(g[keep].sort_values('NES_ATAC_trended'))
s=pd.concat(parts)
out=F/'previews';out.mkdir(exist_ok=True)
s.to_csv(out/'ht1080_promoter_gsea_all_collections_review.csv',index=False)
fig=plt.figure(figsize=(210/25.4,(42+len(s)*6.1)/25.4))
ax=fig.add_axes([.53,.10,.24,.82]);ax.axvline(0,color='#CCCCCC',lw=.7)
y=0
for group in parts:
 coll=group.collection.iloc[0]
 ax.text(-2.7,y-.9,coll,fontsize=8,weight='semibold',ha='right')
 for r in group.itertuples():
  name=r.pathway
  for prefix in ['HALLMARK_','REACTOME_','GOBP_']:name=name.removeprefix(prefix)
  name=name.replace('_',' ').lower()
  ax.text(-2.7,y,'\n'.join(textwrap.wrap(name,60)),ha='right',va='center',fontsize=5.5)
  ax.plot([r.NES_ATAC_median,r.NES_ATAC_trended],[y,y],color='#CCCCCC',lw=.7)
  for v,q,c,m in [(r.NES_ATAC_median,r.padj_ATAC_median,'#555555','o'),(r.NES_ATAC_trended,r.padj_ATAC_trended,'#00A087','s')]:
   ax.scatter(v,y,s=20,marker=m,facecolor=c if q<.05 else 'white',edgecolor=c,lw=.9,zorder=3)
  ax.text(2.8,y,f'{r.padj_ATAC_median:.3f}     {r.padj_ATAC_trended:.3f}',va='center',fontsize=5.5)
  y+=1
 y+=1.6
ax.set_ylim(y-.5,-1.5);ax.set_xlim(-2.5,2.5);ax.set_yticks([]);ax.spines['left'].set_visible(False);ax.set_xlabel('Promoter-accessibility NES')
fig.text(.035,.97,'Promoter GSEA across all three collections',fontsize=10,weight='semibold')
fig.text(.035,.945,'Significant in either analysis, plus five lowest corrected FDRs per collection',fontsize=7)
fig.text(.81,.935,'BH FDR\nPrimary   Corrected',fontsize=6)
fig.legend(handles=[Line2D([],[],marker='o',color='#555555',ls='',label='Primary'),Line2D([],[],marker='s',color='#00A087',ls='',label='Trend-corrected')],loc='lower left',bbox_to_anchor=(.025,.015),ncol=2,frameon=False,fontsize=7)
fig.text(.56,.035,'Filled symbols: FDR <0.05',fontsize=7)
for ext in ['pdf','png']:fig.savefig(out/f'ht1080_promoter_gsea_all_collections_review.{ext}',dpi=220,bbox_inches='tight')
print('Preview saved:',len(s),'sets')
