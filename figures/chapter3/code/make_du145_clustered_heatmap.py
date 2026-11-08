"""Alternative 114-gene heatmap: conditions grouped, clustered genes and samples within condition."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.cluster.hierarchy import linkage,dendrogram,leaves_list
from scipy.spatial.distance import pdist
from du145_rna_common import BLUE,RED,export,T
from du145_reporting import annotate
F=Path(__file__).resolve().parents[1]
de=annotate(pd.read_csv(T/'DE_paired_voom_primary.csv'))
ids=de.loc[de.reporting_subset,'gene_id']
x=pd.read_csv(T/'core_logCPM_all_genes.csv').pivot(index='gene_id',columns='sample_key',values='log2_CPM').loc[ids]
z=x.sub(x.mean(axis=1),axis=0).div(x.std(axis=1,ddof=1),axis=0)
assert z.shape==(114,6) and np.isfinite(z.to_numpy()).all()
# Correlation distance clusters genes by their six-sample expression patterns.
row_tree=linkage(pdist(z.to_numpy(),metric='correlation'),method='average',optimal_ordering=True)
order=leaves_list(row_tree)
col_order=[];col_trees=[]
for group in ['sgCTRL','LOY']:
 cols=[c for c in z.columns if '-'+group+'-' in c]
 tree=linkage(pdist(z[cols].to_numpy().T,metric='euclidean'),method='average',optimal_ordering=True)
 col_trees.append(tree);col_order.extend([cols[i] for i in leaves_list(tree)])
zz=z.iloc[order][col_order];genes=de.set_index('gene_id').gene_name.reindex(zz.index)
fig=plt.figure(figsize=(6.2,11.2))
ax=fig.add_axes([.18,.055,.54,.84]);rd=fig.add_axes([.025,.055,.145,.84])
dendrogram(row_tree,orientation='left',no_labels=True,ax=rd,color_threshold=0,above_threshold_color='#999999')
rd.invert_yaxis();rd.set_axis_off()
im=ax.imshow(zz,aspect='auto',cmap=LinearSegmentedColormap.from_list('rowz',[BLUE,'#F7F7F7',RED]),vmin=-2,vmax=2,interpolation='nearest')
ax.set_yticks(np.arange(114),genes,fontsize=5.4);ax.yaxis.tick_right();ax.tick_params(length=0)
ax.set_xticks(range(6),[c.rsplit('-',1)[-1] for c in col_order],fontsize=7)
for tick,col in zip(ax.get_xticklabels(),[BLUE]*3+[RED]*3):tick.set_color(col)
ax.axvline(2.5,color='white',lw=3)
for spine in ax.spines.values():spine.set_visible(False)
for j,tree in enumerate(col_trees):
 da=fig.add_axes([.18+j*.27,.900,.27,.04])
 dendrogram(tree,ax=da,no_labels=True,color_threshold=0,above_threshold_color=BLUE if j==0 else RED)
 da.set_axis_off()
 fig.text(.315+j*.27,.947,['sgCTRL','Y-depleted'][j],ha='center',color=[BLUE,RED][j],fontsize=9)
cax=fig.add_axes([.37,.973,.20,.007]);cb=fig.colorbar(im,cax=cax,orientation='horizontal',ticks=[-2,0,2]);cb.ax.tick_params(labelsize=6);cb.ax.set_title('Row z-score',fontsize=7,pad=3)
export(fig,'ch03_supp_du145_rna_clustered_heatmap',{'expression_log2CPM':x.loc[zz.index,col_order].reset_index(),'row_z_scores':zz.reset_index(),'gene_order':pd.DataFrame({'gene_id':zz.index,'gene_name':genes.to_numpy()}),'row_linkage':pd.DataFrame(row_tree,columns=['left','right','distance','n'])},{'inputs':[T/'DE_paired_voom_primary.csv',T/'core_logCPM_all_genes.csv'],'notes':'114 genes: FDR<0.05, abs(log2FC)>=0.5, same direction in three passage pairs. Frozen TMM log2CPM, row z-score with sample SD across six libraries. Average-linkage gene clustering with correlation distance. sgCTRL and LOY remain grouped; samples clustered within each group using Euclidean distance on row z-scores. Display clipped at +/-2. No passage adjustment; descriptive clustering of genes selected from this same comparison.'})
