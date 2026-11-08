"""Combined F01: sampling, Y expression, PCA, DE, six genes and PAR expression."""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy.cluster.hierarchy import dendrogram

FIG=Path(__file__).resolve().parents[1]
# Reuse the unchanged DE drawing and label-selection logic, with new axes only.
source=(FIG/'code/make_ht1080_rna_selected_main.py').read_text(encoding='utf-8')
source=source[:source.index("export(fig,'ch03_main_ht1080_rna_expression_changes'")]
source=source.replace("texts=[a.text", "selected=pd.concat([selected,d[d.gene_name.isin(['UNC5A','ROR2','SDK2','COL1A2','FUCA2','MAGEA1'])]]).drop_duplicates('gene_name')\ntexts=[a.text",1)
source=source.replace("fig,(a,b)=plt.subplots(1,2,figsize=(183/25.4,125/25.4))\nfig.subplots_adjust(left=.10,right=.97,bottom=.15,top=.83,wspace=.40)","fig=plt.figure(figsize=(11.8,7.6))\na=fig.add_axes([.07,.085,.28,.39]);b=fig.add_axes([.41,.085,.25,.39])")
source=source.replace("[(a,'a','Overall expression differences'),(b,'b','Agreement across sampling periods')]","[(a,'d','Overall expression differences'),(b,'e','Agreement across sampling periods')]")
source=source.replace("pos.x0-.06,.905","pos.x0-.035,.50").replace("pos.x0,.905","pos.x0,.50")
source=source.replace("fig.text(.10,.97,","fig.text(.07,.535,")
exec(compile(source,str(FIG/'code/make_ht1080_rna_selected_main.py'),'exec'))
for gene,position in [('COL1A2',(-4.9,3.55)),('FUCA2',(-4.3,2.5))]:
    for label in texts:
        if label.get_text()==gene:
            label.set_visible(False)
    row=d.loc[d.gene_name.eq(gene)].iloc[0]
    a.annotate(gene,xy=(row.logFC,-np.log10(row.FDR)),xytext=position,
               fontsize=5.5,color=blue,ha='right',
               arrowprops=dict(arrowstyle='-',color='#888888',lw=.4))
DATA=FIG/'source_data/ht1080_rna_global_structure'
COL={'CTRL':blue,'LOY':red,'PARENTAL':'#888888'}
SHAPE={'AUG_2025':'o','MAY_2026':'^'}
design=pd.read_csv(DATA/'Fig01_design_QC_d.csv')
pc=pd.read_csv(DATA/'panel_a_adjusted_all15.csv')
variance=pd.read_csv(DATA/'panel_a_variance.csv').query("view == 'adjusted'").set_index('component').variance_fraction

# a: compact sampling scheme, keeping the original population/occasion assignments.
ax=fig.add_axes([.09,.68,.23,.22])
for treatment,y in [('CTRL',1),('LOY',0)]:
    sub=design[design.treatment==treatment].sort_values('occasion')
    ax.plot(sub.occasion,[y]*len(sub),color=COL[treatment],lw=1)
    for q in sub.itertuples():
        ax.scatter(q.occasion,y,color=COL[treatment],marker=SHAPE[q.batch],s=28)
        ax.text(q.occasion,y+.16,str(q.passage_number),ha='center',fontsize=7,color=COL[treatment])
ax.set(xlim=(.5,7.5),ylim=(-.4,1.6))
ax.set_xticks([2,6],['Aug 2025','May 2026'],fontsize=8)
ax.set_yticks([1,0],['sgCTRL','Y-depleted'],fontsize=8)
ax.tick_params(length=0)
for y,label in [(1,'2\u00d7'),(0,'2\u00d7 CP + Q53')]:
    ax.annotate(label,(0,y),xycoords=ax.get_yaxis_transform(),xytext=(-3,-12),textcoords='offset points',ha='right',fontsize=6.5,color='#777777')
ax.spines[:].set_visible(False)
ax.text(.5,-.35,'Same established populations\nThree successive passages per period\n12 libraries; six passage matches',transform=ax.transAxes,ha='center',fontsize=7.5,linespacing=1.4)
fig.text(.025,.925,'a',fontsize=12,weight='semibold')

# b: Y-linked expression across all 15 libraries; former b is now S01g.
from ht1080_expression_heatmaps import draw_y_heatmap
draw_y_heatmap(fig, FIG)

# c: saved all-15 period-adjusted PCA.
ax=fig.add_axes([.715,.655,.26,.26])
for q in pc.itertuples():
    ax.scatter(q.PC1,q.PC2,c=COL[q.treatment],marker=SHAPE[q.batch],s=28,edgecolor='white',lw=.4)
    dx,dy,ha=4,4,'left'
    if q.treatment=='PARENTAL' and q.passage_number==2:dx,dy,ha=-5,6,'right'
    if q.batch=='AUG_2025' and q.treatment=='CTRL' and q.passage_number in [2,3]:dx,dy,ha=-5,-7,'right'
    if q.batch=='MAY_2026' and q.treatment=='LOY' and q.passage_number==3:dx,dy,ha=-5,-7,'right'
    ax.annotate(str(q.passage_number),(q.PC1,q.PC2),xytext=(dx,dy),textcoords='offset points',fontsize=7,color=COL[q.treatment],ha=ha)
ax.set(xlabel=f'PC1 ({100*variance.PC1:.1f}%)',ylabel=f'PC2 ({100*variance.PC2:.1f}%)',xlim=(-110,94),ylim=(-95,92))
ax.set_xticks([-100,-50,0,50]);ax.set_yticks([-75,-50,-25,0,25,50,75])
ax.text(.02,.02,'Period-adjusted; 15 libraries',transform=ax.transAxes,fontsize=7)
fig.text(.665,.925,'c',fontsize=12,weight='semibold')
handles=[Patch(facecolor=COL[k],label=v) for k,v in [('CTRL','sgCTRL'),('LOY','Y-depleted'),('PARENTAL','Parental')]]
handles += [Line2D([],[],color='#555555',marker=m,lw=0,label=t) for m,t in [('o','Aug 2025'),('^','May 2026')]]
fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.52,.995),ncol=5,fontsize=8,frameon=False)
fig.set_size_inches(11.8,9.8)
for oldax in fig.axes:
    pos=oldax.get_position()
    oldax.set_position([pos.x0,.29+.71*pos.y0,pos.width,.71*pos.height])
for txt in fig.texts:
    x,y=txt.get_position()
    txt.set_position((x,.29+.71*y))
# Keep the square period comparison at the full middle-row height.
b.set_position([.405,.35035,.26,.2769])
expr=pd.read_csv(FIG/'source_data/ht1080_rna_expression_changes/six_gene_sample_expression.csv')
genes=['UNC5A','ROR2','SDK2','COL1A2','FUCA2','MAGEA1']
for i,gene in enumerate(genes):
    ax=fig.add_axes([.065+i*.157,.085,.115,.175])
    sub=expr[expr.gene.eq(gene)]
    order=['PARENTAL','CTRL','LOY']
    boxes=ax.boxplot([sub.loc[sub.treatment.eq(t),'log2_CPM'] for t in order],positions=[0,1,2],widths=.56,patch_artist=True,showfliers=False,whis=1.5,medianprops=dict(color='#222222',lw=1))
    for box,t in zip(boxes['boxes'],order):
        box.set_facecolor(COL[t]);box.set_alpha(.25)
    for k,t in enumerate(order):
        for batch in ['AUG_2025','MAY_2026']:
            points=sub[sub.treatment.eq(t)&sub.batch.eq(batch)].sort_values('passage')
            offset=-.075 if batch=='AUG_2025' else .075
            jitter=np.linspace(-.035,.035,len(points))
            ax.scatter(k+offset+jitter,points.log2_CPM,s=22,c=COL[t],marker=SHAPE[batch],edgecolor='white',lw=.35,zorder=4)
    ax.set_xticks([0,1,2],['Parental','sgCTRL','Y-depleted'],rotation=45,ha='right',fontsize=7)
    ax.set_xlim(-.55,2.55)
    ax.set_ylabel('log2 CPM',fontsize=7)
    ax.tick_params(axis='y',labelsize=7)
    ax.set_title(gene,fontsize=9,fontstyle='italic',pad=5)
    ax.text(-.28,1.05,'ghijkl'[i],transform=ax.transAxes,fontsize=12,weight='semibold')
fig.text(.065,.012,'g–l: unadjusted all-15 TMM log2 CPM; parental n = 3 (August), sgCTRL and Y-depleted n = 6 each (August + May).',fontsize=7,color='#555555')

# f: compact PAR barplot beside d and e; gene boxes are g-l.
from ht1080_par_expression_panel import add_par_panel
par_data, par_summary = add_par_panel(fig, FIG, blue, red)
export(fig,'ch03_main_ht1080_rna_expression_changes',{'all_genes':d,'labelled_genes':selected,'panel_b_labelled_genes':agree,'six_gene_sample_expression':expr,'panel_f_PAR_expression':par_data,'panel_f_PAR_regional_summary':par_summary},'Combined F01: b shows 14 Y-linked genes plus the top three qualifying genes per direction across 15 samples, clustered on both axes by average linkage and Pearson correlation distance; former correlation b moved to S01g. Other panels retained, with d/e beside f, the PAR gene barplot. Gene expression boxes are g-l. Panel f shows saved within-period gene ratios, region-coloured bars and combined PAR1+PAR2 sums; solid/hatched bars denote August/May. Named genes require sgCTRL mean CPM >=5 in the displayed period; CSF2RA is omitted from individual bars. All annotated regional features, including CSF2RA, contribute to the combined sum. g-l use all-15 contextual CPM; f uses saved primary core-population CPM. No differential-expression analysis refitted. Reproduce with make_ht1080_rna_combined.py.')
