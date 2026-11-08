"""Condensed cross-line figures from saved results, 8 October 2026.

No model, enrichment test, correlation or P value is recomputed. Main panels
show Y inventory, X paralogs, saved effect correlations and seven named gene
sets. The former external-comparison panel becomes a separate supplement.
"""
from pathlib import Path
import argparse
import sys
import shutil
import hashlib
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

ROOT = Path(__file__).resolve().parents[3]
FIG = ROOT / 'figures/chapter3'
CT = ROOT.parent / 'RNAseq/analysis/crossline_SEPT_2026/outputs/tables'
DATA = FIG / 'source_data/crossline_concise'
DATA.mkdir(parents=True,exist_ok=True)
sys.path.insert(0,str(ROOT/'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
apply_style(8)
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--main-only', action='store_true', help='Update only the four-panel main figure.')
args = parser.parse_args()
LINES = ['HT1080','A549','DU145']
COLOURS = {'HT1080':'#00A087','A549':'#7E6148','DU145':'#8491B4'}
SIGNED = LinearSegmentedColormap.from_list('signed',['#3C5488','#F7F7F7','#E64B35'])
SIGNED.set_bad('white')
manifest=[]

def frozen(name):
    src,dst=CT/f'{name}.csv',DATA/f'{name}.csv'
    if not dst.exists():shutil.copy2(src,dst)
    assert src.read_bytes()==dst.read_bytes(),f'Source changed: {name}'
    manifest.append({'file':dst.relative_to(ROOT).as_posix(),'source':str(src),'sha256':hashlib.sha256(dst.read_bytes()).hexdigest()})
    return pd.read_csv(dst,low_memory=False)

inventory=frozen('Y_inventory_by_line')
regions=frozen('Y_PAR_region_sums_by_line')
paralogs=frozen('XY_paralogs_by_line')
correlations=frozen('gene_level_correlations')
agreement=frozen('gene_level_three_way_agreement')
shared=frozen('gene_sets_shared_divergent')
external=frozen('external_replication_sets')
registry=frozen('external_dataset_registry_used')
FOCUS=[
    ('HALLMARK_INTERFERON_ALPHA_RESPONSE','Interferon-α response'),
    ('HALLMARK_INTERFERON_GAMMA_RESPONSE','Interferon-γ response'),
    ('HALLMARK_INFLAMMATORY_RESPONSE','Inflammatory response'),
    ('REACTOME_RRNA_PROCESSING','rRNA processing'),
    ('REACTOME_M_PHASE','M phase'),
    ('HALLMARK_G2M_CHECKPOINT','G2/M checkpoint'),
    ('HALLMARK_E2F_TARGETS','E2F targets'),
]
all_source=CT/'gene_sets_all_common.csv'
all_sets=pd.read_csv(all_source,low_memory=False)
focus=all_sets.set_index('pathway').loc[[key for key,label in FOCUS]].reset_index()
focus.to_csv(DATA/'panel_d_selected_gene_sets.csv',index=False)
manifest.append({'source':str(all_source),'sha256':hashlib.sha256(all_source.read_bytes()).hexdigest(),'display_rows':[key for key,label in FOCUS],'display_table':'panel_d_selected_gene_sets.csv'})
primary=all_sets[(all_sets.collection=='H')|((all_sets.collection=='C2')&(all_sets.subcollection=='CP:REACTOME'))|((all_sets.collection=='C5')&(all_sets.subcollection=='GO:BP'))]
assert len(primary)==4634 and sum(primary['class']=='shared')==4
assert shared['class'].value_counts().to_dict()=={'shared':46,'divergent':23}
assert len(correlations)==3 and set(correlations.genes)=={16165}

def save(fig,stem):
    for ext in ['pdf','svg','png']:
        stage=FIG/f'{stem}.staging.{ext}';target=FIG/f'{stem}.{ext}'
        fig.savefig(stage,dpi=400,facecolor='white');stage.replace(target)
    fig.savefig(FIG/'previews'/f'{stem}_preview.png',dpi=160)
    plt.close(fig)
    print('written',stem)

def letter(fig,x,y,text):
    fig.text(x,y,text,fontsize=12,weight='semibold',va='bottom')

fig=plt.figure(figsize=(9.5,5.8))
a=fig.add_axes([.075,.59,.49,.34])
b=fig.add_axes([.77,.59,.205,.34])
c=fig.add_axes([.13,.09,.25,.29])
d=fig.add_axes([.755,.085,.185,.305])
genes=['RPS4Y1','ZFY','PRKY','USP9Y','DDX3Y','UTY','TXLNGY','KDM5D','EIF1AY']
panel_a_rows=[]
bar_floor=.07
for k,line in enumerate(LINES):
    vals=inventory[inventory.line==line].set_index('gene').reindex(genes).percent_of_control
    a.bar(np.arange(len(genes))+(k-1)*.24,vals.clip(lower=.1)-bar_floor,
          width=.22,bottom=bar_floor,color=COLOURS[line],label=line,zorder=3,linewidth=0)
    # Label the underlying percentage, including values below the plotting floor.
    for x,value in zip(np.arange(len(genes))+(k-1)*.24,vals):
        label=f'{value:.2f}' if 0 < value < 1 else f'{value:.1f}'
        if value == 0: label='0'
        a.annotate(label,(x,max(value,.1)),xytext=(0,2),textcoords='offset points',
                   ha='center',va='bottom',fontsize=5.5,color='#444444',alpha=.72,
                   annotation_clip=False,gid='panel-a-percentage')
    panel_a_rows.extend(dict(line=line,gene=gene,percent_of_control=value,
                             displayed_percent=max(value,.1),bar_base_percent=bar_floor)
                        for gene,value in vals.items())
a.axhline(100,ls='--',lw=.7,color='#999999')
a.axvspan(-.5,1.5,color='#00A087',alpha=.06,lw=0)
a.text(.5,145,'distal Yp',ha='center',fontsize=8,color='#555555')
a.set_yscale('log');a.set_ylim(.07,300)
a.set_yticks([.1,1,10,100],['≤0.1','1','10','100'])
a.set_xticks(range(len(genes)),genes,rotation=45,ha='right',fontsize=8)
a.set(xlim=(-.6,len(genes)-.4),ylabel='Y-depleted, % of sgCTRL')
a.legend(fontsize=8,loc='upper right',ncol=3,frameon=False,handlelength=1.2,handletextpad=.4,columnspacing=.7,borderaxespad=.1)
par=regions[regions.region=='PAR1'].set_index('line').percent_of_control
a.text(8.45,37,'PAR1: '+ '   '.join(f'{line} {par[line]:.0f}%' for line in LINES),ha='right',va='top',fontsize=7.5,color='#555555')
ok=paralogs.groupby('X').X_CTRL_CPM.min()>=1
pairs=[gene for gene in ['DDX3X','EIF1AX','KDM5C','KDM6A','USP9X','ZFX','RPS4X','PRKX','TMSB4X'] if ok.get(gene,False)]
ys=np.arange(len(pairs))[::-1]
panel_b_rows=[]
for k,line in enumerate(LINES):
    data=paralogs[paralogs.line==line].set_index('X')
    vals=data.reindex(pairs).X_mean_difference
    b.barh(ys-(k-1)*.24,vals,height=.20,color=COLOURS[line],zorder=3,linewidth=0)
    b.axvline(data.chrX_median_shift.iloc[0],color=COLOURS[line],lw=.9,ls=':',zorder=4)
    panel_b_rows.extend(dict(line=line,gene=gene,mean_log2_CPM_difference=value,
                             chrX_median_shift=data.chrX_median_shift.iloc[0])
                        for gene,value in vals.items())
b.axvline(0,color='#999999',lw=.6)
b.set_yticks(ys,pairs,fontsize=8);b.set_xlim(-.45,.45)
b.set_xlabel('X-paralog expression difference\n(mean log₂ CPM)',fontsize=8)

matrix=np.full((3,3),np.nan)
for row in correlations.itertuples():
    i,j=LINES.index(row.line_1),LINES.index(row.line_2)
    matrix[i,j]=matrix[j,i]=row.Spearman
im=c.imshow(matrix,cmap=SIGNED,vmin=-1,vmax=1,aspect='auto')
for i in range(3):
    for j in range(3):
        c.text(j,i,'—' if i==j else f'{matrix[i,j]:+.2f}',ha='center',va='center',fontsize=10)
c.set_xticks(range(3),LINES,fontsize=8);c.xaxis.tick_top()
c.set_yticks(range(3),LINES,fontsize=8)
c.tick_params(length=0)
c.set_xlabel('Expression-effect correlation\nSpearman ρ; 16,165 genes',fontsize=8,labelpad=7)
for spine in c.spines.values():spine.set_visible(False)
cax=fig.add_axes([.397,.09,.011,.29])
cb=fig.colorbar(im,cax=cax,ticks=[-1,0,1]);cb.ax.tick_params(labelsize=8)

nes=focus[[f'NES_{line}' for line in LINES]].to_numpy()
bh=focus[[f'BH_{line}' for line in LINES]].to_numpy()
im2=d.imshow(nes,cmap=SIGNED,vmin=-2.6,vmax=2.6,aspect='auto')
for i in range(7):
    for j in range(3):
        d.text(j,i,f'{nes[i,j]:+.1f}'+('*' if bh[i,j]<.05 else ''),ha='center',va='center',fontsize=8)
d.set_xticks(range(3),LINES,fontsize=8);d.xaxis.tick_top()
d.set_yticks(range(7),[label for key,label in FOCUS],fontsize=8)
d.tick_params(length=0)
d.axhline(3.5,color='white',lw=2)
for spine in d.spines.values():spine.set_visible(False)
d.set_xlabel('* FDR <0.05',fontsize=8,labelpad=7)
cax2=fig.add_axes([.954,.085,.011,.305])
cb2=fig.colorbar(im2,cax=cax2,ticks=[-2,0,2]);cb2.ax.set_title('NES',fontsize=8,pad=7)
letter(fig,.015,.955,'a');letter(fig,.66,.955,'b')
letter(fig,.015,.425,'c');letter(fig,.47,.425,'d')
# Stagger neighbouring percentage labels only where their text would overlap.
fig.canvas.draw()
renderer=fig.canvas.get_renderer()
placed=[]
for label in sorted((t for t in a.texts if t.get_gid()=='panel-a-percentage'),
                    key=lambda t: (t.xy[1],t.xy[0])):
    bounds=label.get_window_extent(renderer).expanded(1.08,1.15)
    while any(bounds.overlaps(other) for other in placed):
        label.set_position((0,label.get_position()[1]+1))
        bounds=label.get_window_extent(renderer).expanded(1.08,1.15)
    placed.append(bounds)
assert len(placed)==27
save(fig,'ch03_main_crossline_rna_comparison')
pd.DataFrame(panel_a_rows).to_csv(DATA/'panel_a_bar_display.csv',index=False)
pd.DataFrame(panel_b_rows).to_csv(DATA/'panel_b_bar_display.csv',index=False)
(DATA/'bar_display_notes.json').write_text(json.dumps({
    'date':'2026-10-10',
    'change':'Panel a bars annotated with underlying percentages in small, semi-transparent grey text.',
    'panel_a':{'scale':'logarithmic','bar_base_percent':bar_floor,'lower_display_limit_percent':.1,
               'interpretation':'Bar tops show residual expression; bar lengths are not linear percentages.'},
    'panel_b':{'scale':'linear','bar_base':0,'order_top_to_bottom':LINES},
    'unchanged':'Source values, categories, colours, reference lines, panel c and panel d.'
},indent=2)+'\n',encoding='utf-8')
if args.main_only:
    print('Updated main figure only; external supplement and gene-set tables unchanged.')
    raise SystemExit(0)

# Preserve all eight selected sets and all 21 rows from the former main panel e.
cols=['HALLMARK_INTERFERON_ALPHA_RESPONSE','HALLMARK_INTERFERON_GAMMA_RESPONSE','HALLMARK_INFLAMMATORY_RESPONSE','ONDER_CDH1_TARGETS_2_DN','BMI1_DN.V1_UP','MEL18_DN.V1_UP','HALLMARK_G2M_CHECKPOINT','HALLMARK_E2F_TARGETS']
col_labels=['IFN-α','IFN-γ','Inflammatory','CDH1 targets','BMI1_DN','MEL18_DN','G2/M','E2F']
labels={'HT1080':'HT1080','A549':'A549','DU145':'DU145','A549_clones_LOYvROY':'A549 LOY clones','A549_multiome_pb_LOYvROY':'A549 clones, multiome*','GSE313572_HT1080_ctrl_YnegvYpos':'HT1080 spheroids','Tahoe_A549_DMSO_YnegvYpos':'A549 (Tahoe)','Tahoe_J82_DMSO_YnegvYpos':'J82 (Tahoe)','Tahoe_SKMEL2_DMSO_YnegvYpos':'SK-MEL-2 (Tahoe)','Tahoe_HCT15_DMSO_YnegvYpos':'HCT15 (Tahoe, 6 wells)','MB49_cellline_YnegvYpos':'MB49 cells (mouse)','MB49_tumor_YnegvYpos':'MB49 tumours (mouse)','GSE126712_lymph_XOvXY':'Lymphocytes','GSE207112_troph_XOvXY':'iPSC trophoblast','GSE264742_NCC_XOvXY':'iPSC neural crest','ERP156313_fetal_XOvXY':'Fetal tissues','PageLab_Fib_X0vXY':'Fibroblasts','PageLab_LCL_X0vXY':'Lymphoblastoid lines','GSE301297_undiff_UTYKOvWT':'hESC UTY-KO','GSE301297_undiff_DKOvWT':'hESC UTX/UTY-DKO','XYY_integrated_flipped':'XYY vs XY (reversed)'}
group_labels={'in_house':'Thesis models','A_somatic_LOY':'Somatic LOY','A2_single_cell_mosaic':'Natural Y mosaicism\n(single-cell)','B_45X_vs_46XY':'45,X vs 46,XY','C_not_LOY':'Related Y-dosage\ncomparisons'}
rows=list(registry.dataset)
A=external.pivot(index='dataset',columns='gene_set',values='AUC').reindex(index=rows,columns=cols)
P=external.pivot(index='dataset',columns='gene_set',values='empirical_P').reindex(index=rows,columns=cols)
assert A.shape==(21,8) and not A.isna().any().any()
A.to_csv(DATA/'external_AUC_display.csv');P.to_csv(DATA/'external_P_display.csv')
fig=plt.figure(figsize=(9.3,5.8))
ax=fig.add_axes([.245,.055,.535,.835])
im=ax.imshow(A.to_numpy()-.5,cmap=SIGNED,vmin=-.3,vmax=.3,aspect='auto')
for i in range(len(rows)):
    for j in range(len(cols)):
        if P.iat[i,j]<.05:ax.text(j,i,'•',ha='center',va='center',fontsize=9,color='#222222')
ax.set_xticks(range(len(cols)),col_labels,fontsize=8,rotation=30,ha='left');ax.xaxis.tick_top()
ax.set_yticks(range(len(rows)),[labels[r] for r in rows],fontsize=8)
ax.tick_params(length=0)
for spine in ax.spines.values():spine.set_visible(False)
groups=dict(zip(registry.dataset,registry.group))
bounds=[i for i in range(1,len(rows)) if groups[rows[i]]!=groups[rows[i-1]]]
for boundary in bounds:ax.axhline(boundary-.5,color='white',lw=2)
for start,end in zip([0]+bounds,bounds+[len(rows)]):
    ax.text(7.65,(start+end-1)/2,group_labels[groups[rows[start]]],ha='left',va='center',fontsize=8,color='#555555')
ax.axvline(5.5,color='white',lw=2)
cax=fig.add_axes([.942,.055,.014,.835])
cb=fig.colorbar(im,cax=cax,ticks=[-.3,0,.3]);cb.ax.set_title('AUC − 0.5',fontsize=7,pad=8)
save(fig,'ch03_supp_crossline_external_comparison')

table_dir=ROOT/'tables/chapter3/crossline_comparison'
table_dir.mkdir(parents=True,exist_ok=True)
shared.to_csv(table_dir/'shared_and_divergent_gene_sets.tsv',sep='\t',index=False)
primary[['pathway','collection','subcollection']+[name for line in LINES for name in [f'NES_{line}',f'BH_{line}']]+['class']].to_csv(table_dir/'primary_collection_results.tsv',sep='\t',index=False)
(DATA/'input_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
print('Saved complete gene-set tables and frozen-input provenance.')
