"""Review-only six-panel DU145 draft with the original GRHL2/TAF9/PXDN plots; no refits.

Run from the thesis root with .tmp/qi_figure_python on PYTHONPATH.
Does not alter selected figures, captions, numbering or manuscript sources.
"""
from pathlib import Path
import hashlib,json,shutil,sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib import font_manager
import numpy as np
import pandas as pd

FIG=Path(__file__).resolve().parents[1]
ROOT=FIG.parents[1]
sys.path.insert(0,str(ROOT/'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
font=apply_style(8)
BLUE,RED,DARK='#3C5488','#E64B35','#333333'
COND={'CTRL':BLUE,'LOY':RED}
MARKER={'P12':'o','P13':'s','P14':'^'}
OUT=FIG/'review/du145_global_pathways_20261030_v2_genes'
SRC=FIG/'source_data/du145_global_pathways_genes_draft_20261030'
OUT.mkdir(parents=True,exist_ok=True)
(SRC/'inputs').mkdir(parents=True,exist_ok=True)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
paths={
 'pca':FIG/'source_data/du145_global_expression_combined/main/panel_a_pca.csv',
 'variance':FIG/'source_data/du145_global_expression_combined/main/panel_a_variance.csv',
 'associations':FIG/'source_data/du145_global_expression_combined/main/panel_a_design_associations.csv',
 'genes':FIG/'source_data/du145_global_expression_combined/main/panel_b_c_gene_level.csv',
 'pathways':FIG/'source_data/du145_rna_programmes_regulators/panel_a_b_pathways.csv',
 'pathway_counts':FIG/'source_data/du145_rna_programmes_regulators/panel_c_theme_status_counts.csv',
 'grhl2':FIG/'source_data/du145_global_expression_combined/main/panel_f_GRHL2_CPM.csv',
 'taf9_pxdn':FIG/'source_data/du145_global_expression_combined/main/panel_g_h_TAF9_PXDN_logCPM.csv',
}
preserve=[ROOT/'draftv2/chapter3/chapter3_results.md',ROOT/'Chapter3_figure_captions.md',ROOT/'draftv2/chapter3/chapter3_figure_numbering.md']
preserve.extend((FIG/'review/du145_global_pathways_20261030_v1').glob('*'))
preserve.append(FIG/'code/draft_du145_global_pathways_20261030.py')
for folder in ['CH3_S3_4_6_3_F01_DU145_RNA_global','CH3_S3_4_6_3_F04_DU145_RNA_programmes']:
    preserve.extend((FIG/'selected'/folder).glob('*'))
protected={str(p.relative_to(ROOT)):sha(p) for p in preserve if p.is_file()}
input_hashes={str(p.relative_to(ROOT)):sha(p) for p in paths.values()}
tables={}
for key,p in paths.items():
    dest=SRC/'inputs'/p.name
    if not dest.exists(): shutil.copy2(p,dest)
    assert sha(dest)==sha(p),('Frozen input differs',p)
    tables[key]=pd.read_csv(dest)
pca,var,assoc,de,display,counts=(tables[k] for k in ['pca','variance','associations','genes','pathways','pathway_counts'])
assert len(pca)==6 and len(de)==18767 and len(display)==19
sig=de.FDR<.05
subset=(de.FDR<.05)&(de.logFC.abs()>=.5)&de.pairs_agreeing_with_logFC.eq(3)
assert np.array_equal(subset.to_numpy(),de.reporting_subset.to_numpy())
assert sig.sum()==147 and subset.sum()==114
cn_stable=(subset&de.CN_class.eq('CN_stable')).sum()
assert cn_stable==88 and counts.retained.sum()==139 and counts.sensitive.sum()==15
assert display.pathway.is_unique and np.isfinite(display[['NES','fgsea_FDR','CAMERA_FDR']]).all().all()
assert (display.loc[display.fgsea_FDR<.05,'survives_CN_screen']).all()

# All original pathway rows, their ordering, and both FDR columns are retained.
display.to_csv(SRC/'panel_f_pathways.csv',index=False)
pd.testing.assert_frame_equal(pd.read_csv(SRC/'panel_f_pathways.csv'),display,check_exact=False,rtol=1e-14,atol=1e-25)

fig=plt.figure(figsize=(7.2,9.0))
fig.subplots_adjust(0,0,1,1)
a=fig.add_axes([.105,.760,.32,.205])
b=fig.add_axes([.555,.760,.415,.205])
variance=var.set_index('PC').core_variance
associations=assoc.set_index('axis')
for passage,grp in pca.groupby('passage_block'):
    grp=grp.set_index('treatment')
    a.plot(grp.loc[['CTRL','LOY'],'PC1'],grp.loc[['CTRL','LOY'],'PC2'],color='#AAAAAA',ls='--',lw=.7,zorder=1)
for row in pca.itertuples():
    a.scatter(row.PC1,row.PC2,s=37,color=COND[row.treatment],marker=MARKER[row.passage_block],edgecolors='white',linewidths=.35,zorder=3)
a.set(xlim=(-48,80),ylim=(-57,58),xlabel=f'PC1 ({100*variance[1]:.1f}%)',ylabel=f'PC2 ({100*variance[2]:.1f}%)')
a.set_xticks([-40,0,40,80]);a.set_yticks([-40,0,40])
a.text(.025,.97,f"PC1 ~ passage: R² = {associations.loc['PC1','passage_R2']:.2f}\nPC2 ~ condition: R² = {associations.loc['PC2','treatment_R2']:.2f}",transform=a.transAxes,va='top',fontsize=7,color=DARK)
handles=[Line2D([],[],ls='',marker='o',color=BLUE,label='sgCTRL',ms=4),Line2D([],[],ls='',marker='o',color=RED,label='Y-depleted',ms=4)]
handles += [Line2D([],[],ls='',marker=m,color='#666666',label=p,ms=3.7) for p,m in MARKER.items()]
a.legend(handles=handles,loc='upper left',bbox_to_anchor=(-.04,-.255),ncol=5,fontsize=6.5,handlelength=.6,columnspacing=.65,handletextpad=.25,borderaxespad=0)
fig.text(.025,.978,'a',fontsize=13,weight='semibold')

cn_concordant=de.CN_class.eq('CN_changed_concordant')
categories=[(~sig,'#D5D5D5',2.4,'Not significant'),(sig&~subset&~cn_concordant,'#F0B27A',12,'Other FDR < 0.05'),(sig&~subset&cn_concordant,BLUE,12,'CN-concordant'),(subset,RED,12,'Consistent subset')]
for mask,colour,size,label in categories:
    b.scatter(de.loc[mask,'logFC'],-np.log10(de.loc[mask,'PValue']),s=size,color=colour,lw=0,rasterized=True,zorder=2 if label=='Not significant' else 3)
for x in [-.5,.5]: b.axvline(x,lw=.6,ls='--',color='#BBBBBB',zorder=1)
b.set(xlim=(-5.4,3.4),ylim=(0,14.5),xlabel='log₂ fold change (Y-depleted / sgCTRL)',ylabel='−log₁₀ P')
b.set_xticks([-4,-2,0,2]);b.set_yticks([0,5,10])
b.text(.985,.975,'147 genes at FDR < 0.05\n114 in consistent subset\n88/114 (77%) CN-stable',transform=b.transAxes,ha='right',va='top',fontsize=7,color=DARK)
label_positions={'TACSTD2':(-2.85,10.1),'AP1M2':(-2.85,9.0),'CLDN7':(-2.85,7.9),'GRHL2':(-5.05,6.55),'TAF9':(1.7,9.4),'COL4A2':(1.7,8.3),'CREB3L1':(1.7,7.2),'F13A1':(1.7,6.1)}
label_table=de[de.gene_name.isin(label_positions)].copy()
assert len(label_table)==len(label_positions) and label_table.reporting_subset.all()
for row in label_table.itertuples():
    b.annotate(row.gene_name,(row.logFC,-np.log10(row.PValue)),xytext=label_positions[row.gene_name],fontsize=6.8,ha='left',va='center',color=DARK,bbox=dict(boxstyle='square,pad=.04',fc='white',ec='none',alpha=.9),arrowprops=dict(arrowstyle='-',color='#999999',lw=.45,shrinkA=2,shrinkB=2))
label_table.to_csv(SRC/'panel_b_labelled_genes.csv',index=False)
b.legend(handles=[Line2D([],[],ls='',marker='o',color=colour,label=label,ms=3.5) for _,colour,_,label in categories[1:]],loc='upper left',bbox_to_anchor=(-.04,-.255),ncol=3,fontsize=6.3,handlelength=.6,handletextpad=.3,columnspacing=.55,borderaxespad=0)
fig.text(.485,.978,'b',fontsize=13,weight='semibold')

# Middle row: original Fig 3.12f–h, using the saved sample values and model results.
expression_records=[]
expression_artist_values=[]
expression_stat_artists=[]
for letter,gene,xpos in [('c','GRHL2',.105),('d','TAF9',.425),('e','PXDN',.745)]:
    ax=fig.add_axes([xpos,.504,.225,.141])
    result=de.loc[de.gene_name.eq(gene)].iloc[0]
    if gene=='GRHL2':
        values=tables['grhl2'].copy()
        values['gene_name']=gene
        values['passage']=values.sample_key.str.extract(r'(P1[234])$')[0]
        measure='CPM'
    else:
        values=tables['taf9_pxdn'].loc[tables['taf9_pxdn'].gene_name.eq(gene)].copy()
        measure='log2_CPM'
    assert len(values)==6 and values.gene_id.eq(result.gene_id).all()
    for i,treatment in enumerate(['CTRL','LOY']):
        group=values.loc[values.treatment.eq(treatment)].sort_values('passage')
        assert group.passage.tolist()==['P12','P13','P14']
        colour=COND[treatment]
        if gene=='GRHL2':
            ax.plot([i-.22,i+.22],[group.CPM.mean()]*2,color=colour,lw=1.3)
        else:
            ax.boxplot([group[measure]],positions=[i],widths=.55,patch_artist=True,
                       showfliers=False,whis=1.5,
                       boxprops=dict(facecolor=colour,alpha=.18,edgecolor=colour,linewidth=.8),
                       medianprops=dict(color=colour,linewidth=1.2),
                       whiskerprops=dict(color=colour,linewidth=.8),capprops=dict(color=colour,linewidth=.8))
        dots=ax.scatter(i+np.array([-.09,0,.09]),group[measure],s=22,color=colour,edgecolors='white',linewidths=.3,zorder=3)
        assert np.array_equal(np.asarray(dots.get_offsets())[:,1],group[measure].to_numpy())
        expression_artist_values.extend(np.asarray(dots.get_offsets())[:,1].tolist())
        for row in group.itertuples():
            expression_records.append({'panel':letter,'gene':gene,'sample_key':row.sample_key,'passage':row.passage,'treatment':treatment,'measure':measure,'value':getattr(row,measure),'model_log2FC':result.logFC,'model_FDR':result.FDR})
    ax.set_xlim(-.55,1.55)
    if gene=='GRHL2':
        ax.set_ylim(0,4.6);ax.set_yticks([0,1,2,3,4])
    else:
        span=max(values[measure].max()-values[measure].min(),.5)
        ax.set_ylim(values[measure].min()-.15*span,values[measure].max()+.5*span)
        ax.set_yticks([2,3,4] if gene=='TAF9' else [-2,-1,0,1,2])
    ax.set_xticks([0,1],['sgCTRL','Y-depleted'],fontsize=7)
    ax.tick_params(axis='y',labelsize=7)
    ax.set_ylabel(gene+(' (CPM)' if gene=='GRHL2' else ' (log₂ CPM)'),fontsize=8)
    stat_text=f'FDR = {result.FDR:.1e}\nlog₂FC = {result.logFC:.2f}'
    stat_artist=ax.text(.5,.975,stat_text,transform=ax.transAxes,ha='center',va='top',fontsize=7,color=DARK)
    expression_stat_artists.append((stat_artist,stat_text))
    fig.text(xpos-.07,.660,letter,fontsize=13,weight='semibold')
assert len(expression_records)==18 and len(expression_artist_values)==18
pd.DataFrame(expression_records).to_csv(SRC/'panels_c_d_e_expression.csv',index=False)

# Bottom panel: NES bars replace the old colour strip; FDRs remain row-aligned.
c=fig.add_axes([.505,.055,.195,.335])
fd=fig.add_axes([.745,.055,.225,.335])
nrows=int(display.row.max())+1
c.set(ylim=(nrows-.5,-.5),xlim=(-2.6,2.6),xlabel='Normalised enrichment score (NES)')
c.set_xticks([-2,-1,0,1,2]);c.set_yticks(display.row,display.display_label,fontsize=7.1)
c.tick_params(axis='y',length=0,pad=4)
c.tick_params(axis='x',labelsize=7,pad=2)
c.xaxis.label.set_fontsize(7)
c.axvline(0,color='#A0A0A0',lw=.7,zorder=1)
c.spines['left'].set_visible(False)
for row in display.itertuples():
    colour=BLUE if row.NES<0 else RED
    c.barh(row.row,row.NES,height=.7,color=colour if row.fgsea_FDR<.05 else 'white',edgecolor=colour,lw=.7,zorder=2)
for theme,group in display.groupby('theme',sort=False):
    c.text(-2.77,group.row.min()-1,theme,fontsize=7.2,weight='semibold',ha='right',va='center',color=DARK,clip_on=False)
fd.set(xlim=(-.5,1.5),ylim=(nrows-.5,-.5));fd.axis('off')
fdr_artists=[]
for row in display.itertuples():
    for col,value in enumerate([row.fgsea_FDR,row.CAMERA_FDR]):
        artist=fd.text(col,row.row,f'{value:.1e}',fontsize=7.3,ha='center',va='center',weight='semibold' if value<.05 else 'normal',color=DARK)
        fdr_artists.append((artist,f'{value:.1e}',row.pathway,col))
for col,label in enumerate(['fgsea','CAMERA']): fd.text(col,-1.45,label,fontsize=8,ha='center',va='bottom')
fig.text(.025,.425,'f',fontsize=13,weight='semibold')
fig.text(.079,.425,'23,252 sets tested; 154 significant; 139 retain significance after CN exclusion',fontsize=7.5)
fig.text(.855,.425,'Global FDR',ha='center',fontsize=8,weight='semibold')
# Explanatory footer text belongs in the separate caption, not on the canvas.

# Check the values actually attached to text artists, plus page-edge clipping.
fig.canvas.draw()
renderer=fig.canvas.get_renderer()
for artist,expected,pathway,col in fdr_artists:
    assert artist.get_text()==expected,(pathway,col)
    bbox=artist.get_window_extent(renderer)
    assert bbox.x0>0 and bbox.x1<fig.bbox.width and bbox.y0>0 and bbox.y1<fig.bbox.height
for artist,expected in expression_stat_artists:
    assert artist.get_text()==expected
for p,digest in input_hashes.items(): assert sha(ROOT/p)==digest,p
for p,digest in protected.items(): assert sha(ROOT/p)==digest,p

stem='ch03_main_du145_global_pathways_genes_draft'
for ext in ['pdf','svg','png']:
    fig.savefig(OUT/f'{stem}.{ext}',dpi=600)
fig.savefig(OUT/f'{stem}_preview.png',dpi=220)
digest=sha(OUT/f'{stem}_preview.png')[:12]
preview=OUT/f'{stem}_preview_{digest}.png'
if not preview.exists(): shutil.copy2(OUT/f'{stem}_preview.png',preview)
plt.close(fig)

caption='''# Draft main figure — combined DU145 expression, examples and pathways

**Transcriptional differences and coordinated programmes in Y-depleted DU145 cells.** Three matched passages of one sgCTRL and one Y-depleted population. **a,** PCA of 18,767 autosomal/non-PAR X genes; shapes denote passages and dashed lines connect matched passages. **b,** Passage-blocked voom–limma results (residual df = 2; BH correction across 18,767 genes). Red marks 114 genes with FDR <0.05, |log₂ fold change| ≥0.5 and consistent direction across passages; blue/orange distinguish other significant genes with/without copy-number concordance. Of these 114 genes, 88 have |lpWGS difference| <0.2. **c–e,** GRHL2 CPM with group means (c), and TAF9/PXDN TMM-normalised log₂ CPM, prior count 2 (d,e). Points show passage samples; boxes show medians and interquartile ranges, with whiskers to observations within 1.5 interquartile ranges. FDRs and fold changes come from b. **f,** The same 19 selected gene sets as current Fig 3.13a,b, with saved fgsea NES and fgsea/CAMERA FDRs. Passage-blocked CAMERA used fixed correlation 0.01; each method retains BH correction across 23,252 sets. Filled bars and bold FDRs mark significance for fgsea and the corresponding method. Of all 154 significant fgsea sets, 139 retain significance and direction after excluding genes with |lpWGS difference| ≥0.2. CGP up/down suffixes describe reference signatures. Enrichment measures coordinated expression, not pathway activity.

CN, copy number; H, Hallmark; R, Reactome; CGP, perturbation signatures; MF/CC, GO molecular function/cellular component.

This six-panel review draft adds original Fig 3.12f–h as c–e and moves the pathway summary to f. It retains all 19 pathway rows, both FDR columns and all 18 expression measurements. The three-panel draft, selected figures and manuscript remain unchanged. The volcano retains every tested gene, with eight direct labels for readability. Explanatory caption text remains separate from the figure. The former passage-consistency, copy-number and pathway-theme count panels remain in their original figures pending a decision on supplementary placement.
'''
(OUT/'caption.md').write_text(caption,encoding='utf-8')
validation={'status':'review draft; not adopted','dimensions_inches':[7.2,9.0],'font':font,'font_file':font_manager.findfont(font),'pca_samples':len(pca),'genes_plotted':len(de),'significant_genes':int(sig.sum()),'consistent_subset_genes':int(subset.sum()),'copy_number_stable_subset_genes':int(cn_stable),'pathway_rows_retained':len(display),'fdr_values_retained':len(fdr_artists),'all_displayed_fdr_strings_match_saved_values':True,'pathway_bh_family_each_method':23252,'input_hashes':input_hashes,'protected_original_hashes':protected,'all_protected_originals_unchanged':True,'model_refits':False,'preview':str(preview),'outputs':{p.name:sha(p) for p in OUT.glob(stem+'.*')},'scope':'Draft PCA + volcano + all 19 pathway NES/FDR rows; no manuscript or selected-figure adoption.'}
validation.update({'expression_genes':['GRHL2','TAF9','PXDN'],'expression_measurements_plotted':len(expression_records),'expression_values_verified_against_scatter_artists':True,'expression_statistics_from_existing_model':True,'prior_three_panel_draft_preserved':True,'scope':'Six-panel review draft: PCA, volcano, GRHL2/TAF9/PXDN expression and all 19 pathway NES/FDR rows. No manuscript or selected-figure adoption.'})
(SRC/'validation.json').write_text(json.dumps(validation,indent=2),encoding='utf-8')
(SRC/'README.md').write_text('# DU145 combined figure draft — 30 October 2026\n\nReview-only combination of current draftv2 Figs 3.12 and 3.13. Inputs are frozen copies of the previously rendered tables. All PCA coordinates, gene-level effects/P values/FDRs, pathway NES, pathway order and both pathway FDR columns are unchanged. Statistical testing families were not recalculated.\n\nReproduce with `figures/chapter3/code/draft_du145_global_pathways_genes_20261030.py`, using bundled Python and `.tmp/qi_figure_python` on PYTHONPATH. Exports and caption: `figures/chapter3/review/du145_global_pathways_20261030_v2_genes`. No main/supplementary reassignment has been applied. The validation record includes original input/selected-figure/manuscript hashes. Original sample-key prefixes are historical identifiers; the existing DU145 sample mapping is retained.\n',encoding='utf-8')
print(json.dumps({k:validation[k] for k in ['status','dimensions_inches','genes_plotted','pathway_rows_retained','fdr_values_retained','preview']},indent=2),flush=True)
with (SRC/'README.md').open('a',encoding='utf-8') as f:
    f.write('\nThis version adds original Fig 3.12f–h as c–e and places the pathway panel at f. GRHL2 retains the original linear CPM measurements and mean lines; TAF9 and PXDN retain log2 CPM values (prior count 2), boxes and individual passage measurements. All 18 measurements, three saved gene FDRs and three model log2 fold changes are retained. The figure is 7.2 × 9.0 inches; fonts keep their original physical sizes. The previous three-panel draft is protected and unchanged. No explanatory caption text is placed on the figure.\n')
