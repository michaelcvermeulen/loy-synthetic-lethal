"""Display-only split of Fig 4.3 using frozen analysis outputs."""
import hashlib
import json
import make_ht1080_candidates as C
from ht1080_saga_subunits import draw_saga_membership
from ht1080_method_comparison_panels import draw_comparisons, hashes as method_hashes
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.text import Text

MAIN=C.SUPP
SUPPLEMENT='ch04_supp_ht1080_method_and_d0_support'
DATA=C.OUT/'source_data/ht1080_support_split_20261028'
DATA.mkdir(exist_ok=True)
GAMETOLOG_ROWS=[['DDX3X','EIF1AX','KDM5C'],['ZFX','RPS4X','KDM6A']]
D0_JACKPOTTED_GUIDES={('DDX3X',2),('KDM6A',4)}
D0_GUIDE_LABEL='Brunello sgRNA (1–4 per gene)'
D0_FOOTNOTE='† D0 PCR-jackpotted guide (masked in IsogenicZ)   * Zero reads in all six samples'


def input_hashes():
    paths=[C.DATA/'guide_profiles.csv',C.DATA/'all_gene_scores.csv',C.DATA/'saga_reference_scores.csv']
    return {**{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},**method_hashes()}


def check_text(fig, stem):
    fig.canvas.draw(); renderer=fig.canvas.get_renderer()
    texts=list(fig.texts)
    for ax in fig.axes:
        texts.extend(ax.texts)
        if ax.axison:
            texts.extend([ax.title,ax.xaxis.label,ax.yaxis.label])
            texts.extend(ax.get_xticklabels()+ax.get_yticklabels())
    for legend in fig.legends:texts.extend(legend.get_texts())
    texts=[t for t in texts if t.get_visible() and t.get_text().strip()]
    # Annotation.get_window_extent includes its leader line; check text alone.
    entries=[(t,Text.get_window_extent(t,renderer=renderer)) for t in texts]
    overlaps=[]
    for i,(a,ab) in enumerate(entries):
        for b,bb in entries[i+1:]:
            if min(ab.x1,bb.x1)-max(ab.x0,bb.x0)>.4 and min(ab.y1,bb.y1)-max(ab.y0,bb.y0)>.4:
                overlaps.append([a.get_text(),b.get_text()])
    outside=[t.get_text() for t,b in entries if b.x0<-.5 or b.y0<-.5 or b.x1>fig.bbox.width+.5 or b.y1>fig.bbox.height+.5]
    report={'displayed_text_items':len(entries),'overlaps':overlaps,'outside_page':outside}
    (DATA/f'{stem}_text_layout.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    assert not outside,(stem,'outside',outside)
    assert not overlaps,(stem,'overlap',overlaps)
    return report


def baseline_heatmap(fig,genes,bounds,panel):
    data=C.PROFILES.loc[C.PROFILES.timepoint.eq('D0') & C.PROFILES.gene.isin(genes)].copy()
    lookup=data.set_index(['gene','genotype','guide_number'])
    matrix=np.array([[lookup.loc[(gene,genotype,guide),'log2_cpm_plus_1']
                     for genotype in ['CTRL','LOY'] for guide in range(1,5)] for gene in genes])
    assert matrix.shape==(len(genes),8) and np.isfinite(matrix).all()
    ax=fig.add_axes(bounds)
    mesh=ax.pcolormesh(np.arange(9),np.arange(len(genes)+1),matrix,
                      cmap='magma',vmin=0,vmax=15,edgecolors='white',linewidth=.55)
    ax.set_xlim(0,8);ax.set_ylim(len(genes),0)
    ax.set_yticks(np.arange(len(genes))+.5,['KDM6A (UTX)' if g=='KDM6A' else g for g in genes],fontsize=8)
    ax.set_xticks(np.arange(8)+.5,list('12341234'),fontsize=7.5)
    ax.tick_params(length=0,pad=3)
    for spine in ax.spines.values():spine.set_visible(False)
    ax.axvline(4,color='white',lw=2.8)
    for x,label,col in [(2,'CTRL',C.CTRL),(6,'LOY',C.LOY)]:
        ax.text(x,-.20,label,ha='center',va='bottom',fontsize=8.5,color=col)
    for i,gene in enumerate(genes):
        for j in range(8):
            value=matrix[i,j]
            zero=bool(lookup.loc[(gene,['CTRL','LOY'][j//4],j%4+1),'all_samples_zero'])
            jackpot=(gene,j%4+1) in D0_JACKPOTTED_GUIDES
            ax.text(j+.5,i+.5,f'{value:.1f}'+('†' if jackpot else '')+('*' if zero else ''),ha='center',va='center',
                    fontsize=7.5,color=C.DARK if value>10 else 'white')
    ax.set_xlabel(D0_GUIDE_LABEL,fontsize=8,labelpad=5)
    data['plotted_log2_cpm_plus_1']=data.log2_cpm_plus_1
    data.to_csv(DATA/f'supp_panel_{panel}_d0.csv',index=False)
    return mesh


def make_support():
    before=input_hashes(); height=8.6
    fig=plt.figure(figsize=(9.6,height))
    C.letter(fig,'a',.018,8.48/height)
    C.culture_legend(fig,(.095,8.49/height),size=7.5,ncol=4)
    for i,gene in enumerate(C.EXTENDED):
        row,col=divmod(i,5)
        ax=fig.add_axes([.075+col*.18,(7.15 if row==0 else 5.72)/height,.147,.94/height])
        C.guide_panel(ax,gene,show_y=col==0,small=True,scale='cpm_plus_1',linear_max=60)
        ax.set_yticks([0,20,40,60])
    fig.text(.018,7.00/height,'Endpoint CPM + 1',rotation=90,va='center',fontsize=8)
    fig.text(.51,5.38/height,'Brunello sgRNA',ha='center',fontsize=8)
    C.letter(fig,'b',.018,5.15/height)
    C.culture_legend(fig,(.095,5.23/height),size=7.5,ncol=4)
    grey_box=[.062,.83/height,.299,4.15/height]
    fig.add_artist(Rectangle(grey_box[:2],grey_box[2],grey_box[3],transform=fig.transFigure,
                            facecolor='#F3F3F3',edgecolor='#BFBFBF',linewidth=.7,zorder=-1))
    for row,genes in enumerate(GAMETOLOG_ROWS):
        for col,gene in enumerate(genes):
            ax=fig.add_axes([.085+col*.15,(3.53 if row==0 else 1.04)/height,.115,1.25/height])
            C.guide_panel(ax,gene,show_y=col==0,small=True,scale='cpm_plus_1',linear_max=60)
            ax.set_yticks([0,20,40,60])
            if col<2:ax.set_facecolor('none')
            if gene=='KDM6A':ax.set_title('KDM6A (UTX)',pad=5,fontsize=8.5)
    fig.text(.018,2.65/height,'Endpoint CPM + 1',rotation=90,va='center',fontsize=8)
    fig.text(.29,.55/height,'Brunello sgRNA',ha='center',fontsize=8)
    fig.text(.075,.20/height,'* Zero reads in all six libraries',fontsize=7)
    C.letter(fig,'c',.505,5.15/height)
    saga=pd.read_csv(C.DATA/'saga_reference_scores.csv')
    schematic=fig.add_axes([.535,3.07/height,.45,1.98/height])
    schematic_meta=draw_saga_membership(schematic,saga,DATA)
    for text in schematic.texts:
        if text.get_text()=='DUB':text.set_position((.04,.22))
    saga_columns=[saga.iloc[:13],saga.iloc[13:]]
    for col,rows in enumerate(saga_columns):
        ax=fig.add_axes([.62+col*.23,.73/height,.10,2.10/height])
        C.score_plot(ax,rows.source_gene.tolist(),(-5.6,2.8),font_size=7.5,
                     row_labels=rows.gene.tolist(),separators=[8.5] if col==0 else [3.5,4.5,6.5],tick_values=[-4,-2,0,2])
        for label,fdr in zip(ax.get_yticklabels(),rows.fdr_synth):
            label.set_fontweight('semibold' if fdr<.05 else 'normal')
        ax.set_ylim(12.55,-.65)
        if col==1:ax.axhspan(6.5,10.5,color='#F0F0F0',zorder=-1)
    C.fdr_legend(fig,(.55,.17/height),size=6.4,ncol=3)
    layout=check_text(fig,MAIN)
    for panel,genes in [('a',C.EXTENDED),('b',C.SUPP_ENDPOINT)]:
        points=C.PROFILES.loc[C.PROFILES.timepoint.eq('Endpoint') & C.PROFILES.gene.isin(genes)].copy()
        points['plotted_cpm_plus_1']=points.cpm+1
        points.to_csv(DATA/f'main_panel_{panel}_endpoint.csv',index=False)
    C.export(fig,MAIN)
    assert input_hashes()==before
    record=C.RENDER_RECORD[MAIN]
    record.update({'role':'main','figure_number':'Fig 4.3','panels':list('abc'),
        'panel_mapping':{'a':'former Fig 4.3d','b':'former Fig 4.3g','c':'former Fig 4.3f'},
        'panel_a_genes':C.EXTENDED,'panel_a_rows':[C.EXTENDED[:5],C.EXTENDED[5:]],
        'panel_a_points':160,'panel_b_rows':GAMETOLOG_ROWS,'panel_b_points':96,
        'panel_b_legend':{'colours':{'CTRL':C.CTRL,'LOY':C.LOY},'markers':{'R1':'circle','R2':'triangle'}},
        'panel_b_grey_box_genes':C.SUPP_ENDPOINT[:4],'endpoint_scale':'linear CPM + 1',
        'endpoint_ylim':[0,60],'zero_raw_counts_at':1,'panel_c_schematic':schematic_meta,
        'panel_c_score_columns':[s.gene.tolist() for s in saga_columns],
        'source_hashes':before,'source_hashes_unchanged':True,'text_layout':layout})
    return record


def make_method_d0_support():
    before=input_hashes(); height=7.0
    fig=plt.figure(figsize=(9.6,height))
    C.letter(fig,'a',.018,6.85/height);C.letter(fig,'b',.53,6.85/height)
    comparisons=draw_comparisons(fig,height,
        axes_bounds=[[.075,4.05/height,.35,2.35/height],[.585,4.05/height,.35,2.35/height]],
        notes_positions=[[.075,3.55/height],[.585,3.55/height]])
    C.fdr_legend(fig,(.63,6.89/height),size=7,ncol=3)
    C.letter(fig,'c',.018,2.95/height);C.letter(fig,'d',.53,2.95/height)
    mesh=baseline_heatmap(fig,C.PANEL_B,[.18,.65/height,.31,1.98/height],'c')
    baseline_heatmap(fig,C.SUPP_ENDPOINT,[.68,1.55/height,.29,1.08/height],'d')
    bar=fig.add_axes([.68,.72/height,.29,.09/height])
    cb=fig.colorbar(mesh,cax=bar,orientation='horizontal',ticks=[0,5,10,15])
    cb.outline.set_visible(False);cb.ax.tick_params(labelsize=7,length=2,pad=2)
    cb.set_label('D0 log₂(CPM + 1); panels c,d',fontsize=8,labelpad=4)
    fig.text(.18,.10/height,D0_FOOTNOTE,fontsize=7)
    layout=check_text(fig,SUPPLEMENT)
    C.export(fig,SUPPLEMENT)
    assert input_hashes()==before
    record=C.RENDER_RECORD[SUPPLEMENT]
    record.update({'role':'supp','figure_number':'Supp Fig 4.3','panels':list('abcd'),
        'panel_mapping':{'a':'former Fig 4.3a','b':'former Fig 4.3b','c':'former Fig 4.3c','d':'former Fig 4.3e'},
        'method_comparisons':comparisons,'panel_c_genes':C.PANEL_B,'panel_d_genes':C.SUPP_ENDPOINT,
        'panel_c_cells':88,'panel_d_cells':48,'baseline_scale':'log2(CPM + 1)',
        'baseline_guide_label':D0_GUIDE_LABEL,'baseline_footnote':D0_FOOTNOTE,
        'dagger_guides':[{'gene':gene,'guide_number':guide} for gene,guide in sorted(D0_JACKPOTTED_GUIDES)],
        'baseline_clim':[0,15],'source_hashes':before,'source_hashes_unchanged':True,'text_layout':layout})
    return record


def main():
    records={MAIN:make_support(),SUPPLEMENT:make_method_d0_support()}
    manifest=C.DATA/'render_manifest.json'
    existing=json.loads(manifest.read_text(encoding='utf-8')) if manifest.exists() else {}
    existing.update(records)
    manifest.write_text(json.dumps(existing,indent=2)+'\n',encoding='utf-8')
    (DATA/'render_manifest.json').write_text(json.dumps(records,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({name:{'inches':r['inches'],'panels':r['panels'],'text_layout':r['text_layout']} for name,r in records.items()},indent=2))
