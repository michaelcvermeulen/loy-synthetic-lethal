"""Candidate-discovery displays of frozen drugZ, MAGeCK and guide-count data."""
from pathlib import Path
import os, sys, json
THESIS = Path(__file__).resolve().parents[3]
DEPS = THESIS.parent / '.figure_python_deps'
if DEPS.is_dir(): sys.path.insert(0, str(DEPS))
OUT = THESIS / 'figures/chapter4'
DATA = OUT / 'source_data/ht1080_candidate_discovery'
os.environ.setdefault('MPLCONFIGDIR', str(DATA / 'mpl_cache'))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.text import Text
from matplotlib.ticker import NullLocator, ScalarFormatter
sys.path.insert(0, str(THESIS / 'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style, condition_colours
apply_style(9)
COLOURS = condition_colours()
CTRL, LOY = '#3C5488', '#E64B35'
STRONG, EXPLORATORY, NS = '#7A3E85', '#BC8A26', '#B8B8B8'
DARK = '#222222'
MAIN = 'ch04_main_ht1080_candidate_discovery'
# Promoted to main Fig 4.3; SUPP remains an internal compatibility alias.
SUPP = 'ch04_main_ht1080_candidate_support'
GUIDES = 'ch04_supp_ht1080_gametolog_guides'
GAMETOLOGS = ['DDX3X','EIF1AX','RPS4X','ZFX']
FOLLOWUP = ['DDX3X', 'SREBF1', 'TADA2B', 'TAF6L']
# Supp2b D0 heatmap: ten top-ranked drugZ sensitizers plus TADA2B (author request, 22 Oct 2026).
PANEL_B = ['DDX3X', 'TIMM8A', 'MED12', 'SCLY', 'WTAP', 'SREBF2', 'SREBF1', 'LLGL1', 'TAF6L', 'UBE2G2', 'TADA2B']
SUPP_ENDPOINT = ['DDX3X', 'EIF1AX', 'ZFX', 'RPS4X', 'KDM5C', 'KDM6A']
EXTENDED = ['SREBF1', 'SREBF2', 'TADA2B', 'TAF6L', 'SUPT20H', 'ATXN7L3', 'TADA1', 'TIMM8A', 'MED12', 'EIF4A3']
OTHER = ['EIF1AX', 'ZFX', 'RPS4X']
DZ = pd.read_csv(DATA/'all_gene_scores.csv')
INDEXED = DZ.set_index('GENE')
PROFILES = pd.read_csv(DATA/'guide_profiles.csv')
SUMMARY = json.loads((DATA/'summary.json').read_text(encoding='utf-8'))
RENDER_RECORD = {}

def letter(fig, value, x, y):
    fig.text(x, y, value, fontsize=13, fontweight='semibold', va='top')

def fdr_text(value):
    if value < .001: return f'{value:.3g}'
    if value < .01: return f'{value:.3g}'
    return f'{value:.3f}'

def colour(fdr):
    return STRONG if fdr < .05 else EXPLORATORY if fdr < .25 else NS

def fdr_legend(fig, anchor, size=8, ncol=3):
    handles=[Line2D([],[],color=colour(f),marker='o',ms=4,lw=0,label=label)
        for f,label in [(.01,'FDR < 0.05'),(.1,'0.05 ≤ FDR < 0.25'),(.9,'FDR ≥ 0.25')]]
    fig.legend(handles=handles,loc='upper left',bbox_to_anchor=anchor,ncol=ncol,
        fontsize=size,handletextpad=.35,columnspacing=1.2,borderaxespad=0)

def culture_legend(fig, anchor, size=8, ncol=4):
    handles=[Line2D([],[],color=c,marker='o',ms=4,lw=0,label=s)
             for s,c in [('CTRL',CTRL),('LOY',LOY)]]
    handles += [Line2D([],[],color='#666666',marker=m,ms=4,lw=0,label=s)
                for s,m in [('R1','o'),('R2','^')]]
    fig.legend(handles=handles,loc='upper left',bbox_to_anchor=anchor,ncol=ncol,
               fontsize=size,handletextpad=.3,columnspacing=1,borderaxespad=0)

def score_plot(ax, genes, limits, fdr_field='fdr_synth', font_size=9, separators=(),
               row_labels=None, tick_values=None):
    subset = INDEXED.loc[genes]
    y=np.arange(len(genes))
    vals=subset.normZ.to_numpy()
    fdrs=subset[fdr_field].to_numpy()
    ax.axvline(0,color='#777777',lw=.65,zorder=0)
    ax.hlines(y,0,vals,color='#D5D5D5',lw=.8,zorder=1)
    ax.scatter(vals,y,s=27,c=[colour(f) for f in fdrs],edgecolor='white',lw=.3,zorder=3)
    ax.set_yticks(y, row_labels if row_labels is not None else genes,fontsize=font_size)
    ax.tick_params(axis='y',length=0,pad=5)
    ax.set_ylim(len(genes)-.45,-.65)
    ax.set_xlim(*limits)
    if tick_values is not None: ax.set_xticks(tick_values)
    ax.tick_params(axis='x',labelsize=font_size)
    ax.spines['left'].set_visible(False)
    ax.set_xlabel('drugZ normZ',fontsize=font_size,labelpad=3)
    for yy in separators: ax.axhline(yy,color='#DDDDDD',lw=.65)
    for yy,fdr in zip(y,fdrs):
        ax.text(1.035,yy,fdr_text(fdr),transform=ax.get_yaxis_transform(),
                va='center',fontsize=font_size,clip_on=False)
    ax.text(1.035,1.055,'FDR',transform=ax.transAxes,fontsize=font_size,ha='left')
    return subset

def guide_panel(ax,gene,show_y=True,small=False,scale='log2_cpm_plus_1',linear_max=50):
    assert scale in ['log2_cpm_plus_1','cpm_plus_1']
    sub=PROFILES.loc[PROFILES.gene.eq(gene) & PROFILES.timepoint.eq('Endpoint')]
    offsets={('CTRL',1):-.22,('CTRL',2):-.075,('LOY',1):.075,('LOY',2):.22}
    for (genotype,rep),rows in sub.groupby(['genotype','replicate']):
        values=rows.cpm+1 if scale=='cpm_plus_1' else rows.log2_cpm_plus_1
        ax.scatter(rows.guide_number+offsets[(genotype,int(rep))],values,
            color=CTRL if genotype=='CTRL' else LOY,marker='o' if rep==1 else '^',
            s=17 if small else 22,edgecolors='white',linewidths=.3,zorder=3)
    zeros=set(sub.loc[sub.all_samples_zero,'guide_number'])
    ax.set_xticks([1,2,3,4],[f'{i}*' if i in zeros else str(i) for i in range(1,5)])
    ax.set_xlim(.5,4.5)
    if scale=='cpm_plus_1':
        ax.set_ylim(0,linear_max); ax.set_yticks(np.arange(0,linear_max+1,10))
        ax.axhline(1,color='#DDDDDD',lw=.6,zorder=0)
    else:
        ax.set_ylim(-.17,6.25); ax.set_yticks([0,2,4,6])
        ax.axhline(0,color='#DDDDDD',lw=.6,zorder=0)
    ax.set_title(gene,pad=5,fontsize=8.5 if small else 9)
    ax.tick_params(labelsize=8 if small else 8.5,length=2.5,pad=2)
    if not show_y: ax.tick_params(axis='y',labelleft=False)

def export(fig,stem):
    review=OUT/'review'; review.mkdir(exist_ok=True)
    fig.savefig(OUT/f'{stem}.pdf',dpi=300)
    fig.savefig(OUT/f'{stem}.svg',dpi=300)
    fig.savefig(OUT/f'{stem}.png',dpi=600)
    fig.savefig(review/f'{stem}_preview.png',dpi=160)
    RENDER_RECORD[stem]={'inches':list(fig.get_size_inches()),'axes':len(fig.axes),'font':'Open Sans',
                         'png_dpi':600,'preview_dpi':160}
    plt.close(fig)

def resize_height_from_bottom(fig,new_height):
    """Add/remove bottom canvas while preserving existing physical panel sizes."""
    width,old_height=fig.get_size_inches()
    shift=new_height-old_height
    boxes=[(ax,ax.get_position().bounds) for ax in fig.axes]
    legends=[(legend,legend.get_bbox_to_anchor().transformed(fig.transFigure.inverted()).bounds)
             for legend in fig.legends]
    texts=[(t,t.get_position()) for t in fig.texts]
    lines=[(a,np.asarray(a.get_ydata(),dtype=float)) for a in fig.artists if isinstance(a,Line2D)]
    fig.set_size_inches(width,new_height,forward=False)
    for ax,(x,y,w,h) in boxes:
        ax.set_position([x,(y*old_height+shift)/new_height,w,h*old_height/new_height])
    for legend,(x,y,w,h) in legends:
        legend.set_bbox_to_anchor((x,(y*old_height+shift)/new_height,w,h*old_height/new_height),transform=fig.transFigure)
    for t,(x,y) in texts: t.set_position((x,(y*old_height+shift)/new_height))
    for line,y in lines: line.set_ydata((y*old_height+shift)/new_height)

def main_figure():
    # Four panels: rank overview, guide distributions and complete STRING network.
    from make_ht1080_candidate_layouts import make_main
    RENDER_RECORD[MAIN] = make_main()

def supplementary_figure():
    from make_ht1080_support_compact import make_support
    RENDER_RECORD[SUPP] = make_support()

def gametolog_figure():
    fig=plt.figure(figsize=(9.6,7.6))
    table=pd.read_csv(DATA/'gametolog_guide_table.csv')
    culture_legend(fig,(.34,.674),size=8)
    fig.text(.65,.98,'One D0 sample per genotype',fontsize=8)
    for letter_text,yy in [('a',.985),('b',.663),('c',.347)]: letter(fig,letter_text,.018,yy)
    for i,gene in enumerate(GAMETOLOGS):
        left=.09+i*.225
        p=fig.add_axes([left,.735,.177,.20])
        sub=PROFILES.loc[PROFILES.gene.eq(gene)&PROFILES.timepoint.eq('D0')]
        for genotype,rows in sub.groupby('genotype'):
            p.scatter(rows.guide_number+(-.1 if genotype=='CTRL' else .1),rows.log2_cpm_plus_1,
                color=CTRL if genotype=='CTRL' else LOY,s=23,edgecolors='white',linewidths=.3)
        p.set_ylim(-.25,15.25);p.set_yticks([0,5,10,15]);p.set_xlim(.5,4.5)
        p.set_xticks([1,2,3,4],['1','2','3','4*' if gene=='DDX3X' else '4'])
        p.set_title(gene,fontsize=9,pad=5);p.tick_params(labelsize=8,length=2.5,pad=2)
        if i: p.tick_params(axis='y',labelleft=False)
        p=fig.add_axes([left,.42,.177,.20]);guide_panel(p,gene,show_y=i==0)
        p=fig.add_axes([left,.105,.177,.20])
        rows=table.loc[table.Gene.eq(gene) & table.mageck_LFC_displayed]
        p.axhline(0,color='#888888',lw=.65)
        p.vlines(rows.guide_number,0,rows.mageck_reported_LFC,color='#C4C4C4',lw=.8)
        p.scatter(rows.guide_number,rows.mageck_reported_LFC,color=DARK,s=22,zorder=3)
        p.set_xlim(.5,4.5);p.set_ylim(-7.4,2.8);p.set_yticks([-6,-4,-2,0,2])
        p.set_xticks([1,2,3,4],['1','2','3','4*' if gene=='DDX3X' else '4'])
        p.set_title(gene,fontsize=9,pad=5);p.tick_params(labelsize=8,length=2.5,pad=2)
        if gene=='DDX3X': p.text(4,-3.3,'No reads',rotation=90,ha='center',va='center',fontsize=7.5,color='#777777')
        if i: p.tick_params(axis='y',labelleft=False)
    for yy,label in [(.835,'D0 log₂(CPM + 1)'),(.52,'Endpoint log₂(CPM + 1)'),(.205,'MAGeCK log₂FC (LOY / CTRL)')]:
        fig.text(.018,yy,label,rotation=90,va='center',fontsize=8.5)
    for yy in [.698,.37,.055]: fig.text(.51,yy,'Guide',ha='center',fontsize=8.5)
    fig.text(.94,.025,'* DDX3X guide 4: zero reads in all six samples; LFC omitted',ha='right',fontsize=7.5)
    export(fig,GUIDES)

if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser()
    group=parser.add_mutually_exclusive_group()
    group.add_argument('--main-only',action='store_true')
    group.add_argument('--main-and-support',action='store_true')
    group.add_argument('--support-only',action='store_true')
    args=parser.parse_args()
    manifest=DATA/'render_manifest.json'
    if manifest.exists(): RENDER_RECORD.update(json.loads(manifest.read_text(encoding='utf-8')))
    if not args.support_only: main_figure()
    if not args.main_only:
        supplementary_figure()
        if not (args.main_and_support or args.support_only): gametolog_figure()
    manifest.write_text(json.dumps(RENDER_RECORD,indent=2),encoding='utf-8')
    print(json.dumps(RENDER_RECORD,indent=2))
