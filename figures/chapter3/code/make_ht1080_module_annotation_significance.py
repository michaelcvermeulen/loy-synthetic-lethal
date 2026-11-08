"""Display frozen top-annotation tests with descriptive names and tested terms."""
from pathlib import Path
import hashlib
import json
import re
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.transforms import blended_transform_factory
from matplotlib.font_manager import FontProperties
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd

FIG=Path(__file__).resolve().parents[1]
ROOT=FIG.parents[1]
SOURCE=FIG/'source_data/ht1080_expanded_string_modules'
sys.path.insert(0,str(ROOT/'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
apply_style(7)
STEM='ch03_supp_ht1080_module_annotation_significance'
inputs=[SOURCE/'annotation_significance_plot_data.csv',SOURCE/'module_labels_reviewed.json']
hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}
data=pd.read_csv(inputs[0],float_precision='round_trip')
labels=json.loads(inputs[1].read_text(encoding='utf-8'))
data['module_label']=data.module_id.map(labels)
assert len(data)==42 and data.module_label.notna().all()
assert data.groupby('direction').size().to_dict()=={'DOWN':24,'UP':18}
assert data.background_genes.eq(14564).all()

def probability_text(value):
    return f'{value:.2e}' if value<.001 else f'{value:.3g}'

def display_term(pathway):
    source,term=('GO:BP',pathway.removeprefix('GOBP_')) if pathway.startswith('GOBP_') else ('Reactome',pathway.removeprefix('REACTOME_'))
    term=term.replace('_',' ').lower().replace('signaling','signalling')
    for word in ['MAPK','L1','GPCR','IGF','IGFBPS','BAF','ESBAF','BMP','WNT','MITF','K']:
        formatted={'IGFBPS':'IGFBPs','ESBAF':'esBAF'}.get(word,word)
        term=re.sub(r'\b'+word.lower()+r'\b',formatted,term)
    term=term.replace('pirna','piRNA').replace('p75ntr','p75NTR')
    return source,term

data['annotation_source'],data['annotation_display']=zip(*data.pathway.map(display_term))
WIDTH_PT=183/25.4*72
FONT_MODULE=6.0
FONT_TERM=5.6
LINE_PT=6.7
LABEL_X=.025
TERM_X=.036
TEXT_RIGHT=.565
fig=plt.figure(figsize=(183/25.4,260/25.4),facecolor='white')
fig.canvas.draw()
measure=fig.canvas.get_renderer()
term_font=FontProperties(family='Open Sans',size=FONT_TERM)

def wrap_term(text):
    width=(TEXT_RIGHT-TERM_X)*fig.bbox.width
    lines=[]
    current=''
    for word in text.split():
        candidate=(current+' '+word).strip()
        if current and measure.get_text_width_height_descent(candidate,term_font,False)[0]>width:
            lines.append(current)
            current=word
        else:current=candidate
    lines.append(current)
    return lines

data['annotation_lines']=[wrap_term(f'{r.annotation_source}: {r.annotation_display}') for r in data.itertuples()]
data['row_height_pt']=data.annotation_lines.map(lambda lines:LINE_PT*(1+len(lines))+1.7)
FIRST_TOP=50
PANEL_GAP=40
FOOT=37
HEIGHT_PT=FIRST_TOP+data.row_height_pt.sum()+PANEL_GAP+FOOT
fig.set_size_inches(WIDTH_PT/72,HEIGHT_PT/72)
label_texts=[]
row_texts=[]

def panel(table,top_pt,color,letter,title):
    total=table.row_height_pt.sum()
    rect=[.625,1-(top_pt+total)/HEIGHT_PT,.160,total/HEIGHT_PT]
    ax=fig.add_axes(rect)
    ax.set_ylim(total,0)
    ax.set_xlim(-.6,27)
    ax.set_xticks([0,10,20])
    ax.tick_params(axis='x',labelsize=6,length=2.5,width=.5,color='#58636B')
    ax.set_yticks([])
    for side in ['left','right','top']:ax.spines[side].set_visible(False)
    ax.spines['bottom'].set_linewidth(.5)
    ax.spines['bottom'].set_color('#58636B')
    ax.grid(axis='x',color='#E8ECEF',linewidth=.45)
    ax.set_axisbelow(True)
    ax.axvline(-np.log10(.05),color='#747F85',linewidth=.6,linestyle=(0,(3,2)),zorder=1)
    trans=blended_transform_factory(fig.transFigure,ax.transData)
    cursor=0
    for i,row in enumerate(table.itertuples()):
        y=cursor+row.row_height_pt/2
        if i%2==0:
            fig.add_artist(Rectangle((.020,1-(top_pt+cursor+row.row_height_pt)/HEIGHT_PT),.970,row.row_height_pt/HEIGHT_PT,transform=fig.transFigure,facecolor='#F5F7F8',edgecolor='none',zorder=-1))
            ax.axhspan(cursor,cursor+row.row_height_pt,color='#F5F7F8',zorder=0)
        ax.plot([row.neglog10_q,row.neglog10_P],[y,y],color=color,alpha=.6,linewidth=.7,zorder=2)
        ax.scatter(row.neglog10_P,y,s=16,facecolors='white',edgecolors=color,linewidths=.8,marker='o',zorder=4)
        ax.scatter(row.neglog10_q,y,s=13,facecolors=color,edgecolors=color,linewidths=.4,marker='D',zorder=5)
        label=f'{row.display_id}  {row.module_label}'
        t=ax.text(LABEL_X,cursor+3.35,label,transform=trans,fontsize=FONT_MODULE,va='center',ha='left',color='#202B32',clip_on=False)
        label_texts.append(t)
        for j,line in enumerate(row.annotation_lines):
            t=ax.text(TERM_X,cursor+3.35+LINE_PT*(j+1),line,transform=trans,fontsize=FONT_TERM,va='center',ha='left',color='#52616C',clip_on=False)
            label_texts.append(t)
        ax.text(.602,y,str(row.module_genes),transform=trans,fontsize=6,va='center',ha='right',color='#566169',clip_on=False)
        ax.text(.881,y,probability_text(row.P_value),transform=trans,fontsize=6,va='center',ha='right',clip_on=False)
        ax.text(.984,y,probability_text(row.BH_all_modules_and_terms),transform=trans,fontsize=6,va='center',ha='right',clip_on=False)
        row_texts.append((row.module_id,label_texts[-(1+len(row.annotation_lines)):],cursor,row.row_height_pt,ax))
        cursor+=row.row_height_pt
    fig.text(.006,1-(top_pt-21)/HEIGHT_PT,letter,fontsize=9,weight='semibold',color='#202B32')
    fig.text(.025,1-(top_pt-21)/HEIGHT_PT,title,fontsize=7,weight='semibold',color=color)
    for x,label,align in [(.025,'Module / top GO:BP or Reactome term','left'),(.602,'Genes','right'),(.881,'P value','right'),(.984,'BH FDR','right')]:
        fig.text(x,1-(top_pt-6)/HEIGHT_PT,label,fontsize=6,weight='semibold',ha=align,color='#39464E')
    return ax

handles=[Line2D([0],[0],marker='o',markersize=4,markerfacecolor='white',markeredgecolor='#293A44',linestyle='none'),
         Line2D([0],[0],marker='D',markersize=3.5,markerfacecolor='#293A44',markeredgecolor='#293A44',linestyle='none'),
         Line2D([0],[0],color='#747F85',linewidth=.6,linestyle=(0,(3,2)))]
fig.legend(handles,['Annotation P value','Annotation FDR','P or FDR = 0.05'],loc='upper left',bbox_to_anchor=(.02,1-3/HEIGHT_PT),ncol=3,frameon=False,fontsize=6.5,columnspacing=2,handletextpad=.6)
down=data[data.direction.eq('DOWN')]
up=data[data.direction.eq('UP')]
panel(down,FIRST_TOP,'#0072B2','a','Lower expression · 24 modules')
panel(up,FIRST_TOP+down.row_height_pt.sum()+PANEL_GAP,'#D55E00','b','Higher expression · 18 modules')
fig.text(.705,15/HEIGHT_PT,'−log₁₀(P or FDR)',fontsize=6,ha='center')
fig.text(.025,3/HEIGHT_PT,'Annotation BH adjustment: 42 modules × 6,203 terms; expressed STRING-mapped background.',fontsize=5.8,color='#53616A')
fig.canvas.draw()
renderer=fig.canvas.get_renderer()
width_limit=fig.transFigure.transform((TEXT_RIGHT+.006,0))[0]
assert all(t.get_window_extent(renderer).x1<width_limit for t in label_texts),[(t.get_text(),t.get_window_extent(renderer).x1) for t in label_texts if t.get_window_extent(renderer).x1>=width_limit]
for mid,texts,start,height,ax in row_texts:
    lower=ax.transData.transform((0,start+height))[1]
    upper=ax.transData.transform((0,start))[1]
    assert all(t.get_window_extent(renderer).y0>=lower and t.get_window_extent(renderer).y1<=upper for t in texts),mid
    for t,u in zip(texts,texts[1:]):
        assert not t.get_window_extent(renderer).overlaps(u.get_window_extent(renderer)),(mid,t.get_text(),u.get_text())
for ext in ['pdf','svg','png']:fig.savefig(FIG/f'{STEM}.{ext}',dpi=600)
fig.savefig(FIG/'previews'/f'{STEM}_preview.png',dpi=180)
plt.close(fig)
data.drop(columns=['annotation_lines','row_height_pt']).to_csv(SOURCE/'annotation_significance_reviewed_plot_data.csv',index=False,float_format='%.17g')
assert all(hashlib.sha256(p.read_bytes()).hexdigest()==hashes[p.name] for p in inputs)
(SOURCE/'annotation_label_review_provenance.json').write_text(json.dumps({'date':'2026-10-09','input_sha256':hashes,'module_count':42,'changes':'Reviewed descriptive labels; NLGN4X dagger and explanatory note removed at author request; exact saved top annotation below each label, with GO:BP/Reactome source. All saved probabilities, ranks and row memberships retained. Display changes only: underscores replaced by spaces, British signalling and acronym capitalisation.','family':'42 modules x 6203 terms = 260526 tests','font':'Open Sans','size_mm':[183,HEIGHT_PT/72*25.4],'analysis_recomputed':False},indent=2)+'\n',encoding='utf-8')
print(f'S17: rendered 42 unchanged annotation tests with top terms; 183 x {HEIGHT_PT/72*25.4:.1f} mm.')
