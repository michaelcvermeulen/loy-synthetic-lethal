"""Editable full-candidate STRING draft; frozen data, no network requests."""
from pathlib import Path
import sys, os, json, hashlib
THESIS = Path(__file__).resolve().parents[3]
sys.path.insert(0,str(THESIS.parent/'.figure_python_deps'))
OUT = THESIS/'figures/chapter4'
DATA = OUT/'source_data/ht1080_all_sl_string'
os.environ.setdefault('MPLCONFIGDIR',str(DATA/'mpl_cache'))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import to_rgb
from matplotlib.patches import Ellipse
sys.path.insert(0,str(THESIS/'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
apply_style(9)

STEM = 'ch04_main_ht1080_all_sl_string_draft'
NODES = pd.read_csv(DATA/'all_nodes.csv').set_index('GENE',drop=False)
FUNC = pd.read_csv(DATA/'functional_edges.csv')
PHYS = pd.read_csv(DATA/'physical_edges.csv').set_index('pair')
STRONG, EXPLORATORY = '#7A3E85', '#BC8A26'
CONNECTED = set(FUNC.gene_a)|set(FUNC.gene_b)
ISOLATED = sorted(set(NODES.index)-CONNECTED)

# Fixed display coordinates spread labels and preserve every database edge.
# Positions are illustrative: no formal community or enrichment inference.
POS = {
    'TAF5L':(.80,5.45), 'TAF6L':(2.10,5.45),
    'SUPT20H':(.20,4.55), 'TADA2B':(2.65,4.55),
    'ATXN7L3':(.75,3.65), 'TADA1':(2.05,3.65),
    'TP53':(3.30,3.10), 'CREBBP':(4.60,4.30),
    'MED12':(4.00,5.55), 'SREBF1':(5.65,5.25),
    'SREBF2':(6.80,4.30), 'SIRT6':(7.00,5.55),
    'AGPAT2':(6.50,3.10), 'ACSL4':(7.30,2.25),
    'KMT2D':(4.95,2.65), 'BCORL1':(6.10,1.60),
    'SMARCA4':(3.35,1.50), 'PELO':(3.20,.35),
    'AMBRA1':(4.65,.35), 'TIMM8A':(6.20,.35),
    'DDX3X':(1.75,2.25), 'EIF4A3':(.15,2.10),
    'ZC3H18':(.70,1.00), 'WTAP':(.00,.20),
    'KIF2A':(2.00,.65), 'RPGRIP1L':(1.40,-.05),
}
assert set(POS)==CONNECTED

def colour(g): return STRONG if NODES.loc[g,'fdr_synth']<.05 else EXPLORATORY
def pale(c): return tuple(.89+.11*v for v in to_rgb(c))
def edgestyle(pair):
    if pair not in PHYS.index:
        return {'color':'#BFC4C7','lw':.8,'ls':(0,(3,3)),'alpha':.85}
    high = PHYS.loc[pair,'score'] >= .7
    return {'color':'#59636A' if high else '#A6AEB3',
            'lw':1.5 if high else .8,'ls':'-','alpha':.86 if high else .8}

def draw_network(ax, fontsize=8.5, node_size=85, offset=9):
    """Draw every retained edge and connected node; markers are true circles."""
    ax.set(xlim=(-.6,7.85),ylim=(-.75,6.25)); ax.axis('off')
    ax.add_patch(Ellipse((1.40,4.55),3.4,2.75,facecolor='#F5F3F6',edgecolor='none',zorder=0))
    ax.text(1.40,6.02,'SAGA',ha='center',fontsize=fontsize-.5,color='#666666')
    order=sorted(FUNC.to_dict('records'),key=lambda e:(e['pair'] in PHYS.index,
        PHYS.loc[e['pair'],'score'] if e['pair'] in PHYS.index else 0))
    for e in order:
        x1,y1=POS[e['gene_a']];x2,y2=POS[e['gene_b']]
        style=edgestyle(e['pair'])
        style['lw']*=.85
        ax.plot([x1,x2],[y1,y2],zorder=1,**style)
    labels=[];positions=[]
    below={'ATXN7L3','TADA1','CREBBP','SREBF2','RPGRIP1L','PELO','AMBRA1','TIMM8A'}
    for gene,(x,y) in POS.items():
        ax.scatter([x],[y],s=node_size,c=colour(gene),marker='o',edgecolor='white',linewidth=.7,zorder=4)
        dy=-offset if gene in below else offset
        label=ax.annotate(gene,(x,y),xytext=(0,dy),textcoords='offset points',
            ha='center',va='top' if dy<0 else 'bottom',fontsize=fontsize,color='#252525',
            bbox={'facecolor':'white','edgecolor':'none','alpha':.9,'boxstyle':'square,pad=.03'},zorder=5)
        labels.append(label)
        positions.append({'GENE':gene,'group':'connected','x':x,'y':y,'label_dy_points':dy,'node_shape':'circle'})
    return labels,positions

def draw_isolated(side, fontsize=8):
    side.axis('off');side.set(xlim=(0,2),ylim=(-.4,11.6))
    labels=[];positions=[]
    for i,gene in enumerate(ISOLATED):
        col,row=divmod(i,11);x=col*1.04;y=10.7-row
        side.scatter(x+.04,y,s=19,c=colour(gene),marker='o',linewidth=0)
        labels.append(side.text(x+.15,y,gene,va='center',fontsize=fontsize))
        positions.append({'GENE':gene,'group':'unconnected','x':x,'y':y,'node_shape':'circle'})
    return labels,positions

def edge_legend(fig, anchor, fontsize=8):
    handles=[Line2D([],[],color='#59636A',lw=1.3,label='Physical score ≥0.7'),
        Line2D([],[],color='#A6AEB3',lw=.7,label='Physical score 0.4–<0.7'),
        Line2D([],[],color='#BFC4C7',lw=.7,ls=(0,(3,3)),label='Other functional association ≥0.4')]
    return fig.legend(handles=handles,loc='lower left',bbox_to_anchor=anchor,
        frameon=False,ncol=3,fontsize=fontsize,handlelength=2.2,columnspacing=1.7)

def check_labels(fig,labels):
    from matplotlib.text import Text
    fig.canvas.draw();renderer=fig.canvas.get_renderer();overlaps=[]
    for i,a in enumerate(labels):
        for b in labels[i+1:]:
            if Text.get_window_extent(a,renderer).overlaps(Text.get_window_extent(b,renderer)):
                overlaps.append([a.get_text(),b.get_text()])
    assert not overlaps,overlaps
    assert {a.get_text() for a in labels}==set(NODES.index) and len(labels)==48
    return overlaps

def make_standalone():
    fig=plt.figure(figsize=(9.6,6.35))
    ax=fig.add_axes([.025,.145,.735,.745])
    fig.text(.025,.975,'48 drugZ depleted candidates · FDR < 0.25',fontsize=10,fontweight='semibold',va='top')
    fig.text(.025,.939,'26 connected proteins · 66 associations · STRING v12.5',fontsize=8.5,color='#555555',va='top')
    labels,positions=draw_network(ax)
    side=fig.add_axes([.795,.16,.195,.73])
    fig.add_artist(Line2D([.777,.777],[.165,.927],transform=fig.transFigure,color='#E0E0E0',lw=.7))
    fig.text(.795,.939,'22 without a connection',fontsize=8.5,va='top')
    fig.text(.795,.914,'within this set at score ≥0.4',fontsize=7.7,color='#666666',va='top')
    other,other_positions=draw_isolated(side);labels+=other;positions+=other_positions
    handles=[Line2D([],[],color=STRONG,marker='o',ms=5,lw=0,label='drugZ FDR < 0.05 (28)'),
        Line2D([],[],color=EXPLORATORY,marker='o',ms=5,lw=0,label='0.05 ≤ FDR < 0.25 (20)')]
    fig.legend(handles=handles,loc='lower left',bbox_to_anchor=(.018,.068),frameon=False,ncol=2,fontsize=8,handletextpad=.5,columnspacing=1.5)
    edge_legend(fig,(.018,.028))
    fig.text(.99,.014,'No added proteins; positions are illustrative.',ha='right',fontsize=7.5,color='#666666')
    overlaps=check_labels(fig,labels)
    pd.DataFrame(positions).to_csv(DATA/'display_positions.csv',index=False)
    review=OUT/'review';review.mkdir(exist_ok=True)
    for ext in ['pdf','svg','png']:fig.savefig(OUT/f'{STEM}.{ext}',dpi=600)
    fig.savefig(review/f'{STEM}_preview.png',dpi=160);plt.close(fig)
    record={'stem':STEM,'inches':[9.6,6.35],'font':'Open Sans','png_dpi':600,
        'input_nodes':48,'drawn_labels':48,'node_shape':'circle','edges':len(FUNC),'physical_edges':len(PHYS),
        'physical_high_confidence_edges':int((PHYS.score>=.7).sum()),
        'functional_only_edges':int((~FUNC.pair.isin(PHYS.index)).sum()),
        'label_overlap_pairs':overlaps,'position_method':'Fixed manual display layout; no community inference',
        'input_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest()
          for p in [DATA/'all_nodes.csv',DATA/'functional_edges.csv',DATA/'physical_edges.csv']}}
    (DATA/'render_manifest.json').write_text(json.dumps(record,indent=2),encoding='utf-8')
    print(json.dumps(record,indent=2))

if __name__=='__main__':make_standalone()
