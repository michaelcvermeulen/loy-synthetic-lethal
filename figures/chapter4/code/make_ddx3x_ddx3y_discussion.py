"""Chapter 4 discussion figure; preserved a–c plus divergent-region panel d."""
from pathlib import Path
import sys, json, csv
THESIS=Path(__file__).resolve().parents[3]
PROJECT=THESIS.parent
ROOT=THESIS/'figures/chapter4/source_data/ddx3x_ddx3y_inhibitor_site_conservation'
sys.path.insert(0,str(PROJECT/'.figure_python_deps'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from matplotlib.backends.backend_pdf import PdfPages

for font in ['OpenSans-Regular.ttf','OpenSans-SemiBold.ttf']:
    font_manager.fontManager.addfont(str(THESIS/'skills/thesis-figure-style/assets/fonts'/font))
plt.rcParams.update({'font.family':'Open Sans','font.size':8,'axes.linewidth':.55,
    'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white'})
INK='#203345'; MUTED='#5D6B76'; TEAL='#267F8A'; TEAL_LIGHT='#DBEEF0'
PURPLE='#796296'; PURPLE_LIGHT='#EDE7F2'; ORANGE='#CC713E'; GAP='#B33D59'
GREY='#E6EAED'; Q='#B79A4C'; LIGHT='#F5F7F8'; DARK='#2D5964'

def read(name):
    with (ROOT/'data'/name).open(encoding='utf-8') as f: return list(csv.DictReader(f,delimiter='\t'))
rows=read('alignment_residue_map.tsv')
for r in rows:
    for k in ['x_pos','y_pos','alignment_column']: r[k]=int(r[k]) if r[k] else None
xmap={r['x_pos']:r for r in rows if r['x_pos']}
regions=read('region_similarity.tsv'); sites=read('binding_site_similarity.tsv')
audit=json.loads((ROOT/'data/analysis_audit.json').read_text())
OUT=THESIS/'figures/chapter4'
selected=json.loads((ROOT/'divergent_regions/selection_audit.json').read_text())['selected']

def text(ax,x,y,s,size=8,**kw):
    return ax.text(x,y,s,fontsize=size,color=kw.pop('color',INK),**kw)
def panel(fig,y,letter):
    fig.text(.035,y,letter,fontsize=11,fontweight='semibold',color=INK,va='top')
def box(ax,x,y,w,h,c,**kwargs):
    ax.add_patch(Rectangle((x,y),w,h,facecolor=c,edgecolor=kwargs.pop('edgecolor','none'),**kwargs))

fig=plt.figure(figsize=(183/25.4,192/25.4))
panel(fig,.965,'a')
ax=fig.add_axes([.16,.714,.795,.25]);ax.set_xlim(1,662);ax.set_ylim(0,3.3);ax.axis('off')
domain_ax=ax
for yv,g,length in [(2.43,'DDX3X',662),(1.42,'DDX3Y',660)]:
    box(ax,1,yv-.19,661,.38,GREY)
    for start,end,c in [(180,208,Q),(211,403,TEAL),(414,575,PURPLE)]:
        box(ax,start,yv-.23,end-start+1,.46,c)
    text(ax,-8,yv,g,8.2,ha='right',va='center',fontweight='semibold',clip_on=False)
    text(ax,670,yv,str(length),7.2,va='center',clip_on=False,color=MUTED)
    offset=0 if g=='DDX3X' else 2
    for start,end in [(211,403),(414,575)]:
        text(ax,(start+end)/2,yv,f'{start-offset}–{end-offset}',7.6,ha='center',va='center',color='white')
    for w in selected:
        box(ax,w['x_start'],yv-.28,w['columns'],.56,'none',edgecolor=ORANGE,lw=.85)
text(ax,85,3.01,'N-terminal region',7.4,ha='center')
text(ax,307,3.01,'Helicase ATP-binding',7.4,ha='center',color=TEAL)
text(ax,495,3.01,'Helicase C-terminal',7.4,ha='center',color=PURPLE)
text(ax,620,3.01,'C-terminal',7.1,ha='center')
text(ax,194,2.80,'Q',6.7,ha='center',color=Q,fontweight='semibold')
for g,yv in [('DDX3X',2.43),('DDX3Y',1.42)]:
    for s,e in [(200,207),(224,231)]: box(ax,s,yv+.26,e-s+1,.065,INK)
    for p in [274,276,302,323,326,357,503]: ax.plot(p,yv+.39,marker='o',ms=1.9,color=INK,markeredgewidth=0)
    for r in rows:
        if g=='DDX3Y' and r['x_pos'] and r['status']=='gap':
            ax.plot(r['x_pos'],yv,marker='|',ms=6,color='white',mew=1)
text(ax,348,1.94,'DEAD 347–350',6.7,ha='center',color=TEAL)
ax.plot([347,350],[2.16,2.16],color=TEAL,lw=1.2)
for r in rows:
    if r['status']!='identical':
        ax.vlines(r['x_pos'],.59,.87,color=GAP if r['status']=='gap' else ORANGE,lw=1)
text(ax,-8,.73,'Differences',7.5,ha='right',va='center',clip_on=False)
ax.plot([1,662],[.38,.38],color=MUTED,lw=.55)
for p in [1,100,200,300,400,500,600,662]:
    ax.plot([p,p],[.32,.38],color=MUTED,lw=.55);text(ax,p,.12,str(p),6.8,ha='center',color=MUTED)
fig.text(.16,.701,'DDX3X residue',fontsize=6.8,color=MUTED)
fig.text(.16,.677,'Orange: substitution    Red: gap    Bars: annotated ATP segments    Dots: modeled RNA-site contacts',fontsize=6.6,color=MUTED)

panel(fig,.627,'b')
ax=fig.add_axes([.375,.403,.43,.215]);ax.set_xlim(0,100);ax.set_ylim(-.6,4.75)
names=['N-terminal region','Q motif','Helicase ATP-binding domain','Helicase C-terminal domain','C-terminal region']
label_names=['N-terminal region','Q motif','ATP-binding domain','C-terminal helicase domain','C-terminal region']
for i,(name,label) in enumerate(zip(names,label_names)):
    yy=4-i;r=next(r for r in regions if r['region']==name)
    ident=float(r['identity_pct']);similar=float(r['similarity_pct'])
    c=TEAL if i==2 else PURPLE if i==3 else Q if i==1 else DARK
    ax.barh(yy,100,color=LIGHT,height=.57,zorder=0)
    ax.barh(yy,similar,color='#C5D8DB',height=.57,zorder=1)
    ax.barh(yy,ident,color=c,height=.57,zorder=2)
    text(ax,-4,yy,label+f"  ({r['x_start']}–{r['x_end']})",7.1,ha='right',va='center',clip_on=False)
    text(ax,104,yy,f'{ident:.1f}',7.3,ha='left',va='center',fontweight='semibold',clip_on=False)
    text(ax,121,yy,f'{similar:.1f}',7.3,ha='left',va='center',color=MUTED,clip_on=False)
text(ax,104,4.67,'Identity',6.6,ha='left',clip_on=False,color=MUTED)
text(ax,121,4.67,'Similar',6.6,ha='left',clip_on=False,color=MUTED)
ax.set_xticks([0,25,50,75,100]);ax.tick_params(axis='x',labelsize=6.8,length=2,color=MUTED,labelcolor=MUTED,pad=2)
ax.set_yticks([])
for s in ['top','left','right']:ax.spines[s].set_visible(False)
ax.spines['bottom'].set_color(MUTED)
fig.text(.375,.351,'Alignment columns (%)',fontsize=6.8,color=MUTED)

panel(fig,.320,'c')
left=fig.add_axes([.055,.035,.435,.265]);right=fig.add_axes([.535,.035,.435,.265])
for a in [left,right]:a.set_xlim(0,100);a.set_ylim(0,100);a.axis('off')
text(left,0,96,'ATP-binding cleft',8.5,fontweight='semibold',color=TEAL)
text(left,0,85,'RK-33-associated region',7.0)
text(left,18,73,'200–207',6.8,color=MUTED,ha='left')
text(left,63,73,'224–231',6.8,color=MUTED,ha='left')
for s,startx in [(200,18),(224,63)]:
    for j,p in enumerate(range(s,s+8)):
        xx=startx+j*4.2
        for yv,g in [(58,'x_aa'),(44,'y_aa')]:
            box(left,xx-1.8,yv-5,3.8,10,TEAL_LIGHT)
            text(left,xx,yv,xmap[p][g],8.0,ha='center',va='center',fontfamily='DejaVu Sans Mono')
        if p==230:
            for yv in [58,44]:box(left,xx-1.8,yv-5,3.8,10,'none',edgecolor=ORANGE,lw=1)
    text(left,startx,31,f'{s-2}–{s+5}',6.8,color=MUTED)
text(left,0,58,'X',7.4,va='center',fontweight='semibold');text(left,0,44,'Y',7.4,va='center',fontweight='semibold')
text(left,0,16,'ATP segments: 16/16 identical',7.6,fontweight='semibold',color=TEAL)
text(left,0,5,'2I4I AMP neighbors (≤4 Å): 14/14 identical',6.8,color=MUTED)

text(right,0,96,'RNA-binding surface',8.5,fontweight='semibold',color=PURPLE)
text(right,0,85,'FH-1321 / 16d and BA103 (modeled contacts)',7.0)
positions=[274,276,302,323,326,357,503]
for j,p in enumerate(positions):
    xx=18+j*12.5
    text(right,xx,73,str(p),6.7,ha='center',color=MUTED)
    for yv,g in [(58,'x_aa'),(44,'y_aa')]:
        box(right,xx-4,yv-5,8,10,PURPLE_LIGHT)
        text(right,xx,yv,xmap[p][g],8.0,ha='center',va='center',fontfamily='DejaVu Sans Mono')
    text(right,xx,31,str(xmap[p]['y_pos']),6.7,ha='center',color=MUTED)
text(right,0,58,'X',7.4,va='center',fontweight='semibold');text(right,0,44,'Y',7.4,va='center',fontweight='semibold')
text(right,0,16,'16d: 6/6 identical  |  BA103: 4/4 identical',7.4,fontweight='semibold',color=PURPLE)

# Preserve the original a–c physical sizing while adding a 47-mm bottom strip.
scale=145/192;offset=47/192
for a in fig.axes:
    p=a.get_position();a.set_position([p.x0,offset+scale*p.y0,p.width,scale*p.height])
for t in fig.texts:
    x,y=t.get_position();t.set_position((x,offset+scale*y))

panel(fig,.235,'d')
for i,w in enumerate(selected):
    a=fig.add_axes([.055+i*.48,.025,.435,.195]);a.set_xlim(0,100);a.set_ylim(0,100);a.axis('off')
    label='N-terminal IDR' if i==0 else 'RecA2 variable segment'
    text(a,0,96,label,8.5,fontweight='semibold',color=ORANGE)
    text(a,0,82,f"X {w['x_start']}–{w['x_end']}  /  Y {w['y_start']}–{w['y_end']}",7.2,color=MUTED)
    for j,p in enumerate(range(w['x_start'],w['x_end']+1)):
        xx=15+j*4.25;different=xmap[p]['status']!='identical'
        for yv,key in [(61,'x_aa'),(44,'y_aa')]:
            box(a,xx-1.92,yv-6,3.84,12,'#FBE8DF' if different else LIGHT)
            text(a,xx,yv,xmap[p][key],8.0,ha='center',va='center',fontfamily='DejaVu Sans Mono',color=ORANGE if different else INK)
    text(a,0,61,'X',7.4,va='center',fontweight='semibold')
    text(a,0,44,'Y',7.4,va='center',fontweight='semibold')
    text(a,0,24,f"{w['identical']}/20 identical ({w['identity_pct']:.0f}%)  |  {w['similarity_pct']:.0f}% similar",7.2,fontweight='semibold',color=ORANGE)
    text(a,0,7,'Disordered; binding untested' if i==0 else 'Resolved in 2I4I; binding untested',6.8,color=MUTED)

for ext in ['svg','pdf','png']:
    fig.savefig(OUT/f'ch04_main_ddx3x_ddx3y_inhibitor_site_conservation.{ext}',dpi=450)
plt.close(fig)
print('Exported four-panel Chapter 4 discussion figure as PDF, SVG and PNG.')
