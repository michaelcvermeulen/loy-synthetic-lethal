"""Standalone HT1080 pooled CRISPR-screen schematic with editable vector cells."""
from pathlib import Path
import hashlib
import json
import os
import sys

THESIS = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(THESIS.parent / '.figure_python_deps'))
OUT = THESIS / 'figures/chapter4'
DATA = OUT / 'source_data/ht1080_screen_schematic'
DATA.mkdir(parents=True, exist_ok=True)
os.environ.setdefault('MPLCONFIGDIR', str(DATA / 'mpl_cache'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.path import Path as MplPath
from matplotlib.patches import PathPatch, Circle, Ellipse, Arc, FancyArrowPatch, Rectangle
import pymupdf

sys.path.insert(0, str(THESIS / 'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style, condition_colours
colours = condition_colours()
BLUE, RED = colours['Y-retaining'], colours['Y-depleted']
INK, GREY, LIGHT = '#222222', '#777E86', '#D9DDE2'
GUIDES = ['#599E8B', '#B79549', '#8473AC', '#66A5BB', '#CB8B91', '#8DA663', '#BF8B5D', '#728DAB']
# Each compact population has one cell of each colour, with a distinct arrangement.
POOL_COLOUR_ORDERS = {
    'CTRL_selected': (0,3,5,7,1,6,4,2),
    'LOY_selected': (6,2,0,4,7,3,1,5),
    'CTRL_R1': (3,1,6,0,2,5,7,4),
    'CTRL_R2': (5,7,2,6,4,0,3,1),
    'LOY_R1': (1,4,7,2,5,6,0,3),
    'LOY_R2': (7,5,4,1,3,2,6,0),
}
STEM = 'ch04_main_ht1080_crispr_screen_schematic_draft'
W, H = 183, 92
def draw_schematic(fig, rect=(0,0,1,1), font_scale=1.0, include_note=True, compact=False):
    ax = fig.add_axes(rect, xlim=(0,W), ylim=(0,64 if compact else H))
    ax.set_aspect('equal')
    ax.axis('off')
    texts = []
    
    def text(x,y,label,size=7.5,color=INK,weight='normal',ha='center',va='center'):
        t=ax.text(x,y,label,fontsize=size*font_scale,color=color,fontweight=weight,
                  ha=ha,va=va,linespacing=1.05 if compact else 1.2)
        texts.append(t)
        return t
    
    def arrow(x1,y1,x2,y2,color=GREY,width=.8):
        ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle='-|>',
                     mutation_scale=8,lw=width,color=color,shrinkA=0,shrinkB=0))
    
    def cell(x,y,r,colour,angle=0):
        # Rounded, irregular adherent-cell silhouette; all geometry is schematic.
        n=12
        angles=np.linspace(0,2*np.pi,n,endpoint=False)+angle
        radii=r*np.array([1.28,.60,1.1,.68,1.31,.65,1.07,.67,1.26,.61,1.09,.64])
        points=np.column_stack([x+radii*np.cos(angles),y+radii*np.sin(angles)])
        mid=(points+np.roll(points,1,axis=0))/2
        vertices=[mid[0]]
        codes=[MplPath.MOVETO]
        for i in range(n):
            vertices.extend([points[i],mid[(i+1)%n]])
            codes.extend([MplPath.CURVE3,MplPath.CURVE3])
        vertices.append(mid[0]);codes.append(MplPath.CLOSEPOLY)
        ax.add_patch(PathPatch(MplPath(vertices,codes),fc=colour,ec='white',lw=.35,zorder=3))
        ax.add_patch(Ellipse((x-.12*r,y+.05*r),.72*r,.51*r,
                             angle=np.degrees(angle),fc='white',alpha=.48,ec='none',zorder=4))
    
    def population(x,y,scale=1,library=False,colour_order=None):
        offsets=[(-4.7,3.5),(-.2,4.8),(4.5,3.0),(-6.0,-1.1),(-1.1,.2),(4.4,-1.8),(-3.1,-4.8),(1.5,-5.1)]
        order=range(len(GUIDES)) if colour_order is None else colour_order
        for i,(dx,dy) in enumerate(offsets):
            cell(x+dx*scale,y+dy*scale,2.85*scale,
                 GUIDES[order[i]] if library else '#B9BFC5',angle=i*.42)
    
    def library(x,y,scale=1):
        for dx,dy,i in [(-3.8,0,0),(0,1.4,2),(3.8,-.2,1)]:
            ax.add_patch(Circle((x+dx*scale,y+dy*scale),2.25*scale,fill=False,ec=GREY,lw=.65))
            ax.add_patch(Arc((x+dx*scale,y+dy*scale),4.5*scale,4.5*scale,theta1=25,theta2=135,ec=GUIDES[i],lw=1.7))
    
    def tube(x,y,colour,scale=1):
        verts=[(x+dx*scale,y+dy*scale) for dx,dy in
               [(-1.9,3.3),(1.9,3.3),(1.9,-.5),(.8,-2.2),(-.8,-2.2),(-1.9,-.5),(-1.9,3.3)]]
        ax.add_patch(PathPatch(MplPath(verts,[MplPath.MOVETO]+[MplPath.LINETO]*5+[MplPath.CLOSEPOLY]),
                               fc='white',ec=colour,lw=.75,zorder=3))
        ax.add_patch(Rectangle((x-2.3*scale,y+3.3*scale),4.6*scale,.85*scale,fc=colour,ec='none',zorder=4))
        ax.add_patch(Ellipse((x,y-1.2*scale),1.9*scale,.9*scale,fc=colour,alpha=.45,ec='none',zorder=4))

    if compact:
        # Shorter two-row layout leaves the lower-left area for molecular confirmation.
        for x,label in [(15,'HT1080'),(43,'Brunello library'),(73,'Selected pools'),
                        (112,'Replicate outgrowth'),(161,'Guide readout')]:
            text(x,61,label,8.1,weight='semibold')
        text(43,58.2,'77,441 sgRNAs',7,color=GREY)
        text(112,58.2,'2 cultures / population',7,color=GREY)
        for row_y,genotype,colour in [(48,'CTRL',BLUE),(33,'LOY',RED)]:
            population(15,row_y,scale=.67)
            text(2,row_y,genotype,8,color=colour,weight='semibold',ha='left')
            library(43,row_y+4.8,scale=.8)
            arrow(43,row_y+2.5,43,row_y+.7,width=.65)
            arrow(23,row_y,63,row_y)
            if genotype == 'CTRL':
                text(43,row_y-3.3,'Lentiviral transduction',6.4,color=GREY)
                text(43,row_y-5.8,'Puromycin selection',6.4,color=GREY)
            population(72,row_y,scale=.70,library=True,
                       colour_order=POOL_COLOUR_ORDERS[f'{genotype}_selected'])
            arrow(75,row_y-4.7,80,row_y-5.8,color=colour,width=.65)
            tube(81,row_y-8.3,colour,scale=.6)
            text(84,row_y-8,'D0',7,color=colour,weight='semibold',ha='left')
            text(88.4,row_y-8,'Baseline QC',6.4,color=GREY,ha='left')
            ax.plot([79,87],[row_y,row_y],color=GREY,lw=.75)
            for delta,rep in [(3.75,'R1'),(-3.75,'R2')]:
                arrow(87,row_y,105,row_y+delta,width=.75)
                population(111,row_y+delta,scale=.40,library=True,
                           colour_order=POOL_COLOUR_ORDERS[f'{genotype}_{rep}'])
                text(119,row_y+delta,rep,7,color=colour,weight='semibold')
                ax.plot([122,132],[row_y+delta,row_y+delta],color=LIGHT,lw=.75)
        ax.plot([132,132],[29.25,51.75],color=LIGHT,lw=.8)
        arrow(132,40.5,141,40.5)
        text(113,20.7,'4 passages',7.8,weight='semibold')
        text(113,18,'9.3–9.8 population doublings',6.8,color=GREY)
        text(113,15.3,'32 million cells / culture / passage',6.5,color=GREY)
        text(161,58.2,'D0 + endpoint samples',7.1,color=GREY)
        text(161,54.8,'gDNA extraction',7.5)
        text(161,52,'sgRNA PCR + indexing',7.2)
        arrow(161,49.8,161,48.2,width=.65)
        text(161,46,'Sequencing',7.7,weight='semibold')
        for i in range(6):
            y=42.8-i*1.15
            ax.plot([148,173],[y,y],color=LIGHT,lw=1.2,solid_capstyle='butt')
            x=152+(i%3)*3
            ax.plot([x,x+8],[y,y],color=GUIDES[i],lw=1.9,solid_capstyle='butt')
        text(161,34.3,'sgRNA abundance',7.3)
        arrow(161,31.9,161,30.3,width=.65)
        text(161,28,'drugZ',8.2,weight='semibold')
        text(161,25,'LOY vs CTRL endpoints',7)
        text(161,22.2,'2 cultures per group',6.7,color=GREY)
        return ax,texts
    
    text(15,87,'HT1080',8.5,weight='semibold')
    text(43,87,'Brunello library',8.1,weight='semibold')
    text(43,82.5,'77,441 sgRNAs',7,color=GREY)
    text(73,87,'Selected pools',8.1,weight='semibold')
    text(112,87,'Replicate outgrowth',8.1,weight='semibold')
    text(112,82.5,'2 cultures / population',7,color=GREY)
    text(161,87,'Guide readout',8.1,weight='semibold')
    
    for row_y,genotype,colour in [(65,'CTRL',BLUE),(35,'LOY',RED)]:
        population(15,row_y,library=False)
        text(15,row_y-11.5,genotype,8,color=colour,weight='semibold')
        library(43,row_y+8.5)
        arrow(43,row_y+5.7,43,row_y+1.2,width=.65)
        arrow(27,row_y,59,row_y)
        text(43,row_y-4.6,'Lentiviral transduction',6.4,color=GREY)
        text(43,row_y-8.1,'Puromycin selection',6.4,color=GREY)
        population(72,row_y,library=True)
        # A D0 aliquot from each selected population precedes replicate outgrowth.
        arrow(72,row_y-8.3,72,row_y-12.2,color=colour,width=.65)
        tube(72,row_y-17,colour)
        text(78,row_y-16.1,'D0',7,color=colour,weight='semibold',ha='left')
        text(78,row_y-20,'Baseline QC',6.4,color=GREY,ha='left')
        # R1/R2 are outgrowth cultures from the same selected population.
        ax.plot([82,87],[row_y,row_y],color=GREY,lw=.75)
        for delta,rep in [(7,'R1'),(-7,'R2')]:
            arrow(87,row_y,102,row_y+delta,width=.75)
            population(111,row_y+delta,scale=.58,library=True)
            text(121,row_y+delta,rep,7,color=colour,weight='semibold')
            ax.plot([124,132],[row_y+delta,row_y+delta],color=LIGHT,lw=.75)
    
    ax.plot([132,132],[28,72],color=LIGHT,lw=.8)
    arrow(132,50,141,50)
    text(113,17,'4 passages',7.8,weight='semibold')
    text(113,12.5,'9.3–9.8 population doublings',6.8,color=GREY)
    text(113,8.2,'32 million cells / culture / passage',6.5,color=GREY)
    
    # All six samples are prepared separately; only the four endpoints enter drugZ.
    text(161,78,'D0 + endpoint samples',7.1,color=GREY)
    text(161,72,'gDNA extraction',7.5)
    text(161,67.5,'sgRNA PCR + indexing',7.2)
    arrow(161,64,161,62,width=.65)
    text(161,58.8,'Sequencing',7.7,weight='semibold')
    for i in range(6):
        y=52.2-i*2.25
        ax.plot([148,173],[y,y],color=LIGHT,lw=1.2,solid_capstyle='butt')
        x=152+(i%3)*3
        ax.plot([x,x+8],[y,y],color=GUIDES[i],lw=1.9,solid_capstyle='butt')
    text(161,36,'sgRNA abundance',7.3)
    arrow(161,32.5,161,28.5)
    text(161,24.5,'drugZ',8.2,weight='semibold')
    text(161,18,'LOY vs CTRL endpoints',7)
    text(161,13.5,'2 cultures per group',6.7,color=GREY)
    
    if include_note:
        text(4,4,'Cell colours indicate different sgRNAs; cell numbers and composition are schematic.',6.5,color=GREY,ha='left')
    return ax, texts

def main():
    apply_style(8)
    fig = plt.figure(figsize=(W/25.4, H/25.4))
    ax, texts = draw_schematic(fig)
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer()
    outside=[]
    for t in texts:
        box=t.get_window_extent(renderer)
        if not fig.bbox.contains(box.x0,box.y0) or not fig.bbox.contains(box.x1,box.y1):
            outside.append(t.get_text())
    assert not outside,outside
    for extension in ['pdf','svg','png']:
        fig.savefig(OUT/f'{STEM}.{extension}',dpi=600)
    review=OUT/'review'
    review.mkdir(exist_ok=True)
    fig.savefig(review/f'{STEM}_preview.png',dpi=200)
    plt.close(fig)
    doc=pymupdf.open(OUT/f'{STEM}.pdf')
    page=doc[0]
    page.get_pixmap(matrix=pymupdf.Matrix(2,2),colorspace=pymupdf.csGRAY).save(review/f'{STEM}_greyscale.png')
    sources=[THESIS/'draftv1/chapter4/chapter4_methods.md',THESIS/'draftv1/chapter4/chapter4_results.md']
    record={
        'figure_mm':[W,H], 'all_artwork_vector':True, 'svg_editable_text':'<text' in (OUT/f'{STEM}.svg').read_text(encoding='utf-8'),
        'font_names':[f[3] for f in page.get_fonts()], 'text_outside_canvas':outside,
        'screen':'HT1080 genome-wide Brunello knockout screen',
        'library_guides':77441, 'non_targeting_guides':1000,
        'baseline_samples':2, 'outgrowth_cultures_per_population':2, 'endpoint_samples':4,
        'passages':4, 'cumulative_population_doublings':[9.3,9.8],
        'cells_replated_per_culture_per_passage':32000000,
        'primary_analysis':'drugZ: two LOY endpoints versus two CTRL endpoints; D0 excluded',
        'baseline_role':'Sampled before splitting into outgrowth cultures; secondary QC.',
        'schematic_cell_colours':'Distinct sgRNAs, not fluorescent reporters or quantitative frequencies.',
        'replication':'Two outgrowth cultures derived from one selected population per genotype.',
        'omitted_unresolved_details':['Cas9 delivery route','Library backbone/catalogue confirmation'],
        'reference':'user_layout_reference.png used for visual inspiration; no ES2, metabolic-library or drug-treatment details transferred.',
        'sources':[{'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in sources],
    }
    (DATA/'schematic_metadata.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(record))


if __name__ == '__main__':
    main()
