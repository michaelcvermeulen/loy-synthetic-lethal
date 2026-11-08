"""A549 sampling workflow, using the editable cell/tube style of Fig 4.1a."""
import numpy as np
from matplotlib.path import Path
from matplotlib.patches import PathPatch, Circle, Ellipse, Arc, FancyArrowPatch, Rectangle

BLUE, RED, INK = '#3C5488', '#E64B35', '#222222'
GREY, LIGHT = '#777E86', '#D9DDE2'
GUIDES = ['#599E8B','#B79549','#8473AC','#66A5BB','#CB8B91','#8DA663','#BF8B5D','#728DAB']


def draw_schematic(fig, rect):
    ax=fig.add_axes(rect,xlim=(0,224),ylim=(0,64))
    ax.set_aspect('equal');ax.axis('off')
    def text(x,y,label,size=7.5,color=INK,weight='normal',ha='center'):
        return ax.text(x,y,label,fontsize=size,color=color,fontweight=weight,
                       ha=ha,va='center',linespacing=1.1)
    def arrow(x1,y1,x2,y2,color=GREY,width=.8):
        ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle='-|>',
                     mutation_scale=8,lw=width,color=color,shrinkA=0,shrinkB=0))
    def cell(x,y,r,colour,angle):
        a=np.linspace(0,2*np.pi,12,endpoint=False)+angle
        radii=r*np.array([1.28,.60,1.1,.68,1.31,.65,1.07,.67,1.26,.61,1.09,.64])
        points=np.c_[x+radii*np.cos(a),y+radii*np.sin(a)]
        mid=(points+np.roll(points,1,axis=0))/2
        vertices=[mid[0]];codes=[Path.MOVETO]
        for i in range(12):
            vertices.extend([points[i],mid[(i+1)%12]]);codes.extend([Path.CURVE3,Path.CURVE3])
        vertices.append(mid[0]);codes.append(Path.CLOSEPOLY)
        ax.add_patch(PathPatch(Path(vertices,codes),fc=colour,ec='white',lw=.35,zorder=3))
        ax.add_patch(Ellipse((x-.12*r,y+.05*r),.72*r,.51*r,angle=np.degrees(angle),
                             fc='white',alpha=.48,ec='none',zorder=4))
    def population(x,y,scale=.75,library=False,shift=0):
        offsets=[(-4.7,3.5),(-.2,4.8),(4.5,3),(-6,-1.1),(-1.1,.2),(4.4,-1.8),(-3.1,-4.8),(1.5,-5.1)]
        for i,(dx,dy) in enumerate(offsets):
            cell(x+dx*scale,y+dy*scale,2.85*scale,GUIDES[(i+shift)%8] if library else '#B9BFC5',i*.42)
    def tube(x,y,colour,scale=.68):
        vertices=[(x+dx*scale,y+dy*scale) for dx,dy in [(-1.9,3.3),(1.9,3.3),(1.9,-.5),(.8,-2.2),(-.8,-2.2),(-1.9,-.5),(-1.9,3.3)]]
        ax.add_patch(PathPatch(Path(vertices,[Path.MOVETO]+[Path.LINETO]*5+[Path.CLOSEPOLY]),fc='white',ec=colour,lw=.75,zorder=3))
        ax.add_patch(Rectangle((x-2.3*scale,y+3.3*scale),4.6*scale,.85*scale,fc=colour,ec='none',zorder=4))
        ax.add_patch(Ellipse((x,y-1.2*scale),1.9*scale,.9*scale,fc=colour,alpha=.45,ec='none',zorder=4))

    for x,label in [(14,'A549'),(40,'Brunello library'),(68,'Starting pools'),(100,'Outgrowth'),(137,'P3 collection'),(171,'P4 collection'),(207,'Guide readout')]:
        text(x,61,label,8.0,weight='semibold')
    text(40,57,'77,441 sgRNAs',7,color=GREY)
    text(100,57,'2 intended cultures / pool',6.8,color=GREY)
    for y,genotype,colour,shift in [(45,'CTRL',BLUE,0),(23,'LOY',RED,3)]:
        population(14,y)
        text(14,y-8,genotype,8,color=colour,weight='semibold')
        for dx,dy,i in [(-3.8,0,0),(0,1.4,2),(3.8,-.2,1)]:
            ax.add_patch(Circle((40+dx*.7,y+6+dy*.7),1.6,fill=False,ec=GREY,lw=.65))
            ax.add_patch(Arc((40+dx*.7,y+6+dy*.7),3.2,3.2,theta1=25,theta2=135,ec=GUIDES[i],lw=1.6))
        arrow(40,y+3.5,40,y+.7,width=.65)
        arrow(23,y,59,y)
        population(68,y,library=True,shift=shift)
        arrow(68,y-5.7,68,y-8.0,color=colour,width=.65)
        tube(65.8,y-12,colour,.60);tube(70.2,y-12,colour,.60)
        text(76,y-10.6,'D0',7.5,color=colour,weight='semibold',ha='left')
        text(76,y-14,'2 technical libraries',6.6,color=GREY,ha='left')
        ax.plot([76,82],[y,y],color=GREY,lw=.75)
        for delta,rep in [(4.2,'R1'),(-4.2,'R2')]:
            arrow(82,y,94,y+delta,width=.75)
            population(100,y+delta,scale=.39,library=True,shift=shift+(1 if rep=='R1' else 4))
            text(108,y+delta,rep,7,color=colour,weight='semibold')
            ax.plot([112,119],[y+delta,y+delta],color=LIGHT,lw=.75)
    # Aggregate sample collection avoids asserting a resolved lineage for the
    # P3 culture recorded as CTRL and reassigned to LOY in sensitivity analyses.
    ax.plot([119,119],[18.8,49.2],color=LIGHT,lw=.8)
    arrow(119,34,131,34)
    ax.plot([131,179],[34,34],color=GREY,lw=.8)
    arrow(141,34,165,34)
    for x in [137,171]:
        ax.add_patch(Circle((x,34),.9,fc=INK,ec='none',zorder=3))
        arrow(x,32.5,x,28.5,width=.65)
    text(137,42,'Earlier sample',7,color=GREY)
    text(171,42,'Endpoint',7,color=GREY)
    for x in [131,135,139,143]:tube(x,24,GREY,.63)
    for x,c in [(166.5,BLUE),(171,RED),(175.5,RED)]:tube(x,24,c,.63)
    text(137,17,'4 libraries',7.3)
    text(137,12.9,'PCR jackpotting',7,color=GREY)
    text(137,8.8,'Sensitivity only',7,color=GREY)
    text(171,17,'1 CTRL + 2 LOY',7.3)
    text(171,12.9,'Primary analysis',7.2,weight='semibold')
    text(171,8.8,'1 CTRL unavailable',6.8,color=GREY)
    arrow(179,34,187,34)
    text(207,53,'gDNA extraction',7.3)
    text(207,48.8,'sgRNA PCR + indexing',7.1)
    text(207,42.5,'Sequencing',7.7,weight='semibold')
    for i in range(5):
        yy=38.5-i*1.4
        ax.plot([194,220],[yy,yy],color=LIGHT,lw=1.2)
        xx=198+(i%3)*3
        ax.plot([xx,xx+8],[yy,yy],color=GUIDES[i],lw=1.9)
    text(207,28.5,'Guide abundance',7.2)
    arrow(207,25.8,207,23.8,width=.65)
    text(207,20.8,'drugZ',8.2,weight='semibold')
    text(207,16.5,'LOY vs CTRL at P4',7)
    text(207,12.8,'2 vs 1 cultures',6.8,color=GREY)
    text(4,2,'D0 technical libraries share one starting population per genotype; cell colours represent sgRNAs.',6.7,color=GREY,ha='left')
    return ax
