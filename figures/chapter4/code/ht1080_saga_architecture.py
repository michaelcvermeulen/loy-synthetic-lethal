"""Original editable hSAGA schematic; not an atomic model or density map."""
import numpy as np
import pandas as pd
from matplotlib.path import Path
from matplotlib.patches import FancyBboxPatch, PathPatch

MODULES = {
    'DUB': {'colour':'#8C5291','fill':'#F1E6F2','protein':'#AF72B4',
            'box':(.03,1.43,1.40,.52),'rows':[['ATXN7','ATXN7L3'],['ENY2','USP22']]},
    'HAT': {'colour':'#328A96','fill':'#E5F3F5','protein':'#62BBC7',
            'box':(2.94,1.43,1.43,.52),'rows':[['KAT2A','SGF29'],['TADA2B','TADA3']]},
    'Core': {'colour':'#947E2A','fill':'#FBF6DB','protein':'#E4CC59',
             'box':(2.80,.39,1.57,.84),
             'rows':[['SUPT3H','SUPT7L','SUPT20H'],['TADA1','TAF5L','TAF6L'],['TAF9B','TAF10','TAF12']]},
    'Splicing': {'colour':'#6D7D82','fill':'#EDF1F2','protein':'#BACACF',
                 'box':(.03,.88,1.31,.39),'rows':[['SF3B3','SF3B5']]},
    'TRRAP': {'colour':'#686868','fill':'#EEEEEE','protein':'#A3A3A3',
              'box':(.03,.27,1.31,.52),'rows':[['TRRAP']]},
}


def outline(cx, cy, width, height, phase=0, n=128):
    theta=np.linspace(0,2*np.pi,n,endpoint=False)
    radius=1+.035*np.sin(5*theta+phase)+.020*np.sin(9*theta+phase*.7)
    return np.c_[cx+width/2*radius*np.cos(theta),cy+height/2*radius*np.sin(theta)]


def closed_path(points):
    vertices=np.vstack([points,points[0]])
    codes=np.full(len(vertices),Path.LINETO);codes[0]=Path.MOVETO;codes[-1]=Path.CLOSEPOLY
    return Path(vertices,codes)


def lobe(ax, module, bounds, phase=0, z=3):
    spec=MODULES[module]
    ax.add_patch(PathPatch(closed_path(outline(*bounds,phase)),facecolor=spec['protein'],
                          edgecolor=spec['colour'],lw=.50,zorder=z))


def draw_complex(ax):
    # Large TRRAP cradle beneath the core; all contours are illustrative.
    outer=outline(2.24,.73,1.08,.91,.6)
    inner=outline(2.30,.73,.49,.35,1.2)[::-1]
    ring=Path.make_compound_path(closed_path(outer),closed_path(inner))
    ax.add_patch(PathPatch(ring,facecolor=MODULES['TRRAP']['protein'],
                          edgecolor=MODULES['TRRAP']['colour'],lw=.55,zorder=2))
    for x,y,w,h,p in [(1.87,.94,.35,.44,.8),(2.04,1.08,.42,.30,.2)]:
        lobe(ax,'TRRAP',(x,y,w,h),p,3)
    for theta in np.linspace(-2.55,.5,6):
        r=np.array([.34,.49])
        ax.plot(2.24+r*np.cos(theta),.73+r*.82*np.sin(theta),
                color='#777777',lw=.38,alpha=.65,zorder=3)
    for bounds,phase in [((2.14,1.24,.56,.34),.7),((2.50,1.25,.48,.30),1.4),
                         ((1.95,1.31,.30,.39),2.2)]:
        lobe(ax,'Core',bounds,phase,5)
    ax.plot([1.72,1.96],[1.08,1.22],color=MODULES['Splicing']['colour'],lw=2.8,zorder=1)
    for bounds,phase in [((1.59,1.09,.30,.25),1),((1.74,.91,.26,.29),2),
                         ((1.81,1.16,.22,.22),3)]:
        lobe(ax,'Splicing',bounds,phase,4)
    # HAT and DUB were unresolved in Herbst et al.; dashed lines denote flexibility.
    ax.plot([1.89,1.98],[1.55,1.39],color=MODULES['DUB']['colour'],lw=.8,ls=(0,(2,2)),zorder=4)
    ax.plot([2.27,2.30],[1.48,1.37],color=MODULES['HAT']['colour'],lw=.8,ls=(0,(2,2)),zorder=4)
    lobe(ax,'DUB',(1.89,1.76,.37,.39),.3,5)
    lobe(ax,'DUB',(1.75,1.76,.18,.28),1.7,4)
    lobe(ax,'HAT',(2.25,1.73,.35,.49),1.5,5)
    lobe(ax,'HAT',(2.14,1.61,.21,.30),.2,6)
    # Activator on DNA, providing the TRRAP functional context.
    x=np.linspace(1.46,2.78,100);y=.12+.014*np.sin((x-1.46)*15)
    ax.plot(x,y,color='#B8843F',lw=1.15,zorder=1)
    ax.plot(x,y-.045,color='#B8843F',lw=.65,zorder=1)
    for xx in np.linspace(1.48,2.76,17):
        yy=.12+.014*np.sin((xx-1.46)*15)
        ax.plot([xx,xx],[yy-.045,yy],color='#B8843F',lw=.4,zorder=1)
    ax.add_patch(PathPatch(closed_path(outline(1.97,.23,.23,.26,.8)),
                          facecolor='#EAC99A',edgecolor='#A77738',lw=.5,zorder=4))
    ax.plot([2.03,2.09],[.35,.40],color='#8B744F',lw=.65,zorder=2)
    ax.text(2.86,.10,'DNA',fontsize=6.6,color='#6E5C42',va='center')
    ax.text(1.12,.16,'Activator',fontsize=6.5,color='#6E5C42',ha='center')
    ax.plot([1.38,1.84],[.19,.24],color='#9B8A72',lw=.5,zorder=1)


def draw_saga_membership(ax, reference, data_dir):
    ax.set_xlim(0,4.4);ax.set_ylim(0,2.15);ax.set_axis_off()
    expected=reference.loc[reference.reference_20,['gene','module']].set_index('gene').module.to_dict()
    fdr=reference.set_index('gene').fdr_synth
    significant=set(fdr.index[fdr.lt(.05)])
    draw_complex(ax)
    # Neutral leaders identify modules, not molecular contacts.
    for xs,ys in [([1.43,1.72],[1.69,1.77]),([2.94,2.43],[1.69,1.77]),
                  ([1.34,1.48],[1.08,1.09]),([1.34,1.78],[.50,.64]),
                  ([2.80,2.67],[.98,1.16])]:
        ax.plot(xs,ys,color='#929292',lw=.55,zorder=1)
    ax.text(.73,2.065,'H2BK120 deubiquitination',ha='center',va='center',
            fontsize=7.1,color=MODULES['DUB']['colour'])
    ax.text(3.66,2.065,'H3K9/K14 acetylation',ha='center',va='center',
            fontsize=7.1,color=MODULES['HAT']['colour'])
    records=[]
    for module,spec in MODULES.items():
        x,y,w,h=spec['box']
        ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.012,rounding_size=0.035',
            facecolor=spec['fill'],edgecolor=spec['colour'],linewidth=.55,zorder=8))
        title='Splicing (SPL)' if module=='Splicing' else module
        ax.text(x+w/2,y+h-.085,title,ha='center',va='center',fontsize=7.3,
                fontweight='normal' if module=='TRRAP' else 'semibold',color=spec['colour'],zorder=9)
        for row,genes in enumerate(spec['rows']):
            gy=y+h-.245-row*(.185 if module=='Core' else .155)
            if module=='Splicing': gy=y+.105
            if module=='TRRAP': gy=y+h-.085
            for col,gene in enumerate(genes):
                gx=x+w*(col+.5)/len(genes)
                assert expected[gene]==module,(gene,module)
                if module!='TRRAP':
                    ax.text(gx,gy,gene,ha='center',va='center',fontsize=6.9,color='#333333',
                            fontweight='semibold' if gene in significant else 'normal',zorder=9)
                records.append({'gene':gene,'module':module,'x':gx,'y':gy,
                    'depletion_fdr':float(fdr.loc[gene]),'bold_fdr_below_005':gene in significant,
                    'module_colour':spec['colour'],'source':'Herbst et al. 2021 Fig. 1a',
                    'doi':'10.1038/s41594-021-00682-7'})
        if module=='TRRAP':
            ax.text(x+w/2,y+.27,'Pseudokinase (ΨPIKK)',ha='center',va='center',
                    fontsize=6.4,color='#555555',zorder=9)
            ax.text(x+w/2,y+.095,'Binds c-MYC, E2F1, p53',ha='center',va='center',
                    fontsize=6.4,color='#555555',zorder=9)
        if module=='Core':
            ax.text(x+w/2,y+.065,'Histone-fold scaffold',ha='center',va='center',
                    fontsize=6.4,color=spec['colour'],zorder=9)
    assert len(records)==20 and {r['gene'] for r in records}==set(expected)
    pd.DataFrame(records).to_csv(data_dir/'saga_schematic_members.csv',index=False)
    return {'reference_genes':20,'modules':list(MODULES),
            'bold_genes':sorted(significant),'bold_rule':'drugZ depletion FDR <0.05',
            'description':'Original vector silhouette of illustrative hSAGA architecture and functions; not an atomic model or density map',
            'flexible_modules':['HAT','DUB'],'positions_and_subunit_boundaries_not_experimental':True,
            'functional_labels':['H2BK120 deubiquitination','H3K9/K14 acetylation',
                                 'Histone-fold scaffold','Pseudokinase (ΨPIKK)','Binds c-MYC, E2F1, p53'],
            'activator_source_doi':'10.1128/MCB.22.16.5650-5661.2002',
            'source_doi':'10.1038/s41594-021-00682-7',
            'tbp_in_reference_set':False,'alternative_subunits_in_schematic':False}
