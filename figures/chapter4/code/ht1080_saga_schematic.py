"""Editable module-membership schematic from the frozen Herbst reference set."""
from matplotlib.patches import FancyBboxPatch
import pandas as pd

MODULES = {
    'DUB': {'colour':'#8C5291','fill':'#F1E6F2','box':(.04,.96,1.70,.54),
            'rows':[['ATXN7','ATXN7L3'],['ENY2','USP22']]},
    'HAT': {'colour':'#328A96','fill':'#E5F3F5','box':(2.63,.96,1.70,.54),
            'rows':[['KAT2A','SGF29'],['TADA2B','TADA3']]},
    'Core': {'colour':'#947E2A','fill':'#FBF6DB','box':(1.43,.05,2.90,.75),
             'rows':[['SUPT3H','SUPT7L','SUPT20H'],['TADA1','TAF5L','TAF6L'],['TAF9B','TAF10','TAF12']]},
    'Splicing': {'colour':'#6D7D82','fill':'#EDF1F2','box':(.04,.49,1.19,.33),
                 'rows':[['SF3B3','SF3B5']]},
    'TRRAP': {'colour':'#686868','fill':'#EEEEEE','box':(.04,.05,1.19,.33),
              'rows':[['TRRAP']]},
}

def draw_saga_membership(ax, reference, data_dir):
    """Module links indicate grouping only, not residue contacts or geometry."""
    ax.set_xlim(0,4.4); ax.set_ylim(0,1.55); ax.set_axis_off()
    expected=reference.loc[reference.reference_20,['gene','module']].set_index('gene').module.to_dict()
    fdr=reference.set_index('gene').fdr_synth
    significant=set(fdr.index[fdr.lt(.05)])
    records=[]
    # Abstract connections between the labelled modules and core.
    for xs,ys in [([1.2,1.9],[.96,.8]),([3.4,3.0],[.96,.8]),
                  ([1.23,1.43],[.655,.57]),([1.23,1.43],[.215,.27])]:
        ax.plot(xs,ys,color='#C4C4C4',lw=.8,zorder=0)
    for module,spec in MODULES.items():
        x,y,w,h=spec['box']
        ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.012,rounding_size=0.035',
            facecolor=spec['fill'],edgecolor=spec['colour'],linewidth=.55))
        ax.text(x+w/2,y+h-.085,module,ha='center',va='center',fontsize=8,
                fontweight='semibold',color=spec['colour'])
        nrows=len(spec['rows'])
        for row,genes in enumerate(spec['rows']):
            gy=y+h-.25-row*.165
            if nrows==1: gy=y+.085
            for col,gene in enumerate(genes):
                gx=x+w*(col+.5)/len(genes)
                assert expected[gene]==module,(gene,module)
                ax.text(gx,gy,gene,ha='center',va='center',fontsize=7.5,color='#333333',
                        fontweight='semibold' if gene in significant else 'normal')
                records.append({'gene':gene,'module':module,'x':gx,'y':gy,
                    'depletion_fdr':float(fdr.loc[gene]),'bold_fdr_below_005':gene in significant,
                    'module_colour':spec['colour'],'source':'Herbst et al. 2021 Fig. 1a',
                    'doi':'10.1038/s41594-021-00682-7'})
    assert len(records)==20 and {r['gene'] for r in records}==set(expected)
    pd.DataFrame(records).to_csv(data_dir/'saga_schematic_members.csv',index=False)
    return {'reference_genes':20,'modules':list(MODULES),
            'bold_genes':sorted(significant),'bold_rule':'drugZ depletion FDR <0.05',
            'description':'Module membership schematic; positions and links are illustrative',
            'tbp_in_reference_set':False,'alternative_subunits_in_schematic':False}
