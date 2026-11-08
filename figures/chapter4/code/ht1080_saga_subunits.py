"""Editable subunit cartoon of hSAGA; positions are illustrative, not coordinates."""
import pandas as pd
from matplotlib.patches import Ellipse


MODULES = {
    'HAT': {'fill':'#D2E9ED', 'edge':'#79A9B1', 'hit':'#266A76'},
    'Core': {'fill':'#F3E8B9', 'edge':'#C6B768', 'hit':'#806414'},
    'DUB': {'fill':'#E7D5E9', 'edge':'#AC82B0', 'hit':'#74417B'},
    'Splicing': {'fill':'#D9E2E5', 'edge':'#8F9EA5', 'hit':'#53636B'},
    'TRRAP': {'fill':'#E2E2E2', 'edge':'#A3A3A3', 'hit':'#666666'},
}

# Relative placement follows the author's module/subunit layout. Membership uses
# the frozen Herbst 20-subunit reference, including TAF9B and ATXN7.
# gene, module, centre x/y, width/height in drawing coordinates
NODES = [
    ('KAT2A','HAT',.94,1.77,.76,.40),
    ('TADA2B','HAT',1.51,1.67,.88,.35),
    ('TADA3','HAT',1.00,1.50,.75,.34),
    ('SGF29','HAT',.43,1.45,.76,.35),
    ('SUPT7L','Core',2.10,1.49,.86,.34),
    ('TAF10','Core',2.04,1.27,.70,.32),
    ('TAF9B','Core',2.60,1.29,.74,.32),
    ('SUPT3H','Core',3.26,1.35,.87,.31),
    ('TAF6L','Core',2.25,1.02,.78,.32),
    ('TADA1','Core',2.86,1.02,.76,.32),
    ('TAF5L','Core',2.16,.76,.78,.33),
    ('TAF12','Core',2.84,.76,.72,.32),
    ('SUPT20H','Core',2.31,.50,.97,.34),
    ('ATXN7','DUB',1.51,.78,.77,.33),
    ('ATXN7L3','DUB',.66,.66,1.00,.35),
    ('USP22','DUB',1.23,.51,.73,.32),
    ('ENY2','DUB',.75,.38,.73,.35),
    ('TRRAP','TRRAP',3.20,.33,1.08,.60),
    ('SF3B3','Splicing',3.55,1.09,.81,.33),
    ('SF3B5','Splicing',3.48,.81,.81,.33),
]


def draw_saga_membership(ax, reference, data_dir):
    ax.set_xlim(0,4.8);ax.set_ylim(0,2.15);ax.set_axis_off()
    expected=reference.loc[reference.reference_20,['gene','module']].set_index('gene').module.to_dict()
    fdr=reference.set_index('gene').fdr_synth
    nominated=set(fdr.index[fdr.lt(.25)])
    strong=set(fdr.index[fdr.lt(.05)])
    assert nominated=={'TAF5L','SUPT20H','TAF6L','TADA1','TADA2B','ATXN7L3'}
    records=[]
    # Shapes are drawn before their labels so overlaps never obscure text.
    for gene,module,x,y,w,h in NODES:
        assert expected[gene]==module,(gene,module)
        palette=MODULES[module];hit=gene in nominated
        ax.add_patch(Ellipse((x,y),w,h,facecolor=palette['hit'] if hit else palette['fill'],
                            edgecolor=palette['edge'],lw=.55,zorder=2))
        records.append({'gene':gene,'module':module,'x':x,'y':y,'width':w,'height':h,
                        'depletion_fdr':float(fdr.loc[gene]),'white_bold_fdr_below_025':hit,
                        'score_label_bold_fdr_below_005':gene in strong,
                        'source':'Herbst et al. 2021 Fig. 1a',
                        'doi':'10.1038/s41594-021-00682-7'})
    for gene,module,x,y,w,h in NODES:
        hit=gene in nominated
        ax.text(x,y+(.065 if gene=='TRRAP' else 0),gene,ha='center',va='center',fontsize=7.0,
                fontweight='semibold' if hit else 'normal',
                color='white' if hit else '#333333',zorder=5)
    ax.text(.07,2.075,'HAT',fontsize=8.0,fontweight='semibold',color=MODULES['HAT']['hit'],va='center')
    ax.text(.51,2.075,'H3K9/K14 acetylation',fontsize=6.6,color=MODULES['HAT']['hit'],va='center')
    ax.text(2.25,2.075,'White bold: nominated (FDR <0.25)',fontsize=6.5,color='#333333',va='center')
    ax.text(.08,1.15,'Core scaffold',fontsize=7.5,fontweight='semibold',color=MODULES['Core']['hit'],va='center')
    ax.plot([1.14,1.80],[1.15,1.10],color=MODULES['Core']['edge'],lw=.5,zorder=1)
    ax.text(.04,.16,'DUB',fontsize=7.5,fontweight='semibold',color=MODULES['DUB']['hit'],va='center')
    ax.text(.04,.025,'H2BK120 deubiquitination',fontsize=6.5,color=MODULES['DUB']['hit'],va='center')
    ax.text(4.00,1.00,'Splicing',fontsize=7.0,fontweight='semibold',color='#53636B',va='center')
    ax.text(3.20,-.065,'Activator binding',fontsize=6.7,color='#555555',ha='center',va='center')
    ax.text(3.20,.19,'Pseudokinase',fontsize=6.2,color='#555555',ha='center',va='center',zorder=5)
    assert len(records)==20 and {r['gene'] for r in records}==set(expected)
    pd.DataFrame(records).to_csv(data_dir/'saga_schematic_members.csv',index=False)
    return {'reference_genes':20,'modules':list(MODULES),
            'white_bold_genes':sorted(nominated),'white_bold_rule':'drugZ depletion FDR <0.25',
            'score_bold_genes':sorted(strong),'score_bold_rule':'drugZ depletion FDR <0.05',
            'description':'Editable subunit cartoon inspired by author reference; positions and contacts are illustrative',
            'source_doi':'10.1038/s41594-021-00682-7',
            'reference_name_choices':{'TAF9':'TAF9B','ATXN7L1/L2':'ATXN7'},
            'tbp_binding_not_asserted':True,'alternative_subunits_in_schematic':False}
