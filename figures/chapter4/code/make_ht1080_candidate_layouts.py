"""Guide distributions and matched DepMap/STRING networks in the main figure."""
from pathlib import Path
import json, hashlib
import make_ht1080_candidates as C
import make_ht1080_depmap_string_network as N
import ht1080_main_guide_distributions as G
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import NullLocator

def make_main():
    C.apply_style(9)
    fig=plt.figure(figsize=(9.6,11.3))
    # Overall screen distribution stays separate from the named-hit display.
    C.letter(fig,'a',.018,.985)
    ax=fig.add_axes([.085,.825,.225,.128])
    for cl,col,size in [('FDR ≥ 0.25',C.NS,3),('0.05–0.25',C.EXPLORATORY,10),('FDR < 0.05',C.STRONG,14)]:
        rows=C.DZ.loc[C.DZ.display_class.eq(cl)]
        ax.scatter(rows.rank_synth,rows.normZ,s=size,color=col,lw=0,
            rasterized=cl=='FDR ≥ 0.25',zorder=2 if cl=='FDR ≥ 0.25' else 3)
    ax.axhline(0,color='#888888',lw=.6,zorder=0)
    ax.set(xscale='log',xlim=(.8,25000),ylim=(-12.5,8),yticks=[-12,-8,-4,0,4])
    ax.set_xticks([1,100,10000],['1','100','10,000']);ax.xaxis.set_minor_locator(NullLocator())
    ax.tick_params(labelsize=7.5,length=2.5,pad=2)
    ax.set_xlabel('Depletion rank (log scale)',fontsize=8,labelpad=3)
    ax.set_ylabel('drugZ normZ',fontsize=8,labelpad=3)
    ax.text(.04,.95,'19,112 genes',transform=ax.transAxes,fontsize=7.5,va='top')
    label_positions={'DDX3X':(1.8,-11.1),'TIMM8A':(4.6,-9.0),
                     'SREBF2':(16,-6.8),'TAF6L':(120,-5.1),
                     'MED12':(1.2,-3.0),'SREBF1':(8,1.3),'TADA2B':(55,-1.3)}
    for gene,position in label_positions.items():
        point=C.INDEXED.loc[gene]
        ax.annotate(gene,(point.rank_synth,point.normZ),xytext=position,
                    fontsize=6.4,color='#858585',ha='left',va='center',
                    arrowprops={'arrowstyle':'-','color':'#B0B0B0','lw':.45,'shrinkA':1,'shrinkB':2})
    fig.text(.085,.770,'FDR <0.25: 48 depleted; 9 enriched',fontsize=7.4)

    # Keep the four prespecified gametologs, including nonsignificant results.
    C.letter(fig,'b',.018,.755)
    gametologs=['DDX3X']+C.OTHER

    # Same label rule the author chose: 28 depleted + TAF5L + 2 enriched.
    C.letter(fig,'c',.365,.985)
    saga=pd.read_csv(C.DATA/'saga_reference_scores.csv')
    saga_hits=set(saga.loc[saga.reference_20 & saga.fdr_synth.lt(.25),'source_gene'])
    depleted=C.DZ.loc[C.DZ.normZ.lt(0)&(C.DZ.fdr_synth.lt(.05)|C.DZ.GENE.isin(saga_hits))].sort_values('rank_synth')
    enriched=C.DZ.loc[C.DZ.normZ.gt(0)&C.DZ.fdr_supp.lt(.05)].sort_values('rank_supp')
    named=pd.concat([depleted,enriched],ignore_index=True)
    assert len(depleted)==29 and len(enriched)==2 and len(named)==31
    named['direction']=np.where(named.normZ.lt(0),'depleted','enriched')
    named['SAGA_candidate']=named.GENE.isin(saga_hits)
    named['display_column']=np.where(np.arange(len(named))<16,1,2)
    named['display_row']=np.arange(len(named))%16
    named.to_csv(C.DATA/'main_named_hits.csv',index=False)
    chromosome_source=C.OUT/'source_data/ht1080_chromosome_guide_bias/all_guide_chromosome_values.csv'
    chromosomes=pd.read_csv(chromosome_source,usecols=['count_gene','chromosome'],dtype=str)
    x_linked=set(chromosomes.loc[chromosomes.chromosome.eq('X'),'count_gene']) & set(named.GENE)
    assert x_linked=={'DDX3X','TIMM8A','MED12','BCORL1','ACSL4'}
    x_label_colour='#0072B2'
    selected,grid,curve,guide_manifest=G.prepare(named,gametologs)
    marks=[]
    bx=fig.add_axes([.14,.675,.17,.084])
    G.draw(bx,gametologs,selected,grid,curve,'b',1,marks,font_size=7.5)
    for col,start in enumerate([0,16]):
        part=named.iloc[start:start+16]
        cx=fig.add_axes([.445+col*.33,.655,.165,.295])
        G.draw(cx,part.GENE.tolist(),selected,grid,curve,'c',col+1,marks,saga=saga_hits,slots=16)
        for label in cx.texts:
            if label.get_text().rstrip('*') in x_linked:
                label.set_color(x_label_colour)
    G.finish(marks,guide_manifest)
    C.fdr_legend(fig,(.085,.624),size=7.8)
    fig.legend(handles=[Line2D([],[],color='#333333',marker='|',ms=6,lw=0,label='Guide'),
                        Line2D([],[],color='#AAAAAA',lw=2,label='Non-targeting reference (n = 999)')],
               loc='upper left',bbox_to_anchor=(.53,.620),ncol=2,fontsize=7.1,
               handletextpad=.35,columnspacing=1,borderaxespad=0,handlelength=1.3)

    # Preserve a-c at their reviewed physical size; reuse frozen network positions.
    height=9.4
    recovered_height=.5
    C.resize_height_from_bottom(fig,height)
    # Explain label encodings beside the guide key, above the network panel.
    label_key=[]
    key_x=.365
    for wording,colour,weight in [
        ('Blue labels: X-linked',x_label_colour,'normal'),
        ('Bold: SAGA candidates',C.DARK,'semibold'),
        ('*: 1 of 4 guides undetected in all six libraries (b,c)',C.DARK,'normal')]:
        item=fig.text(key_x,4.65/height,wording,fontsize=7.1,color=colour,fontweight=weight,va='baseline')
        fig.canvas.draw()
        key_x=item.get_window_extent(fig.canvas.get_renderer()).x1/fig.bbox.width+.012
        label_key.append(item)
    assert key_x<=.985, f'Label key exceeds page width: {key_x}'
    groups,module,core,strip,order,pos=N.load_saved_structure()
    C.letter(fig,'d',.018,(5.00-recovered_height)/height)
    fig.text(.075,(4.95-recovered_height)/height,'STRING: score ≥0.4 · 82 associations',fontsize=8)
    fig.text(.555,(4.95-recovered_height)/height,'DepMap: |r| > 0.2 · 41 connections',fontsize=8)
    # Crop unused data-coordinate margins while preserving physical node spacing.
    network_height=3.83/1.20*1.02
    left,ll=N.network_panel(fig,[.025,(1.50-recovered_height)/height,.455,network_height/height],
                           'string',groups,core,pos,fontsize=7.2)
    right,rl=N.network_panel(fig,[.520,(1.50-recovered_height)/height,.455,network_height/height],
                            'depmap',groups,core,pos,fontsize=7.2)
    left.set_ylim(0,1.02);right.set_ylim(0,1.02)
    fig.text(.04,.85/height,'Shared candidates with no displayed edge (18)',fontsize=7.1,color='#666666')
    fig.text(.975,.85/height,'* ZCCHC16: DepMap unavailable',ha='right',fontsize=6.8,color='#666666')
    shared=fig.add_axes([.025,.29/height,.95,.50/height])
    shared.set(xlim=(0,1),ylim=(0,1));shared.axis('off')
    shared.axhline(1,color='#DDDDDD',lw=.65)
    shared_labels=[]
    shared_positions={}
    for i,gene in enumerate(strip):
        xy=((i%9+.5)/9,.78 if i<9 else .25)
        shared_positions[gene]=xy
        # The shared strip encodes screen direction; missing DepMap data are
        # retained as a label note, independently of that direction encoding.
        label=N.node(shared,gene,xy,'string',strip=True,fontsize=7)
        if not N.NODES.loc[gene,'depmap_available']:
            label.set_text(gene+'*')
        shared_labels.append(label)
    assert len(left.lines)==82 and len(right.lines)==41
    assert len(left.collections)==len(right.collections)==len(core)==39
    for i,gene in enumerate(core):
        assert np.array_equal(left.collections[i].get_offsets(),right.collections[i].get_offsets())
        assert np.array_equal(left.collections[i].get_sizes(),right.collections[i].get_sizes())
        assert np.isclose(left.collections[i].get_sizes()[0],N.AREA_FACTOR*abs(N.NODES.loc[gene,'normZ']))
    assert len(shared.collections)==len(strip)==18 and len(set(core)|set(strip))==57
    for i,gene in enumerate(strip):
        assert np.isclose(shared.collections[i].get_sizes()[0],N.AREA_FACTOR*abs(N.NODES.loc[gene,'normZ']))
    pd.DataFrame([{'GENE':gene,'row':i//9+1,'column':i%9+1,
                   'x':shared_positions[gene][0],'y':shared_positions[gene][1],
                   'node_area_pt2':N.AREA_FACTOR*abs(N.NODES.loc[gene,'normZ']),
                   'depmap_available':bool(N.NODES.loc[gene,'depmap_available'])}
                  for i,gene in enumerate(strip)]).to_csv(C.DATA/'main_shared_strip_layout.csv',index=False)
    handles=[Line2D([],[],marker='o',mfc='#EDE1EF',mec=N.DEPLETED,lw=0,ms=6,label='LOY-depleted (48)'),
             Line2D([],[],marker='o',mfc='white',mec=N.ENRICHED,lw=0,ms=6,label='LOY-enriched (9)'),
             Line2D([],[],color=N.POSITIVE,lw=1.3,label='Positive r'),
             Line2D([],[],color=N.NEGATIVE,lw=1.3,ls='--',label='Negative r')]
    fig.legend(handles=handles,loc='center left',bbox_to_anchor=(.032,.14/height),ncol=4,
               fontsize=7.3,columnspacing=1,handletextpad=.4,handlelength=1.5,frameon=False)
    key=fig.add_axes([.72,.04/height,.25,.20/height]);key.axis('off');key.set(xlim=(0,1),ylim=(0,1))
    key.text(0,.5,'Area: |drugZ normZ|',fontsize=7,va='center')
    for x,z in zip([.52,.70,.88],[4,8,12]):
        key.scatter(x,.5,s=N.AREA_FACTOR*z,facecolor='white',edgecolor='#777777',lw=.8,clip_on=False)
        key.text(x+.07,.5,str(z),va='center',fontsize=7.2)
    # Explanations of edge width, shading and asterisks are in the figure caption.
    labels=ll+rl+shared_labels
    fig.canvas.draw();renderer=fig.canvas.get_renderer()
    assert not [(item.get_text(),i) for item in label_key for i,legend in enumerate(fig.legends)
                if item.get_window_extent(renderer).overlaps(legend.get_window_extent(renderer))]
    overlap=[[a.get_text(),b.get_text()] for i,a in enumerate(labels) for b in labels[i+1:]
             if a.get_window_extent(renderer).overlaps(b.get_window_extent(renderer))]
    assert not overlap,overlap
    C.export(fig,C.MAIN)
    record=C.RENDER_RECORD[C.MAIN]
    record.update({'layout':'labelled rank overview / gametolog guide distributions / named-hit guide distributions / matched DepMap and STRING',
        'panels':['a','b','c','d'],'abundance_panel_destination':'Supp Fig 2f',
        'named_hit_panel':'c','named_depleted':len(depleted),'named_enriched':len(enriched),
        'named_rule':'Directional FDR <0.05 plus SAGA depleted candidates at FDR <0.25',
        'additional_saga_labels':['TAF5L'],'network_panel':'d','candidate_set_per_view':57,
        'drawn_network_nodes_per_view':39,'shared_strip_drawn_once':True,
        'shared_strip_rows':2,'shared_strip_columns':9,'shared_strip_gene_order':strip,
        'height_recovered_inches':recovered_height-.3,'label_key_added_height_inches':.3,
        'shared_strip_missing_depmap_label':'ZCCHC16*',
        'network_core_nodes':39,'shared_strip_nodes':18,'depmap_edges':41,'string_edges':82,
        'shared_edges':14,'node_shape':'circle','node_area_pt2_per_abs_normZ':N.AREA_FACTOR,
        'node_colour':'screen direction','shared_positions':True,'layout_source':'frozen node_layout_and_modules.csv',
        'network_view_order':['STRING','DepMap'],'bottom_explanations':'central figure caption',
        'network_label_box_opacity':N.LABEL_BOX_ALPHA,
        'positive_depmap_communities':groups,'matrix_destination':'ch04_supp_ht1080_depmap_string_evidence_matrix_draft',
        'network_label_overlaps':overlap,
        'panel_a_labels':list(label_positions),'guide_panels':['b','c'],'guide_tick_height':'equal',
        'panel_c_x_linked_labels':sorted(x_linked),'panel_c_x_label_colour':x_label_colour,
        'on_figure_label_key':[item.get_text() for item in label_key],
        'chromosome_annotation_source':str(chromosome_source),
        'chromosome_annotation_sha256':hashlib.sha256(chromosome_source.read_bytes()).hexdigest(),
        'guide_reference_n':999,'guide_unique_genes':34,'guide_unique_observations':134,
        'guide_marks_including_repeated_DDX3X':137,
        'network_input_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest()
            for p in [N.DATA/'candidates.csv',N.DATA/'pairwise_evidence.csv',N.DATA/'node_layout_and_modules.csv']}})
    return record

if __name__=='__main__':
    record=make_main()
    manifest=C.DATA/'render_manifest.json'
    existing=json.loads(manifest.read_text(encoding='utf-8')) if manifest.exists() else {}
    existing[C.MAIN]=record
    manifest.write_text(json.dumps(existing,indent=2),encoding='utf-8')
    print(json.dumps(record,indent=2))
