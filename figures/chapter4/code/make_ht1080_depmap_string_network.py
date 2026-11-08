"""Matched DepMap/STRING networks and complete evidence matrix for 57 hits."""
from pathlib import Path
import sys, os, json, hashlib, itertools
THESIS = Path(__file__).resolve().parents[3]
sys.path.insert(0,str(THESIS.parent/'.figure_python_deps'))
OUT=THESIS/'figures/chapter4'
DATA=OUT/'source_data/ht1080_depmap_string_network_57'
os.environ.setdefault('MPLCONFIGDIR',str(DATA/'mpl_cache'))
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Polygon, Rectangle
import pymupdf as fitz
sys.path.insert(0,str(THESIS/'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
apply_style(8)

STEM='ch04_supp_ht1080_depmap_string_network_draft'
W,H=11.2,14.4
DEPLETED,ENRICHED='#7A3E85','#269E95'
POSITIVE,NEGATIVE='#629EBB','#C97769'
AREA_FACTOR=15.0  # Exact plotted circle area (points squared) / |normZ|.
LABEL_BOX_ALPHA=.35  # White label background: 65% transparent.
NODES=pd.read_csv(DATA/'candidates.csv').set_index('GENE',drop=False)
PAIRS=pd.read_csv(DATA/'pairwise_evidence.csv')
DE=PAIRS.loc[PAIRS.depmap_edge].copy()
ST=PAIRS.loc[PAIRS.string_edge].copy()
MODULE_COLS=['#F2ECF4','#EEF4F4','#F8F0E9']


def build_structure():
    g=nx.Graph()
    g.add_nodes_from(NODES.index)
    for e in DE.loc[DE.pearson_r.gt(0)].itertuples():
        g.add_edge(e.gene_a,e.gene_b,weight=e.pearson_r)
    communities=nx.community.louvain_communities(g,weight='weight',resolution=1,seed=17)
    groups=sorted([sorted(s,key=lambda x:-abs(NODES.loc[x,'normZ'])) for s in communities if len(s)>1],
                  key=lambda s:(-len(s),-max(abs(NODES.loc[x,'normZ']) for x in s)))
    assert [len(s) for s in groups]==[8,6,3]
    module={gene:i+1 for i,s in enumerate(groups) for gene in s}
    core=NODES.loc[NODES.depmap_degree.add(NODES.string_degree).gt(0),'GENE'].tolist()
    strip=NODES.loc[NODES.depmap_degree.add(NODES.string_degree).eq(0),'GENE'].tolist()
    assert len(core)==39 and len(strip)==18
    # Shared coordinates reflect both edge sets; community membership is derived
    # only from positive DepMap edges. Negative correlations do not enter Louvain.
    layout=nx.Graph();layout.add_nodes_from(core)
    for e in PAIRS.loc[PAIRS.depmap_edge | PAIRS.string_edge].itertuples():
        weight=(4*e.pearson_r if e.depmap_edge and e.pearson_r>0 else .15 if e.depmap_edge else 0)
        weight+=.35*e.string_score if e.string_edge else 0
        layout.add_edge(e.gene_a,e.gene_b,weight=weight)
    initial=nx.spring_layout(layout,seed=17,k=.65,iterations=700,weight='weight')
    a=np.array([initial[x] for x in core])
    a=(a-a.min(axis=0))/(a.max(axis=0)-a.min(axis=0))
    # Modest geometric repulsion reduces labels colliding at the actual panel size.
    # It changes positions only, never node membership, edge weights or modules.
    rng=np.random.default_rng(17)
    for _ in range(2200):
        delta=np.zeros_like(a)
        for i,j in itertools.combinations(range(len(core)),2):
            d=a[i]-a[j]
            need_x=.135 if max(len(core[i]),len(core[j]))>6 else .115
            need_y=.105
            if abs(d[0])<need_x and abs(d[1])<need_y:
                scaled=d/np.array([need_x,need_y])
                if np.linalg.norm(scaled)<1e-9:scaled=rng.normal(size=2)*.01
                push=scaled/(np.linalg.norm(scaled)+1e-9)*.003
                delta[i]+=push;delta[j]-=push
        a+=delta
        a=np.clip(a,.025,.975)
        if np.max(np.abs(delta))<1e-6:break
    pos={gene:tuple(a[i]) for i,gene in enumerate(core)}
    for i,gene in enumerate(strip):
        row,col=divmod(i,6)
        pos[gene]=((col+.5)/6,(2.5-row)/3)
    order=[gene for group in groups for gene in group]
    order += [x for x in NODES.index if x not in module and x in core]
    order += strip
    assert len(order)==len(set(order))==57
    records=[]
    for gene in NODES.index:
        records.append(dict(GENE=gene,community=module.get(gene,0),community_type='positive-edge module' if gene in module else 'singleton',
            display_region='network' if gene in core else 'shared_no_edge_strip',x=pos[gene][0],y=pos[gene][1],
            matrix_order=order.index(gene),node_area_pt2=AREA_FACTOR*abs(NODES.loc[gene,'normZ'])))
    pd.DataFrame(records).to_csv(DATA/'node_layout_and_modules.csv',index=False)
    (DATA/'communities.json').write_text(json.dumps(dict(algorithm='NetworkX Louvain weighted modularity',
        networkx_version=nx.__version__,input='Positive Pearson edges r >0.2 only',weight='Pearson r',resolution=1,
        seed=17,modularity=nx.community.modularity(g,communities,weight='weight',resolution=1),
        non_singleton_communities=groups,singletons=[x for x in NODES.index if x not in module]),indent=2),encoding='utf-8')
    return groups,module,core,strip,order,pos


def load_saved_structure():
    """Reuse the reviewed coordinates and modules without refitting the network."""
    layout=pd.read_csv(DATA/'node_layout_and_modules.csv').set_index('GENE')
    groups=json.loads((DATA/'communities.json').read_text(encoding='utf-8'))['non_singleton_communities']
    module={gene:i+1 for i,group in enumerate(groups) for gene in group}
    core=layout.index[layout.display_region.eq('network')].tolist()
    strip=layout.index[layout.display_region.eq('shared_no_edge_strip')].tolist()
    order=layout.sort_values('matrix_order').index.tolist()
    pos={gene:(float(row.x),float(row.y)) for gene,row in layout.iterrows()}
    assert len(core)==39 and len(strip)==18 and len(order)==57
    return groups,module,core,strip,order,pos


def hull(points):
    points=sorted(set(map(tuple,points)))
    def cross(o,a,b):return (a[0]-o[0])*(b[1]-o[1])-(a[1]-o[1])*(b[0]-o[0])
    lower=[]
    for p in points:
        while len(lower)>=2 and cross(lower[-2],lower[-1],p)<=0:lower.pop()
        lower.append(p)
    upper=[]
    for p in reversed(points):
        while len(upper)>=2 and cross(upper[-2],upper[-1],p)<=0:upper.pop()
        upper.append(p)
    a=np.array(lower[:-1]+upper[:-1]);c=a.mean(axis=0)
    return c+(a-c)*1.14


def node(ax,gene,xy,panel,strip=False,fontsize=None):
    r=NODES.loc[gene]; colour=DEPLETED if r.direction=='depleted' else ENRICHED
    face='#EDE1EF' if r.direction=='depleted' else 'white'
    if panel=='depmap' and not r.depmap_available:colour,face='#888888','#D6D6D6'
    ax.scatter(*xy,s=AREA_FACTOR*abs(r.normZ),facecolor=face,edgecolor=colour,lw=1.1,zorder=4)
    label=gene+'*' if panel=='depmap' and not r.depmap_available else gene
    radius=np.sqrt(AREA_FACTOR*abs(r.normZ))/2
    return ax.annotate(label,xy,xytext=(0,-radius-2.2 if strip else radius+1.8),textcoords='offset points',
        fontsize=fontsize if fontsize is not None else (7.5 if strip else 8),ha='center',va='top' if strip else 'bottom',color='#252525',
        bbox=dict(facecolor='white',edgecolor='none',alpha=LABEL_BOX_ALPHA,pad=.1),zorder=5)


def network_panel(fig,rect,side,groups,core,pos,fontsize=None):
    ax=fig.add_axes(rect);ax.set(xlim=(-.07,1.07),ylim=(-.08,1.12));ax.axis('off')
    for i,group in enumerate(groups):
        ax.add_patch(Polygon(hull([pos[x] for x in group]),facecolor=MODULE_COLS[i],edgecolor='none',zorder=0))
    edges=DE if side=='depmap' else ST
    for e in edges.sort_values('pearson_r' if side=='depmap' else 'string_score',key=abs).itertuples():
        x,y=zip(pos[e.gene_a],pos[e.gene_b])
        if side=='depmap':
            colour=POSITIVE if e.pearson_r>0 else NEGATIVE
            width=.35+2.3*abs(e.pearson_r);ls='-' if e.pearson_r>0 else (0,(3,2))
        else:colour='#A3A7AA';width=.25+2.1*e.string_score;ls='-'
        ax.plot(x,y,color=colour,lw=width,ls=ls,alpha=.62,zorder=1)
    return ax,[node(ax,g,pos[g],side,fontsize=fontsize) for g in core]


def strip_panel(fig,rect,side,strip,pos,fontsize=None):
    ax=fig.add_axes(rect);ax.set(xlim=(0,1),ylim=(-.12,1.12));ax.axis('off')
    ax.axhline(1.10,color='#DDDDDD',lw=.65)
    return ax,[node(ax,g,pos[g],side,strip=True,fontsize=fontsize) for g in strip]


def matrix_panel(fig,order,groups,rect=None,bar_rect=None,output_path=None):
    codes={'neither':0,'depmap_only':1,'string_only':2,'both':3,'depmap_unavailable':4,'string_with_depmap_unavailable':5}
    colors=['#F6F6F6','#80AFC7','#F6C8B6','#715095','#C7C7C7','#D8B7AA']
    cm=ListedColormap(colors);cm.set_bad('white')
    n=len(order)
    pairs=PAIRS.loc[PAIRS.gene_a.isin(order)&PAIRS.gene_b.isin(order)]
    counts=pairs.evidence_category.value_counts()
    a=np.full((n,n),np.nan);idx={g:i for i,g in enumerate(order)}
    neg=[]
    for r in pairs.itertuples():
        i,j=sorted([idx[r.gene_a],idx[r.gene_b]],reverse=True)
        a[i,j]=codes[r.evidence_category]
        if r.depmap_edge and r.pearson_r<0:neg.append((j,i))
    ax=fig.add_axes(rect if rect is not None else [.12,.075,.83,.445])
    ax.imshow(a,cmap=cm,vmin=-.5,vmax=5.5,interpolation='none',aspect='auto')
    # White cell borders keep sparse evidence readable, without a heavy full grid.
    for i in range(n+1):
        ax.plot([-.5,i-.5],[i-.5,i-.5],color='white',lw=.4,zorder=2)
        ax.plot([i-.5,i-.5],[i-.5,n-.5],color='white',lw=.4,zorder=2)
    for r in pairs.loc[~pairs.depmap_available].itertuples():
        i,j=sorted([idx[r.gene_a],idx[r.gene_b]],reverse=True)
        ax.add_patch(Rectangle((j-.5,i-.5),1,1,facecolor='none',edgecolor='#888888',
                              linewidth=0,hatch='///',zorder=2.5))
    for x,y in neg:
        ax.text(x,y,'−',fontsize=6.4,ha='center',va='center',
                color='white' if a[y,x]==codes['both'] else '#222222',zorder=4)
    labels=[x+'*' if x=='ZCCHC16' else x for x in order]
    ax.set(xticks=range(n),yticks=range(n),xticklabels=labels,yticklabels=labels,
           xlim=(-.5,n-.5),ylim=(n-.5,-.5))
    ax.tick_params(length=0,pad=2,labelsize=7.1)
    plt.setp(ax.get_xticklabels(),rotation=90,ha='center')
    for spine in ax.spines.values():spine.set_visible(False)
    # Module blocks outline comparisons within each positive-edge community.
    start=0
    for group in groups:
        end=start+len(group)
        ax.plot([start-.5,end-.5,end-.5],[end-.5,end-.5,start-.5],color='#B7B7B7',lw=.7,zorder=3)
        start=end
    effect_cm=LinearSegmentedColormap.from_list('drugz',[DEPLETED,'#FFFFFF',ENRICHED])
    lim=float(NODES.normZ.abs().max())
    bar=fig.add_axes(bar_rect if bar_rect is not None else [.12,.528,.83,.010])
    bar.imshow(NODES.loc[order,'normZ'].to_numpy()[None,:],cmap=effect_cm,vmin=-lim,vmax=lim,aspect='auto',interpolation='none')
    bar.set(xticks=[],yticks=[])
    for spine in bar.spines.values():spine.set_visible(False)
    # Use the unused triangle for keys, rather than shrinking the complete matrix.
    handles=[Patch(facecolor=colors[3],label=f'Both sources ({counts.get("both",0)} pairs)'),
             Patch(facecolor=colors[1],label=f'DepMap only ({counts.get("depmap_only",0)})'),
             Patch(facecolor=colors[2],label=f'STRING only ({counts.get("string_only",0)})'),
             Patch(facecolor=colors[0],edgecolor='#DDDDDD',label='Neither above threshold')]
    if (~pairs.depmap_available).any():
        handles.append(Patch(facecolor=colors[4],edgecolor='#888888',hatch='///',label='DepMap unavailable*'))
    ax.legend(handles=handles,loc='upper left',bbox_to_anchor=(.48,.955),frameon=False,
              fontsize=8.5,labelspacing=.8,handlelength=1.5,handleheight=1.1)
    ax.text(.49,.66,'−  Negative DepMap correlation',transform=ax.transAxes,fontsize=8.5)
    ax.text(.49,.605,f'{n} candidates · {len(pairs):,} unique gene pairs',transform=ax.transAxes,fontsize=8.3,color='#555555')
    ax.text(.49,.565,'Order: positive DepMap modules, then remaining hits',transform=ax.transAxes,fontsize=8,color='#555555')
    # Keep the entire colour key and its labels inside the unused upper triangle.
    cbax=ax.inset_axes([.65,.485,.31,.025])
    cb=fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(-lim,lim),cmap=effect_cm),cax=cbax,orientation='horizontal',ticks=[-10,-5,0,5,10])
    cb.outline.set_visible(False);cb.ax.tick_params(labelsize=7.5,length=2)
    cb.set_label('drugZ normZ (annotation strip)',fontsize=8,labelpad=3)
    matrix=pd.DataFrame(a,index=order,columns=order)
    matrix.to_csv(output_path if output_path is not None else DATA/'panel_c_evidence_codes.csv')
    return ax


def main():
    groups,module,core,strip,order,pos=build_structure()
    fig=plt.figure(figsize=(W,H))
    fig.text(.017,.985,'a',fontweight='semibold',fontsize=12,va='top')
    fig.text(.512,.985,'b',fontweight='semibold',fontsize=12,va='top')
    fig.text(.06,.982,'DepMap: |r| > 0.2',fontsize=9.5,va='top')
    fig.text(.555,.982,'STRING: combined score ≥ 0.4',fontsize=9.5,va='top')
    aa,la=network_panel(fig,[.025,.685,.455,.275],'depmap',groups,core,pos)
    ab,lb=network_panel(fig,[.520,.685,.455,.275],'string',groups,core,pos)
    assert len(aa.lines)==41 and len(ab.lines)==82
    assert len(aa.collections)==len(ab.collections)==39
    fig.text(.04,.677,'No displayed edge in either network (18)',fontsize=7.5,color='#666666')
    fig.text(.535,.677,'No displayed edge in either network (18)',fontsize=7.5,color='#666666')
    sa,lsa=strip_panel(fig,[.025,.591,.455,.075],'depmap',strip,pos)
    sb,lsb=strip_panel(fig,[.520,.591,.455,.075],'string',strip,pos)
    assert len(sa.collections)==len(sb.collections)==18
    for left,right,genes in [(aa,ab,core),(sa,sb,strip)]:
        for i,gene in enumerate(genes):
            assert np.array_equal(left.collections[i].get_offsets(),right.collections[i].get_offsets())
            assert np.array_equal(left.collections[i].get_sizes(),right.collections[i].get_sizes())
            assert np.isclose(left.collections[i].get_sizes()[0],AREA_FACTOR*abs(NODES.loc[gene,'normZ']))
    labels=la+lb+lsa+lsb
    # Shared keys distinguish screen direction, correlation sign and score magnitude.
    handles=[Line2D([],[],marker='o',mfc='#EDE1EF',mec=DEPLETED,lw=0,ms=7,label='LOY-depleted (48)'),
             Line2D([],[],marker='o',mfc='white',mec=ENRICHED,lw=0,ms=7,label='LOY-enriched (9)'),
             Line2D([],[],color=POSITIVE,lw=1.3,label='Positive r'),
             Line2D([],[],color=NEGATIVE,lw=1.3,ls='--',label='Negative r')]
    fig.legend(handles=handles,loc='center left',bbox_to_anchor=(.032,.567),ncol=4,fontsize=8,
               columnspacing=1.1,handletextpad=.5,handlelength=1.8,frameon=False)
    key=fig.add_axes([.71,.55,.27,.038]);key.axis('off');key.set(xlim=(0,1),ylim=(0,1))
    key.text(0,.8,'Node area: |drugZ normZ|',fontsize=8)
    for x,z in zip([.16,.47,.79],[4,8,12]):
        key.scatter(x,.28,s=AREA_FACTOR*z,facecolor='white',edgecolor='#777777',lw=.8)
        key.text(x+.08,.28,str(z),va='center',fontsize=7.8)
    fig.text(.04,.548,'Line width: |r| in a; STRING score in b. Shading: positive DepMap modules.',fontsize=7.6,color='#555555')
    fig.text(.017,.537,'c',fontweight='semibold',fontsize=12,va='top')
    matrix_panel(fig,order,groups)
    fig.text(.12,.016,'* ZCCHC16: no DepMap profile; retained in STRING and the complete candidate set.',fontsize=8,color='#555555')
    fig.canvas.draw();renderer=fig.canvas.get_renderer()
    overlaps=[]
    for i,a in enumerate(labels):
        for b in labels[i+1:]:
            if a.get_window_extent(renderer).overlaps(b.get_window_extent(renderer)):
                overlaps.append([a.get_text(),b.get_text()])
    review=OUT/'review';review.mkdir(exist_ok=True)
    for ext in ['pdf','svg','png']:fig.savefig(OUT/f'{STEM}.{ext}',dpi=600 if ext=='png' else 300)
    fig.savefig(review/f'{STEM}_preview.png',dpi=140)
    plt.close(fig)
    with fitz.open(OUT/f'{STEM}.pdf') as doc:
        page=doc[0];text=page.get_text()
        assert all(g in text for g in NODES.index)
        assert all('OpenSans' in f[3] for f in page.get_fonts())
        spans=[s for block in page.get_text('dict')['blocks'] if 'lines' in block for line in block['lines'] for s in line['spans'] if s['text'].strip()]
        outside=[s['text'] for s in spans if not page.rect.contains(fitz.Rect(s['bbox']))]
        page.get_pixmap(matrix=fitz.Matrix(1.25,1.25),colorspace=fitz.csGRAY).save(review/f'{STEM}_greyscale.png')
    record=dict(candidates=57,network_nodes_per_panel=57,core_nodes=39,shared_strip_nodes=18,
        shared_node_positions=True,area_factor_pt2_per_abs_normZ=AREA_FACTOR,
        node_area_range_pt2=[float(AREA_FACTOR*NODES.normZ.abs().min()),float(AREA_FACTOR*NODES.normZ.abs().max())],
        depmap_edges=41,string_edges=82,shared_edges=14,positive_communities=groups,
        matrix_genes=len(order),matrix_unique_pairs=1596,matrix_negative_glyphs=7,
        unavailable_pairs=56,missing_depmap=['ZCCHC16'],font='Open Sans',page_inches=[W,H],png_dpi=600,
        label_overlaps=overlaps,text_outside_page=outside,networkx_version=nx.__version__)
    (DATA/'figure_validation.json').write_text(json.dumps(record,indent=2),encoding='utf-8')
    manifest={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in
        [Path(__file__),DATA/'candidates.csv',DATA/'pairwise_evidence.csv',DATA/'node_layout_and_modules.csv']+
        [OUT/f'{STEM}.{ext}' for ext in ['pdf','svg','png']]}
    (DATA/'figure_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    print(json.dumps(record,indent=2))
    assert not overlaps and not outside, (overlaps,outside)


if __name__=='__main__':main()
