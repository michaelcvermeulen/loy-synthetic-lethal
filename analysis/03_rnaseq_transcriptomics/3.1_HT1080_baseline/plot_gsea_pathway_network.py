"""Publication and offline interactive views of the frozen enrichment maps."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import math
import textwrap
import html

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.collections import LineCollection
from adjustText import adjust_text
import networkx as nx
import numpy as np
import pandas as pd
from PIL import Image
from pypdf import PdfReader

ROOT = Path(__file__).resolve().parent
TABLES = ROOT / 'outputs/tables/gsea_network'
FIGURES = ROOT / 'outputs/figures/gsea_network'
PROV = ROOT / 'provenance/gsea_network'
COLORS = {'DC01':'#0072B2','DC02':'#56B4E9','DC03':'#009E73','DC04':'#6EA684',
          'DC05':'#CCAA45','DC06':'#E69F00','DC07':'#8C6BB1','DC08':'#D55E00','DC63':'#CC79A7',
          'UC01':'#E69F00','UC02':'#0072B2','UC03':'#009E73','UC04':'#8C6BB1',
          'UC05':'#CC79A7','UC06':'#D55E00'}
SHAPES = {'Both periods':'o','August only':'s','May only':'^','Combined only':'D'}
MARKER_AREA = {'o':math.pi/4,'s':1,'^':.5,'D':1}
DISPLAY_COMMUNITIES = {'UP':['UC05','UC02','UC01','UC06','UC03','UC04'],
                      'DOWN':['DC01','DC02','DC03','DC04','DC05','DC06','DC07','DC08','DC63']}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(name):
    return pd.read_csv(TABLES / (name+'.csv'), float_precision='round_trip')


def make_graph(nodes, edges, basis):
    graph = nx.Graph()
    graph.add_nodes_from(nodes.node_id)
    selected = edges[edges[f'{basis}_primary_edge'] & edges.source.isin(nodes.node_id)]
    graph.add_weighted_edges_from(selected[['source','target',f'{basis}_jaccard']].itertuples(index=False,name=None))
    return graph


def component_layout(graph):
    """Pack independently laid out connected components; do not imply distances."""
    groups = sorted((c for c in nx.connected_components(graph) if len(c)>1), key=lambda c:(-len(c),min(c)))
    circles = []; result = {}
    for group in groups:
        radius = math.sqrt(len(group))
        center = np.zeros(2)
        if circles:
            # Deterministic spiral search for an empty component footprint.
            for t in np.arange(.1, 250, .075):
                center = .6*t*np.array([math.cos(t), math.sin(t)])
                if all(np.linalg.norm(center-c) > radius+r+.5 for c,r in circles):
                    break
            else:
                raise RuntimeError('Component packing did not converge')
        circles.append((center, radius))
        ordered = graph.subgraph(sorted(group)).copy()
        local = nx.spring_layout(ordered, seed=42, iterations=300, scale=radius*.8, weight='weight')
        for n in ordered:
            result[n] = center+local[n]
    coords=np.array(list(result.values()))
    lo=coords.min(axis=0); hi=coords.max(axis=0); span=max(hi-lo)
    return {n:(p-(lo+hi)/2)/span*.84+.5 for n,p in result.items()}


def positions(nodes, edges):
    rows=[]
    for direction, subset in nodes.groupby('direction'):
        graph=make_graph(subset,edges,'leading_edge')
        layout=component_layout(graph)
        # Repel glyph footprints at final print dimensions so dense groups remain legible.
        ids=list(layout); points=np.array([layout[n] for n in ids])
        metadata=subset.set_index('node_id')
        footprint=metadata.loc[ids,'period_support'].map({'Both periods':1/math.sqrt(math.pi),
            'August only':math.sqrt(.5),'May only':1,'Combined only':math.sqrt(.5)}).to_numpy()
        radii=(np.sqrt(node_size(metadata.loc[ids,'neglog10_global_BH']))*footprint+.5)*25.4/72/(183*.445)
        minimum=radii[:,None]+radii[None,:]
        for _ in range(1500):
            delta=points[:,None,:]-points[None,:,:]
            distance=np.linalg.norm(delta,axis=2)
            np.fill_diagonal(distance,np.inf)
            overlap=np.maximum(minimum-distance,0)
            if overlap.max()<.0001:
                break
            force=delta/np.maximum(distance[:,:,None],1e-12)*overlap[:,:,None]*.35
            points+=force.sum(axis=1)
            points=np.clip(points,.035,.965)
        assert overlap.max()<.001, f'Unresolved node crowding: {direction}'
        for n,p in zip(ids,points):
            rows.append({'node_id':n,'direction':direction,'area':'connected','x':p[0],'y':p[1]})
        isolated=sorted(nx.isolates(graph))
        columns=35 if direction=='DOWN' else 10
        nrows=math.ceil(len(isolated)/columns)
        for i,n in enumerate(isolated):
            rows.append({'node_id':n,'direction':direction,'area':'isolated',
                         'x':.035+.93*(i%columns+.5)/columns,
                         'y':.94-.87*(i//columns+.5)/max(nrows,1)})
    return pd.DataFrame(rows)


def node_size(logq):
    return 4+4*np.asarray(logq)


def draw_nodes(ax, nodes, xy, size_factor=1):
    for support, marker in SHAPES.items():
        group=nodes[nodes.period_support.eq(support)]
        if group.empty:
            continue
        ax.scatter([xy[n][0] for n in group.node_id],[xy[n][1] for n in group.node_id],
                   s=node_size(group.neglog10_global_BH)*size_factor/MARKER_AREA[marker],marker=marker,
                   c=group.color,edgecolors='#48535B',linewidths=.28,zorder=3,alpha=.92)


def write_figure(nodes, edges, coords, communities, labels, basis):
    fig=plt.figure(figsize=(183/25.4,235/25.4),facecolor='white')
    primary=basis=='leading_edge'
    title='GSEA pathways connected by shared leading-edge genes' if primary else 'GSEA pathways connected by full gene-set overlap'
    fig.text(.035,.978,title,fontsize=8.5,fontweight='bold',va='top')
    sub='All 808 combined global-BH < 0.05 sets; H and C1–C9, MSigDB 2026.1.Hs'
    fig.text(.035,.955,sub,fontsize=6.2,color='#53616A',va='top')
    if not primary:
        fig.text(.035,.937,'Same nodes, positions and leading-edge community colours; only connections change.',fontsize=5.8,color='#53616A')
    for index,direction in enumerate(['UP','DOWN']):
        x0=.025+index*.505; width=.445
        subset=nodes[nodes.direction.eq(direction)]
        graph=make_graph(subset,edges,basis)
        placed=coords[coords.direction.eq(direction)].set_index('node_id')
        globalxy={n:(x0+width*r.x, (.525+.365*r.y) if r.area=='connected' else (.328+.15*r.y)) for n,r in placed.iterrows()}
        # One axes spans both connected and isolate areas so sensitivity edges can link either.
        ax=fig.add_axes([x0,.328,width,.562]); ax.set_xlim(0,1); ax.set_ylim(0,1); ax.axis('off')
        xy={n:np.array([(p[0]-x0)/width,(p[1]-.328)/.562]) for n,p in globalxy.items()}
        lines=[[xy[u],xy[v]] for u,v in graph.edges()]
        ax.add_collection(LineCollection(lines,colors='#A8B0B6',linewidths=[.15+.65*graph[u][v]['weight'] for u,v in graph.edges()],alpha=.42,zorder=1))
        draw_nodes(ax,subset,xy)
        fig.text(x0,.919,'ab'[index],fontsize=8,fontweight='bold')
        heading='Higher RNA in LOY' if direction=='UP' else 'Lower RNA in LOY'
        fig.text(x0+.025,.919,heading,fontsize=7.3,fontweight='bold')
        fig.text(x0+.025,.902,f'{len(subset)} gene sets · {graph.number_of_edges():,} connections',fontsize=6,color='#53616A')
        isolates=int(subset.leading_edge_degree.eq(0).sum())
        note=f'{isolates} sets without a leading-edge connection'
        if not primary:
            note=f'{isolates} originally isolated sets (positions retained)'
        fig.text(x0+.01,.496,note,fontsize=5.7,color='#53616A')
        fig.add_artist(Line2D([x0+.01,x0+width-.01],[.486,.486],transform=fig.transFigure,color='#CBD1D5',linewidth=.4))
        annotations=[]
        for cid in DISPLAY_COMMUNITIES[direction]:
            c=communities.loc[cid]; members=subset[subset.community_id.eq(cid)]
            target=np.mean([xy[n] for n in members.node_id],axis=0)
            annotations.append(ax.text(*target,cid,fontsize=5.5,fontweight='bold',color='#202B32',
                ha='center',va='center',bbox={'facecolor':'white','edgecolor':'none','alpha':.88,'pad':1},zorder=5))
        np.random.seed(42)
        adjust_text(annotations,x=[p[0] for p in xy.values()],y=[p[1] for p in xy.values()],ax=ax,
                    expand=(1.2,1.25),force_static=(.2,.3),force_text=(.2,.3),
                    arrowprops={'arrowstyle':'-','color':'#53616A','lw':.4},iter_lim=300)
        for i,cid in enumerate(DISPLAY_COMMUNITIES[direction]):
            y=.302-i*.0165
            fig.add_artist(Line2D([x0+.014],[y],marker='o',markersize=3,linestyle='none',
                           markerfacecolor=COLORS[cid],markeredgecolor='#53616A',markeredgewidth=.3,transform=fig.transFigure))
            n=int(communities.loc[cid,'pathways'])
            fig.text(x0+.026,y,f'{cid}  {labels[cid]} ({n})',fontsize=5.5,va='center')
        if direction=='UP':
            fig.text(x0+.01,.183,'Each node is one enriched gene set.\nLabels identify selected overlap communities;\ncounts in parentheses give their numbers of sets.',fontsize=5.7,
                     linespacing=1.5,color='#53616A',va='top')
        else:
            fig.text(x0+.01,.143,'Grey: other communities and isolated sets.',fontsize=5.5,color='#53616A')
    shape_handles=[Line2D([0],[0],marker=v,markersize=3.5,linestyle='none',markerfacecolor='#B8BEC3',
                   markeredgecolor='#48535B',markeredgewidth=.35,label=k) for k,v in SHAPES.items()]
    fig.text(.035,.121,'Period support: global BH < 0.05 with the combined NES direction',fontsize=6,fontweight='bold')
    fig.legend(handles=shape_handles,loc='upper left',bbox_to_anchor=(.025,.115),ncol=4,frameon=False,
               fontsize=5.8,columnspacing=1.5,handletextpad=.4)
    size_handles=[plt.scatter([],[],s=node_size(v)/MARKER_AREA['o'],facecolor='#B8BEC3',edgecolor='#48535B',linewidth=.35) for v in [2,4,6]]
    fig.legend(size_handles,['2','4','6'],title='Node area: −log10(global BH)',title_fontsize=5.8,fontsize=5.8,
               loc='upper left',bbox_to_anchor=(.025,.083),ncol=3,frameon=False,columnspacing=1,handletextpad=.5)
    fig.text(.54,.073,'Edge width: Jaccard similarity ≥ 0.25\nand at least 3 shared genes',fontsize=5.8,linespacing=1.5)
    fig.text(.035,.036,'Communities and spacing summarize gene overlap; positions are arbitrary. Colours identify communities within each panel.',fontsize=5.3,color='#53616A')
    fig.text(.035,.022,'Exploratory enrichment in repeatedly sampled populations; estimated-correlation CAMERA has no combined global-BH discoveries.',fontsize=5.3,color='#53616A')
    fig.canvas.draw()
    stem='HT1080_GSEA_leading_edge_network' if primary else 'HT1080_GSEA_full_membership_network'
    for ext in ['pdf','svg','png']:
        fig.savefig(FIGURES/f'{stem}.{ext}',dpi=600,facecolor='white')
    plt.close(fig)
    return stem


def write_html(nodes, edges, coords, communities, labels):
    data=nodes.copy().fillna('')
    points=coords.set_index('node_id')
    for axis in ['x','y']:
        data[axis]=data.node_id.map(points[axis])
    data['area']=data.node_id.map(points.area)
    data['community_label']=data.community_id.map(lambda c:labels.get(c,communities.loc[c,'representative_pathway'].replace('_',' ')))
    # Sensitivity slider uses all eligible edges down to Jaccard 0.20.
    pairs=pd.read_csv(TABLES/'all_nonzero_gene_overlap_pairs.csv.gz',float_precision='round_trip')
    pairs=pairs[((pairs.leading_edge_jaccard>=.2)&(pairs.leading_edge_shared_n>=3))|
                ((pairs.full_set_jaccard>=.2)&(pairs.full_set_shared_n>=3))]
    keep=['source','target','direction','leading_edge_jaccard','full_set_jaccard',
          'leading_edge_shared_n','full_set_shared_n','shared_leading_edge_genes']
    payload={'nodes':data.to_dict('records'),'edges':pairs[keep].fillna('').to_dict('records'),
             'labels':labels,'display':DISPLAY_COMMUNITIES,'colors':COLORS}
    template=(ROOT/'report/gsea_network_template.html').read_text(encoding='utf-8')
    output=template.replace('__NETWORK_DATA__',json.dumps(payload,ensure_ascii=False).replace('</','<\\/'))
    report=ROOT/'outputs/report/gsea_pathway_network.html'
    report.write_text(output,encoding='utf-8')
    return report


def main():
    FIGURES.mkdir(parents=True,exist_ok=True)
    nodes=read('network_nodes'); edges=read('network_edges'); communities=read('communities').set_index('community_id')
    labels=json.loads((ROOT/'config/gsea_network_labels.json').read_text(encoding='utf-8'))
    assert all(cid in communities.index for cid in labels)
    before={p.relative_to(ROOT).as_posix():sha(p) for p in [TABLES/'network_nodes.csv',TABLES/'network_edges.csv',TABLES/'communities.csv']}
    nodes['color']=nodes.community_id.map(COLORS).fillna('#C7CBCE')
    coords=positions(nodes,edges)
    assert coords.node_id.is_unique and set(coords.node_id)==set(nodes.node_id)
    coords.to_csv(TABLES/'plot_coordinates.csv',index=False,float_format='%.17g')
    plt.rcParams.update({'font.family':'Arial','font.size':6,'pdf.fonttype':42,'ps.fonttype':42,
                         'svg.fonttype':'none','svg.hashsalt':'HT1080_GSEA_network_seed42'})
    stems=[write_figure(nodes,edges,coords,communities,labels,b) for b in ['leading_edge','full_set']]
    report=write_html(nodes,edges,coords,communities,labels)
    caption=('**Gene-set similarity networks for the HT1080 LOY transcriptional ranking.** '
        'a, Higher-RNA enrichment (32 sets); b, lower-RNA enrichment (776 sets). Every combined fgsea result with original '
        'global Benjamini–Hochberg q < 0.05 is retained from the 23,696 eligible MSigDB 2026.1.Hs H and C1–C9 tests. '
        'Nodes represent gene sets, and node area is 4 + 4 × [−log10(q)] square points. Connections require at least three '
        'shared combined-ranking leading-edge genes and Jaccard similarity (intersection divided by union) ≥ 0.25; '
        'line width increases with similarity. Only within-direction relationships are displayed. Weighted Louvain communities '
        '(resolution 1, seed 42) define the groups; selected groups have descriptive labels and colours, while other groups are grey. '
        'Community colours are local to each direction. Source cell-type or phenotype labels describe the originating signatures, '
        'not inferred HT1080 cell identity. Shapes show global q < 0.05 and agreement with the combined NES direction in August, '
        'May, both periods or neither period (combined only). These categories describe threshold support, not a statistical test '
        'of differences between periods. All 808 sets have concordant NES signs in both periods. The 19 higher-RNA and 489 lower-RNA '
        'sets without qualifying leading-edge edges occupy separate strips. Component packing, spring layouts and isolate ordering '
        'are display choices; distances have no biological meaning. The companion full-membership figure uses the same nodes, colours '
        'and coordinates, but computes edges from complete set membership intersected with the 25,461-gene ranking universe. '
        'The offline interactive report displays exact names, statistics, memberships and shared genes, and offers Jaccard thresholds '
        '0.20, 0.25, 0.35 and 0.50. No enrichment tests or BH adjustments were rerun. These are exploratory patterns in repeated samples '
        'of the same established populations; estimated-correlation CAMERA has no combined global-BH discoveries. '
        'Neither gene-set similarity nor source chromatin signatures demonstrate protein interaction or measured HT1080 chromatin activity.')
    caption_file=ROOT/'outputs/captions/gsea_pathway_network.md'
    caption_file.write_text('# GSEA pathway networks\n\n'+caption+'\n\n'
        'Method background: [EnrichmentMap similarity definitions](https://enrichmentmap.readthedocs.io/en/latest/Parameters.html) '
        'and [GSEA leading-edge definition](https://enrichmentmap.readthedocs.io/en/docs-2.2/AdditionalFeatures.html). '
        'This extension uses a local NetworkX implementation of an enrichment map.\n\n'
        '[Leading-edge PDF](../figures/gsea_network/HT1080_GSEA_leading_edge_network.pdf) · '
        '[Full-membership PDF](../figures/gsea_network/HT1080_GSEA_full_membership_network.pdf) · '
        '[Interactive report](../report/gsea_pathway_network.html) · '
        '[Node statistics](../tables/gsea_network/network_nodes.csv) · '
        '[Edge similarities and genes](../tables/gsea_network/network_edges.csv)\n',encoding='utf-8')
    checks=[]
    for stem in stems:
        pdf=PdfReader(FIGURES/f'{stem}.pdf'); assert len(pdf.pages)==1
        page=pdf.pages[0]; text=page.extract_text()
        assert abs(float(page.mediabox.width)*25.4/72-183)<.01
        assert abs(float(page.mediabox.height)*25.4/72-235)<.01
        for cid in COLORS:
            assert cid in text, cid
        fonts=[]
        for ref in page['/Resources']['/Font'].values():
            font=ref.get_object()
            for item in font.get('/DescendantFonts',[font]):
                f=item.get_object(); desc=f['/FontDescriptor'].get_object()
                assert any(k in desc for k in ['/FontFile','/FontFile2','/FontFile3'])
                assert 'Arial' in str(f.get('/BaseFont')); fonts.append(str(f.get('/BaseFont')))
        with Image.open(FIGURES/f'{stem}.png') as im:
            assert abs(im.width-183/25.4*600)<1.1 and abs(im.height-235/25.4*600)<1.1
            assert all(abs(d-600)<1 for d in im.info['dpi'])
        checks.append({'figure':stem,'dimensions_mm':[183,235],'png_dpi':600,'embedded_fonts':fonts})
    assert before=={p:sha(ROOT/p) for p in before}
    outputs=[FIGURES/f'{s}.{e}' for s in stems for e in ['pdf','svg','png']]+[report,caption_file,TABLES/'plot_coordinates.csv']
    result={'created_utc':datetime.now(timezone.utc).isoformat(),'source_hashes':before,
            'script_sha256':sha(__file__),'nodes_rendered':len(coords),'figure_checks':checks,
            'outputs_sha256':{p.relative_to(ROOT).as_posix():sha(p) for p in outputs},
            'visual_review':'pending','interactive_review':'pending'}
    (PROV/'figure_completion.json').write_text(json.dumps(result,indent=2))
    print(json.dumps({'figures':stems,'nodes_rendered':len(coords),'report':str(report)},indent=2))


if __name__=='__main__':
    main()
