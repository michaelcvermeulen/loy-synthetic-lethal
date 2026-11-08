"""DU145 counterparts of the selected HT1080 Metascape figures.

Read saved full-input enrichment statistics and the exported term network.
No enrichment tests are rerun. Preserve evidence-filtered export attributes
separately, rather than combining them with full-list hit counts.
"""
from collections import Counter
from pathlib import Path
import hashlib
import json
import sys
import textwrap
import xml.etree.ElementTree as ET

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import MaxNLocator
import networkx as nx
import numpy as np
import pandas as pd

FIG = Path(__file__).resolve().parents[1]
THESIS = FIG.parents[1]
SOURCE = FIG / 'source_data/du145_metascape_current'
SUPP = FIG / 'source_data/du145_metascape_terms_and_gene_modules'
SUPP.mkdir(exist_ok=True)
sys.path.insert(0, str(THESIS / 'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
apply_style(8)
BLUE, RED, DARK = '#3C5488', '#E64B35', '#444444'
PALETTE = [BLUE, RED, '#00A087', '#F39B7F', '#4DBBD5', '#8491B4', '#7E6148']
MM = 1/25.4
NS = {'x': 'http://www.cs.rpi.edu/XGMML'}
THEMES = {
    'HALLMARK EPITHELIAL MESENCHYMAL TRANSITION': 'Matrix and mesenchymal signatures',
    'response to nutrient': 'Nutrient response and migration',
    'HELLEBREKERS SILENCED DURING TUMOR ANGIOGENESIS': 'Angiogenesis and secretion signatures',
    'LENAOUR DENDRITIC CELL MATURATION UP': 'Dendritic-cell and perturbation signatures',
    'intracellular lipid transport': 'Lipid transport and hormone signatures',
    'heparin binding': 'Heparin binding and cell projections',
    'COLDREN GEFITINIB RESISTANCE DN': 'Epithelial and CDH1-associated signatures',
    'SATO SILENCED BY METHYLATION IN PANCREATIC CANCER 1': 'Methylation and Polycomb signatures',
    'positive regulation of interleukin-12 production': 'Cytokine regulation and chemotaxis',
    'NIKOLSKY BREAST CANCER 20P13 AMPLICON': '20p13 positional signature',
    'ROZANOV MMP14 TARGETS UP': 'MMP14 and SATB1 target signatures',
    'Amoebiasis': 'Innate-immune signalling signatures',
}
SHORT = {
    'HALLMARK EPITHELIAL MESENCHYMAL TRANSITION': 'Hallmark epithelial–mesenchymal transition',
    'TURASHVILI BREAST LOBULAR CARCINOMA VS LOBULAR NORMAL DN': 'Turashvili: lobular carcinoma vs normal (down)',
    'HELLEBREKERS SILENCED DURING TUMOR ANGIOGENESIS': 'Hellebrekers: tumour-angiogenesis silencing',
    'YAO TEMPORAL RESPONSE TO PROGESTERONE CLUSTER 16': 'Yao: progesterone response, cluster 16',
    'extracellular matrix structural constituent conferring tensile strength': 'ECM structural constituent: tensile strength',
    'external encapsulating structure organization': 'External encapsulating structure organisation',
    'phenol-containing compound metabolic process': 'Phenol-containing compound metabolism',
    'NABA CORE MATRISOME': 'NABA core matrisome',
    'NABA COLLAGENS': 'NABA collagens',
    'WU CELL MIGRATION': 'Wu: cell migration',
    'LENAOUR DENDRITIC CELL MATURATION UP': 'Lenaour: dendritic-cell maturation (up)',
    'LINDSTEDT DENDRITIC CELL MATURATION D': 'Lindstedt: dendritic-cell maturation (down)',
    'COLDREN GEFITINIB RESISTANCE DN': 'Coldren: gefitinib resistance (down)',
    'ONDER CDH1 TARGETS 2 DN': 'Onder: CDH1 targets 2 (down)',
    'MCBRYAN PUBERTAL BREAST 4 5WK UP': 'McBryan: pubertal breast, 4–5 weeks (up)',
    'CHARAFE BREAST CANCER LUMINAL VS MESENCHYMAL UP': 'Charafe: luminal vs mesenchymal (up)',
    'CHARAFE BREAST CANCER BASAL VS MESENCHYMAL UP': 'Charafe: basal vs mesenchymal (up)',
    'LEE BMP2 TARGETS UP': 'Lee: BMP2 targets (up)',
    'HOLLERN EMT BREAST TUMOR DN': 'Hollern: breast-tumour EMT (down)',
    'AIGNER ZEB1 TARGETS': 'Aigner: ZEB1 targets',
    'SATO SILENCED BY METHYLATION IN PANCREATIC CANCER 1': 'Sato: silenced by methylation in pancreatic cancer',
    'DELYS THYROID CANCER UP': 'Delys: thyroid cancer (up)',
    'GOZGIT ESR1 TARGETS DN': 'Gozgit: ESR1 targets (down)',
    'GOZGIT ESR1 TARGETS UP': 'Gozgit: ESR1 targets (up)',
    'MISSIAGLIA REGULATED BY METHYLATION UP': 'Missiaglia: methylation-regulated (up)',
    'RASHI NFKB1 TARGETS': 'Rashi: NFKB1 targets',
    'NIKOLSKY BREAST CANCER 20P13 AMPLICON': 'Nikolsky: breast-cancer 20p13 amplicon',
    'RASHI RESPONSE TO IONIZING RADIATION 2': 'Rashi: ionising-radiation response 2',
    'ROZANOV MMP14 TARGETS UP': 'Rozanov: MMP14 targets (up)',
    'HAN SATB1 TARGETS UP': 'Han: SATB1 targets (up)',
    'HAN SATB1 TARGETS DN': 'Han: SATB1 targets (down)',
    'positive regulation of interleukin-12 production': 'Positive regulation of IL-12 production',
    'TARTE PLASMA CELL VS PLASMABLAST UP': 'Tarte: plasma cell vs plasmablast (up)',
    'MANALO HYPOXIA UP': 'Manalo: hypoxia (up)',
    'LEE LIVER CANCER CIPROFIBRATE UP': 'Lee: liver-cancer ciprofibrate response (up)',
    'MCBRYAN PUBERTAL BREAST 3 4WK UP': 'McBryan: pubertal breast, 3–4 weeks (up)',
    'ONDER CDH1 TARGETS 3 DN': 'Onder: CDH1 targets 3 (down)',
}


def wrap(s, width=34):
    return '\n'.join(textwrap.wrap(str(s), width=width, break_long_words=False))


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def load(direction):
    full = pd.read_csv(SOURCE / direction / 'Enrichment_GO/GO_AllLists.csv')
    assert full.GO.is_unique
    expected = 31 if direction == 'up' else 83
    assert full['#GeneInHitList'].eq(expected).all()
    root = ET.parse(SOURCE / direction / 'Enrichment_GO/GONetwork.xgmml').getroot()
    original, lookup = [], {}
    for n in root.findall('x:node', NS):
        a = {v.get('name'): v.get('value') for v in n.findall('x:att', NS)}
        original.append(a)
        lookup[n.get('id')] = a['canonicalName']
    pd.DataFrame(original).to_csv(SOURCE / f'{direction}_original_network_attributes.csv', index=False)
    graph = nx.Graph()
    graph.add_nodes_from(lookup.values())
    edges = [(lookup[e.get('source')], lookup[e.get('target')]) for e in root.findall('x:edge', NS)]
    graph.add_edges_from((a, b) for a, b in edges if a != b)
    assert set(graph) <= set(full.GO)
    nodes = full.set_index('GO').loc[list(graph)].reset_index()
    nodes = nodes.rename(columns={'GO':'name', 'Description':'description', 'Hits':'hits'})
    nodes['mlog10p'] = -nodes.LogP
    nodes['log10q'] = nodes['Log(q-value)']
    nodes['n_hits'] = nodes['#GeneInGOAndHitList']
    nodes['display_label'] = nodes.description.map(lambda s: SHORT.get(s,s))
    communities = nx.community.louvain_communities(graph, seed=42)
    # Rank communities by their strongest full-input term, using stable tie-breaks.
    communities = sorted(communities, key=lambda c: (nodes.loc[nodes.name.isin(c), 'LogP'].min(), sorted(c)))
    summaries = []
    for cid, members in enumerate(communities, 1):
        subset = nodes[nodes.name.isin(members)].sort_values(['LogP', 'description'])
        top = subset.description.iloc[0]
        assert top in THEMES, top
        color = '#9B9B9B' if len(members) == 1 else PALETTE[cid-1]
        nodes.loc[nodes.name.isin(members),'community'] = cid
        nodes.loc[nodes.name.isin(members),'color'] = color
        nodes.loc[nodes.name.isin(members),'theme'] = THEMES[top]
        counts = Counter(g for hits in subset.hits for g in set(hits.split('|')))
        genes = [g for g,n in sorted(counts.items(), key=lambda pair: (-pair[1],pair[0]))[:4]]
        summaries.append({'community':cid,'theme':THEMES[top],'color':color,'terms':len(members),
                          'best_mlog10p':-subset.LogP.iloc[0], 'representative_1':subset.display_label.iloc[0],
                          'representative_2':subset.display_label.iloc[1] if len(subset)>1 else '',
                          'recurring_genes':', '.join(genes)})
    nodes['community'] = nodes.community.astype(int)
    assert (nodes.hits.str.split('|').map(len) == nodes.n_hits).all()
    return full, nodes, graph, pd.DataFrame(summaries)


def save(fig, stem):
    exports = []
    for ext in ['pdf','svg','png']:
        target = FIG/f'{stem}.{ext}'
        fig.savefig(target,dpi=600)
        exports.append({'file':target.name,'sha256':sha(target)})
    fig.savefig(FIG/'previews'/f'{stem}_preview.png',dpi=160)
    return exports


def network(direction, nodes, graph, summary):
    # Compact 3 x 2 display anchors; distances carry no biological meaning.
    centers = [(-1.08,.65),(0,.65),(1.08,.65),(-1.08,-.62),(0,-.62),(1.08,-.62)]
    positions = {}
    for row, center in zip(summary.itertuples(),centers):
        ids = nodes.loc[nodes.community.eq(row.community),'name']
        sub = graph.subgraph(ids)
        local = nx.spring_layout(sub,seed=42,k=.65,iterations=700)
        radius = .35 if len(ids)>1 else 0
        for g,xy in local.items():positions[g] = np.asarray(center)+radius*xy
    nodes = nodes.copy()
    nodes['x'] = nodes.name.map(lambda g:positions[g][0])
    nodes['y'] = nodes.name.map(lambda g:positions[g][1])
    fig = plt.figure(figsize=(183*MM,180*MM))
    condition = 'Higher' if direction == 'up' else 'Lower'
    input_count = 31 if direction == 'up' else 83
    fig.text(.035,.977,f'{condition} expression in Y-depleted DU145  |  {input_count} genes  |  {len(nodes)} terms',
             fontsize=7.3,color=DARK,va='top')
    ax = fig.add_axes([.035,.36,.93,.59])
    ax.set_axis_off()
    nx.draw_networkx_edges(graph,positions,ax=ax,width=.43,edge_color='#AAAAAA',alpha=.27)
    nx.draw_networkx_nodes(graph,positions,ax=ax,nodelist=nodes.name.tolist(),
                          node_size=12+4*nodes.mlog10p.to_numpy(),node_color=nodes.color.tolist(),
                          edgecolors='#444444',linewidths=.35)
    for row, (x,y) in zip(summary.itertuples(),centers):
        ax.text(x,y+.51,wrap(row.theme,24),fontsize=7.2,ha='center',va='bottom',weight='semibold')
    ax.set_xlim(-1.62,1.62);ax.set_ylim(-1.12,1.43)
    values = [3,6,9] if direction=='up' else [5,10,15]
    handles = [ax.scatter([],[],s=12+4*v,c='#BDBDBD',edgecolors='#444444',linewidths=.35) for v in values]
    fig.legend(handles,[str(v) for v in values],title=r'$-\log_{10}P$',loc='center',
               bbox_to_anchor=(.5,.329),ncol=3,frameon=False,fontsize=7,title_fontsize=7)
    table = fig.add_axes([.025,.024,.95,.267]);table.set_axis_off()
    table.set_xlim(0,1);table.set_ylim(0,1)
    columns = [.025,.342,.414,.50,.82]
    for x, header in zip(columns,['Theme','Terms','Best\n−log₁₀P','Representative terms','Recurring genes']):
        table.text(x,.99,header,fontsize=6.5,weight='semibold',va='top')
    table.plot([0,1],[.845,.845],color='#555555',lw=.55)
    for i,row in enumerate(summary.itertuples()):
        y=.81-i*.132
        table.scatter(.01,y-.028,s=18,color=row.color)
        terms='\n'.join(textwrap.shorten(s,width=35,placeholder='…') for s in [row.representative_1,row.representative_2] if s)
        vals=[wrap(row.theme,29),str(row.terms),f'{row.best_mlog10p:.1f}',terms,wrap(row.recurring_genes,20)]
        for x,v in zip(columns,vals):table.text(x,y,v,fontsize=6.0,va='top',linespacing=1.15)
    nodes.to_csv(SOURCE/f'{direction}_nodes.csv',index=False)
    summary.to_csv(SOURCE/f'{direction}_communities.csv',index=False)
    pd.DataFrame(graph.edges(),columns=['source','target']).to_csv(SOURCE/f'{direction}_edges.csv',index=False)
    stem=f'ch03_supp_du145_metascape_{direction}_network'
    return fig,save(fig,stem)


def gene_counts(nodes, summary):
    records=[]
    genes=sorted({g for hits in nodes.hits for g in hits.split('|')})
    for gene in genes:
        member=nodes[nodes.hits.str.split('|').map(lambda gs:gene in gs)]
        counts=member.community.value_counts()
        maximum=int(counts.max())
        dominant=int(min(counts[counts.eq(maximum)].index))
        row={'gene':gene,'total_term_hits':len(member),'dominant_module':dominant,
             'dominant_tie':int(counts.eq(maximum).sum())>1}
        for cid in summary.community:row[f'module_{cid}']=int(counts.get(cid,0))
        records.append(row)
    d=pd.DataFrame(records).sort_values(['total_term_hits','gene'],ascending=[False,True])
    d.to_csv(SUPP/'all_up_gene_module_counts.csv',index=False)
    top=d.head(40).sort_values(['dominant_module','total_term_hits','gene'],ascending=[True,False,True]).reset_index(drop=True)
    top.to_csv(SUPP/'displayed_up_gene_module_counts.csv',index=False)
    summary.to_csv(SUPP/'up_module_key.csv',index=False)
    assert sum(d.total_term_hits)==sum(nodes.n_hits)
    return top


def top_terms(full,direction):
    d=full.sort_values(['LogP','Log(q-value)','Description']).head(20).copy()
    d['display_label']=d.Description.map(lambda s: SHORT.get(s,s))
    d.to_csv(SUPP/f'{direction}_top20.csv',index=False)
    return d


def draw_terms(ax,direction,full):
    d=top_terms(full,direction).iloc[::-1]
    yy=np.arange(len(d));score=-d.LogP.to_numpy()
    ax.barh(yy,score,height=.67,color=RED if direction=='up' else BLUE,alpha=.9)
    ax.set_yticks(yy,[wrap(s,49) for s in d.display_label],fontsize=6.2)
    ax.tick_params(axis='y',length=0,pad=4)
    ax.set_xlabel(r'Enrichment $-\log_{10}P$',fontsize=7)
    ax.spines[['top','right','left']].set_visible(False)
    ax.set_xlim(0,score.max()*1.72)
    ax.set_ylim(-.6,len(d)-.3)
    ax.grid(axis='x',color='#E7E7E7',lw=.5);ax.set_axisbelow(True)
    for y,(_,r) in zip(yy,d.iterrows()):
        q=10.**r['Log(q-value)']
        ax.text(-r.LogP+.13,y,f"{int(r['#GeneInGOAndHitList'])}; q≈{q:.2g}",fontsize=5.9,va='center',color=DARK)
    label='Higher expression (31 genes)' if direction=='up' else 'Lower expression (83 genes)'
    ax.text(0,1.022,label,transform=ax.transAxes,fontsize=7.8,weight='semibold',color=RED if direction=='up' else BLUE)


def composite(full_up,full_down,nodes_up,summary_up):
    genes=gene_counts(nodes_up,summary_up)
    fig=plt.figure(figsize=(290*MM,210*MM))
    au=fig.add_axes([.274,.545,.265,.412])
    ad=fig.add_axes([.274,.074,.265,.412])
    ag=fig.add_axes([.635,.182,.345,.775])
    draw_terms(au,'up',full_up);draw_terms(ad,'down',full_down)
    au.set_xlabel('')
    y=[];cursor=.95;bounds=[]
    for cid,g in genes.groupby('dominant_module',sort=True):
        start=cursor
        for _ in range(len(g)):y.append(cursor);cursor+=1
        bounds.append((cid,start,cursor-1));cursor+=1.5
    left=np.zeros(len(genes))
    for row in summary_up.itertuples():
        v=genes[f'module_{row.community}'].to_numpy()
        ag.barh(y,v,left=left,height=.71,color=row.color,edgecolor='white',lw=.25,label=wrap(row.theme,30))
        left+=v
    ag.set_yticks(y,genes.gene,fontsize=7)
    ag.tick_params(axis='y',length=0,pad=3)
    ag.invert_yaxis();ag.set_ylim(max(y)+.7,-.55)
    ag.set_xlim(0,max(left)+4)
    ag.xaxis.set_major_locator(MaxNLocator(integer=True))
    ag.set_xlabel('Number of higher-expression network terms\ncontaining the gene',fontsize=7.3)
    ag.spines[['top','right','left']].set_visible(False)
    ag.grid(axis='x',color='#E7E7E7',lw=.5);ag.set_axisbelow(True)
    for row,yy in zip(genes.itertuples(),y):
        ag.text(row.total_term_hits+.4,yy,str(row.total_term_hits)+('†' if row.dominant_tie else ''),va='center',fontsize=6.3)
    community_key=summary_up.set_index('community')
    for cid,start,end in bounds:
        row=community_key.loc[cid]
        ag.scatter(.3,start-.85,s=15,marker='s',color=row.color,clip_on=False,zorder=4)
        ag.text(.9,start-.85,row.theme,fontsize=6.8,weight='semibold',va='center',color=DARK)
        if end<max(y):ag.axhline(end+.72,color='#CCCCCC',lw=.6)
    ag.text(0,1.022,'Recurring genes, grouped by pathway community',transform=ag.transAxes,fontsize=7.8,weight='semibold')
    for x,ypos,letter in [(.018,.963,'a'),(.018,.492,'b'),(.573,.963,'c')]:fig.text(x,ypos,letter,fontsize=12,weight='semibold')
    handles,labels=ag.get_legend_handles_labels()
    fig.legend(handles,labels,loc='lower left',bbox_to_anchor=(.603,.025),ncol=2,
               frameon=False,fontsize=6.1,handlelength=1.1,columnspacing=1.3,labelspacing=.65)
    return fig,save(fig,'ch03_supp_du145_metascape_terms_and_gene_modules')


def main():
    inputs={};figs=[];exports=[];records={}
    for direction in ['up','down']:
        full,nodes,g,summary=load(direction)
        inputs[direction]=(full,nodes,g,summary)
        fig,files=network(direction,nodes,g,summary)
        figs.append(fig);exports+=files
        records[direction]={'input_genes':31 if direction=='up' else 83,'full_input_terms':len(full),
                            'network_nodes':len(nodes),'network_edges':g.number_of_edges(),
                            'isolates_displayed':len(list(nx.isolates(g))),'communities':len(summary)}
    fup,nup,_,sup=inputs['up'];fdown=inputs['down'][0]
    fig,files=composite(fup,fdown,nup,sup)
    figs.append(fig);exports+=files
    bundle=FIG/'ch03_supp_du145_metascape_review.pdf'
    with PdfPages(bundle) as pdf:
        for fig in figs:pdf.savefig(fig);plt.close(fig)
    exports.append({'file':bundle.name,'sha256':sha(bundle)})
    provenance={'created_date':'2026-10-05','font':'Open Sans','status':'selected as CH3-S3.4.6.3-S06, S07 and S08 on 2026-10-05',
                'selection':'FDR <0.05, abs(log2FC) >=0.5, matching direction in all three passage pairs',
                'statistics':'Full input GO_AllLists.csv for both network node sizes and top-term bars. Source LogP and Log(q-value) are rounded; displayed q values are approximate. No tests rerun.',
                'network':'Original XGMML term selection and all unique non-self edges retained, including isolated nodes. Louvain seed 42; communities ordered by strongest full-input term. Spring layout seed 42 within fixed display anchors.',
                'recurrence':'Up to 40 genes with at least one membership among the 102 UP exported-network terms, using full-input Hits; ranked by membership count then symbol and grouped by dominant community, tie broken by community ID.',
                'background':{'submitted_unique_symbols':18760,'reported_user_background':16520,'enrichment_table_denominator':15651},
                'runs':records,'exports':exports,'script_sha256':sha(Path(__file__)),
                'inputs':[{'file':str(p),'sha256':sha(p)} for direction in ['up','down'] for p in [SOURCE/direction/'Enrichment_GO/GO_AllLists.csv',SOURCE/direction/'Enrichment_GO/GONetwork.xgmml',THESIS/'research'/f'DU145_{direction.upper()}_BACKGROUND.zip']]}
    (SOURCE/'figure_provenance.json').write_text(json.dumps(provenance,indent=2),encoding='utf-8')
    (SUPP/'provenance.json').write_text(json.dumps(provenance,indent=2),encoding='utf-8')
    print(json.dumps(records,indent=2))
    print('UP genes shown:',len(pd.read_csv(SUPP/'displayed_up_gene_module_counts.csv')))
    print(bundle)


if __name__=='__main__':main()
