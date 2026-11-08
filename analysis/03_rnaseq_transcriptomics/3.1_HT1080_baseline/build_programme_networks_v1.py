"""New gene selection and programme networks; preserve all previous analyses."""
from pathlib import Path
import json, hashlib
from itertools import combinations
from collections import Counter
import numpy as np
import pandas as pd
import networkx as nx
import network_modules_analysis as legacy

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'outputs/programme_analysis_v1';T=OUT/'tables'
N=T/'STRING';N.mkdir(exist_ok=True)
def save(d,name): d.to_csv(T/f'{name}.csv',index=False)
def gene_sets(v): return set(str(v).split(';')) if pd.notna(v) and str(v) else set()

def string_modules():
    cfg=json.loads((ROOT/'config/network_modules.json').read_text())
    allgenes=pd.read_csv(T/'gene_consistency.csv')
    reps=pd.read_csv(ROOT/'outputs/tables/symbol_mapping_representatives.csv')
    data=allgenes[allgenes.gene_id.isin(reps.gene_id)].copy()
    assert data.gene_name.is_unique and len(data)==25461
    info=pd.read_csv(legacy.STRING/'9606.protein.info.v12.0.txt.gz',sep='\t')
    data['string_id']=data.gene_name.map(info.set_index('preferred_name')['#string_protein_id'])
    data['direction']=np.where(data.logFC>0,'UP','DOWN')
    # Adapter name required by the reused annotation function; selection is new.
    data['descriptive_candidate']=data.directionally_consistent
    data.to_csv(N/'all_symbol_mapping_and_effects.csv',index=False)
    candidates=data[data.directionally_consistent & data.string_id.notna()]
    symbols=set(candidates.gene_name)
    print(f'STRING: {len(candidates)} mapped candidates',flush=True)
    edge=legacy.load_edges(info,symbols,'9606.protein.links.v12.0.txt.gz',700)
    physical=legacy.load_edges(info,symbols,'9606.protein.physical.links.v12.0.txt.gz',700)
    edge.to_csv(N/'all_candidate_association_edges.csv',index=False)
    physical.to_csv(N/'all_candidate_physical_edges.csv',index=False)
    modules=[];coverage=[];node=[];strict=[]
    for direction in ['DOWN','UP']:
        chosen=set(candidates.loc[candidates.direction==direction,'gene_name'])
        graph=nx.Graph();graph.add_nodes_from(sorted(chosen))
        ed=edge[edge.gene_a.isin(chosen)&edge.gene_b.isin(chosen)]
        graph.add_weighted_edges_from(ed[['gene_a','gene_b','score']].itertuples(index=False,name=None))
        cls,scores=legacy.mcode(graph,cfg['mcode'])
        assigned={}
        for i,c in enumerate(cls,1):
            mid=f'{direction}_{i:02d}';c.update(module_id=mid,direction=direction)
            modules.append(c)
            for gene in c['members']: assigned[gene]=mid
        for gene in graph:
            node.append(dict(gene_name=gene,direction=direction,module_id=assigned.get(gene,''),degree=graph.degree(gene)))
        nx.write_graphml(graph,OUT/'models'/f'STRING_{direction}.graphml')
        pg=nx.Graph();pg.add_nodes_from(sorted(chosen))
        pe=physical[physical.gene_a.isin(chosen)&physical.gene_b.isin(chosen)]
        pg.add_edges_from(pe[['gene_a','gene_b']].itertuples(index=False,name=None))
        pcs,_=legacy.mcode(pg,cfg['mcode'])
        for i,c in enumerate(pcs,1):strict.append({**c,'direction':direction,'module_id':f'{direction}_PHYSICAL_{i:02d}'})
        coverage.append(dict(direction=direction,selected_symbols=int((data.directionally_consistent&data.direction.eq(direction)).sum()),
          mapped=len(chosen),connected=sum(graph.degree(g)>0 for g in graph),edges=graph.number_of_edges(),modules=len(cls),
          physical_edges=pg.number_of_edges(),physical_modules=len(pcs)))
    legacy.TABLES=N
    if modules:
        _,top,_=legacy.annotate_modules(modules,data,cfg)
        md=pd.DataFrame(modules);md['members']=md.members.map(lambda x:';'.join(x))
        md=md.merge(top[['module_id','pathway','BH_all_modules_and_terms','candidate_background_BH']],on='module_id',how='left')
        md.to_csv(N/'modules.csv',index=False)
    if strict:
        sd=pd.DataFrame(strict);sd['members']=sd.members.map(lambda x:';'.join(x));sd.to_csv(N/'physical_modules.csv',index=False)
    pd.DataFrame(node).to_csv(N/'nodes.csv',index=False);pd.DataFrame(coverage).to_csv(N/'coverage.csv',index=False)
    print(json.dumps(coverage),flush=True)

def programme_networks():
    p=pd.read_csv(T/'programme_consistency.csv')
    p=p[p.directionally_consistent].copy()
    p['direction']=np.where(p.NES>0,'UP','DOWN')
    p=p.sort_values(['direction','pathway']).reset_index(drop=True)
    p['node_id']=[f'P{i:04d}' for i in range(1,len(p)+1)]
    genes={r.node_id:gene_sets(r.leadingEdge) for r in p.itertuples()}
    membership=pd.read_csv(ROOT/'resources/msigdb_all_membership.csv.gz')
    universe=set(pd.read_csv(T/'ranking_overall.csv').gene_name)
    membership=membership[membership.pathway.isin(p.pathway)&membership.symbol.isin(universe)]
    full=membership.groupby('pathway').symbol.agg(set).to_dict()
    full={r.node_id:full[r.pathway] for r in p.itertuples()}
    all_edges=[];communities=[];assign={};drivers=[];sens=[]
    for direction,group in p.groupby('direction'):
        for u,v in combinations(group.node_id,2):
            shared=genes[u]&genes[v];fshared=full[u]&full[v]
            if not fshared:continue
            all_edges.append(dict(source=u,target=v,direction=direction,shared_n=len(shared),
              jaccard=len(shared)/len(genes[u]|genes[v]),shared_genes=';'.join(sorted(shared)),
              full_shared_n=len(fshared),full_jaccard=len(fshared)/len(full[u]|full[v])))
    ed=pd.DataFrame(all_edges,columns=['source','target','direction','shared_n','jaccard','shared_genes','full_shared_n','full_jaccard'])
    for direction,group in p.groupby('direction'):
        for basis in ['leading','full']:
            jac='jaccard' if basis=='leading' else 'full_jaccard'
            sh='shared_n' if basis=='leading' else 'full_shared_n'
            for cutoff in [.2,.25,.35,.5]:
                g=nx.Graph();g.add_nodes_from(group.node_id)
                se=ed[(ed.direction==direction)&(ed[jac]>=cutoff)&(ed[sh]>=3)]
                g.add_weighted_edges_from(se[['source','target',jac]].itertuples(index=False,name=None))
                cs=nx.community.louvain_communities(g,weight='weight',resolution=1,seed=42) if g.number_of_edges() else [{n} for n in g]
                cs=sorted(cs,key=lambda a:(-len(a),min(a)))
                sens.append(dict(direction=direction,basis=basis,jaccard=cutoff,nodes=len(g),edges=g.number_of_edges(),communities=len(cs),isolates=len(list(nx.isolates(g)))))
                if basis!='leading' or cutoff!=.25:continue
                nx.write_graphml(g,OUT/'models'/f'programme_{direction}.graphml')
                for i,c in enumerate(cs,1):
                    mid=f'{direction}_C{i:02d}';gg=group[group.node_id.isin(c)].sort_values(['FDR','pathway'])
                    count=Counter(gene for n in c for gene in genes[n])
                    communities.append(dict(community=mid,direction=direction,pathways=len(c),
                      representative_pathway=gg.iloc[0].pathway,min_FDR=gg.FDR.min(),both_periods_supported=int(gg.both_periods_FDR05.sum()),
                      members=';'.join(gg.pathway),leading_genes=';'.join(sorted(count))))
                    for n in c:assign[n]=mid
                    for gene,n in count.items():drivers.append(dict(community=mid,gene_name=gene,contributing_pathways=n))
    p['community']=p.node_id.map(assign)
    save(p,'programme_network_nodes');save(ed,'programme_network_edges');save(pd.DataFrame(communities),'programme_communities')
    save(pd.DataFrame(drivers),'programme_gene_drivers');save(pd.DataFrame(sens),'programme_network_sensitivity')
    # Link inferred factors to programme genes; overlap is evidence mapping, not a new test.
    tf=pd.read_csv(T/'TF_consistency.csv');tf=tf[tf.directionally_consistent]
    net=pd.read_csv(ROOT/'resources/collectri_network.csv')[['source','target','mor']].drop_duplicates()
    net=net[~net.duplicated(['source','target'],keep=False)&np.isfinite(net.mor)]
    de=pd.read_csv(T/'gene_consistency.csv').set_index('gene_name')
    rows=[]
    for c in communities:
        geneset=gene_sets(c['leading_genes'])
        for r in tf.itertuples():
            hits=net[(net.source==r.source)&net.target.isin(geneset)]
            if len(hits):rows.append(dict(community=c['community'],TF=r.source,TF_score=r.score,TF_FDR=r.FDR,
              target_overlap=len(hits),target_genes=';'.join(sorted(hits.target))))
    save(pd.DataFrame(rows,columns=['community','TF','TF_score','TF_FDR','target_overlap','target_genes']),'TF_programme_links')
    print(f'Programme networks: {len(p)} pathways, {len(communities)} communities, {len(rows)} TF/theme links',flush=True)

if __name__=='__main__':
    import sys
    if '--string-only' in sys.argv:string_modules()
    elif '--programmes-only' in sys.argv:programme_networks()
    else:string_modules();programme_networks()
