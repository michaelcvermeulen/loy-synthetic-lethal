"""Build reproducible enrichment maps from frozen fgsea results; no new tests."""
from pathlib import Path
from datetime import datetime, timezone
from itertools import combinations
from collections import Counter
import hashlib
import json

import networkx as nx
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'outputs/tables/gsea_network'
MODEL = ROOT / 'outputs/models/gsea_network'
PROV = ROOT / 'provenance/gsea_network'
CONFIG = ROOT / 'config/gsea_network.json'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def csv(path, **kwargs):
    return pd.read_csv(path, float_precision='round_trip', **kwargs)


def save(data, name):
    data.to_csv(OUT / f'{name}.csv', index=False, float_format='%.17g')


def graph_at(nodes, edges, basis, threshold):
    graph = nx.Graph()
    graph.add_nodes_from(nodes.node_id)
    selected = edges[(edges[f'{basis}_jaccard'] >= threshold) &
                     (edges[f'{basis}_shared_n'] >= 3)]
    graph.add_weighted_edges_from(selected[['source', 'target', f'{basis}_jaccard']].itertuples(index=False, name=None))
    return graph


def main():
    for folder in [OUT, MODEL, PROV]:
        folder.mkdir(parents=True, exist_ok=True)
    cfg = json.loads(CONFIG.read_text())
    paths = [ROOT / f for f in [
        'outputs/tables/fgsea_all_collections.csv',
        'resources/msigdb_all_membership.csv.gz',
        'outputs/tables/pathway_TF_effect_ranking.csv',
        'outputs/tables/camera_all_collections.csv',
        'outputs/tables/camera_estimated_correlation_sensitivity.csv']]
    before = {p.relative_to(ROOT).as_posix(): sha(p) for p in paths}
    fg = csv(paths[0])
    assert not fg.duplicated(['contrast', 'pathway']).any()
    combined = fg[fg.contrast.eq('combined')].copy()
    assert len(combined) == 23696 and combined.test_status.eq('tested').all()
    nodes = combined[combined.BH_across_collections.lt(cfg['global_BH_threshold'])].copy()
    nodes['direction'] = np.where(nodes.NES > 0, 'UP', 'DOWN')
    nodes = nodes.sort_values(['direction', 'pathway']).reset_index(drop=True)
    nodes['node_id'] = nodes.direction.str[0] + (nodes.groupby('direction').cumcount()+1).map(lambda n: f'{n:03d}')
    assert nodes.direction.value_counts().to_dict() == {'DOWN': 776, 'UP': 32}
    rank = csv(paths[2])
    assert rank.gene_name.is_unique and len(rank) == 25461
    universe = set(rank.gene_name)
    member = csv(paths[1])
    member = member[member.pathway.isin(nodes.pathway) & member.symbol.isin(universe)]
    full_by_path = member.groupby('pathway').symbol.agg(set).to_dict()
    full = {}; leading = {}
    for row in nodes.itertuples():
        genes = row.leadingEdge.split(';')
        assert len(set(genes)) == len(genes) and len(genes) > 0
        leading[row.node_id] = set(genes)
        full[row.node_id] = full_by_path[row.pathway]
        assert len(full[row.node_id]) == row.size, row.pathway
        assert leading[row.node_id] <= full[row.node_id] <= universe, row.pathway
    nodes['leading_edge_n'] = nodes.node_id.map(lambda n: len(leading[n]))
    nodes['represented_members'] = nodes.node_id.map(lambda n: ';'.join(sorted(full[n])))
    nodes['neglog10_global_BH'] = -np.log10(nodes.BH_across_collections)
    for period in ['AUG_2025', 'MAY_2026']:
        other = fg[fg.contrast.eq(period)].set_index('pathway')
        nodes[f'{period}_NES'] = nodes.pathway.map(other.NES)
        nodes[f'{period}_global_BH'] = nodes.pathway.map(other.BH_across_collections)
        nodes[f'{period}_same_direction'] = np.sign(nodes[f'{period}_NES']) == np.sign(nodes.NES)
        nodes[f'{period}_support'] = nodes[f'{period}_same_direction'] & nodes[f'{period}_global_BH'].lt(.05)
    a = nodes.AUG_2025_support; m = nodes.MAY_2026_support
    nodes['period_support'] = np.select([a & m, a, m], ['Both periods', 'August only', 'May only'], default='Combined only')
    nodes['both_period_NES_directions_agree'] = nodes.AUG_2025_same_direction & nodes.MAY_2026_same_direction
    for path, prefix in [(paths[3], 'camera_fixed'), (paths[4], 'camera_estimated')]:
        cam = csv(path); cam = cam[cam.contrast.eq('combined')].set_index('pathway')
        assert cam.index.is_unique
        nodes[f'{prefix}_global_BH'] = nodes.pathway.map(cam.BH_across_collections)
        nodes[f'{prefix}_direction'] = nodes.pathway.map(cam.Direction)
        nodes[f'{prefix}_same_direction_BH05'] = nodes[f'{prefix}_global_BH'].lt(.05) & (
            nodes[f'{prefix}_direction'].str.upper() == nodes.direction)
    assert not nodes.camera_estimated_global_BH.lt(.05).any()
    # Store every within-direction pair sharing a gene in either representation.
    records = []
    for direction, subset in nodes.groupby('direction', sort=True):
        for u, v in combinations(subset.node_id, 2):
            inter = full[u] & full[v]
            if not inter:
                continue
            row = {'source': u, 'target': v, 'direction': direction}
            for basis, sets in [('leading_edge', leading), ('full_set', full)]:
                shared = sets[u] & sets[v]
                row[f'{basis}_shared_n'] = len(shared)
                row[f'{basis}_jaccard'] = len(shared)/len(sets[u] | sets[v])
                row[f'{basis}_overlap_coefficient'] = len(shared)/min(len(sets[u]), len(sets[v]))
                # Full membership lists in nodes allow recovery of all shared genes.
                if basis == 'leading_edge':
                    row['shared_leading_edge_genes'] = ';'.join(sorted(shared))
            records.append(row)
    edges = pd.DataFrame(records)
    edges['leading_edge_primary_edge'] = edges.leading_edge_jaccard.ge(cfg['jaccard_threshold']) & edges.leading_edge_shared_n.ge(3)
    edges['full_set_primary_edge'] = edges.full_set_jaccard.ge(cfg['jaccard_threshold']) & edges.full_set_shared_n.ge(3)
    edges.to_csv(OUT / 'all_nonzero_gene_overlap_pairs.csv.gz', index=False, float_format='%.17g', compression={'method':'gzip','mtime':0})
    save(edges[edges.leading_edge_primary_edge | edges.full_set_primary_edge], 'network_edges')
    communities = []; assignments = {}; degrees = {}; sensitivity = []; drivers = []
    stats = {}
    rank_effect = rank.set_index('gene_name').combined.to_dict()
    for direction, subset in nodes.groupby('direction', sort=True):
        pairs = edges[edges.direction.eq(direction)]
        primary = graph_at(subset, pairs, 'leading_edge', cfg['jaccard_threshold'])
        groups = nx.community.louvain_communities(primary, weight='weight', resolution=cfg['louvain_resolution'], seed=cfg['seed'])
        groups = sorted(groups, key=lambda g: (-len(g), min(g)))
        metadata = subset.set_index('node_id')
        for i, group in enumerate(groups, 1):
            cid = f'{direction[0]}C{i:02d}'
            for node in group:
                assignments[node] = cid; degrees[node] = primary.degree(node)
            sg = primary.subgraph(group)
            representative = sorted(group, key=lambda n: (-sg.degree(n, weight='weight'), metadata.loc[n, 'BH_across_collections'], n))[0]
            counts = Counter(g for n in group for g in leading[n])
            top = sorted(counts, key=lambda g: (-counts[g], g))
            terms = metadata.loc[sorted(group)].sort_values('BH_across_collections').pathway.tolist()
            communities.append({'community_id': cid, 'direction': direction, 'pathways': len(group),
                'internal_edges': sg.number_of_edges(), 'representative_node': representative,
                'representative_pathway': metadata.loc[representative, 'pathway'],
                'minimum_global_BH': metadata.loc[sorted(group), 'BH_across_collections'].min(),
                'top_shared_leading_genes': ';'.join(f'{g}:{counts[g]}' for g in top[:15]),
                'leading_edge_union_n': len(counts), 'top_pathways_by_BH': ';'.join(terms[:12]),
                'all_pathways': ';'.join(terms)})
            for gene in top:
                drivers.append({'community_id': cid, 'direction': direction, 'gene': gene,
                    'leading_edge_pathways': counts[gene], 'community_pathways': len(group),
                    'pathway_fraction': counts[gene]/len(group), 'combined_log2_difference': rank_effect[gene]})
        for basis in ['leading_edge', 'full_set']:
            for threshold in cfg['sensitivity_jaccard_thresholds']:
                graph = graph_at(subset, pairs, basis, threshold)
                comps = list(nx.connected_components(graph))
                sensitivity.append({'direction':direction,'gene_basis':basis,'jaccard_threshold':threshold,
                    'minimum_shared_genes':3,'nodes':len(graph),'edges':graph.number_of_edges(),
                    'isolates':len(list(nx.isolates(graph))),'connected_components':len(comps),
                    'largest_component':max(map(len,comps))})
                if threshold == cfg['jaccard_threshold']:
                    attrs = metadata.fillna('').to_dict('index')
                    for n in graph:
                        graph.nodes[n].update({k:(v.item() if hasattr(v,'item') else v) for k,v in attrs[n].items()})
                        graph.nodes[n]['leading_edge_community'] = assignments[n]
                    nx.write_graphml(graph, MODEL / f'{direction}_{basis}.graphml')
        stats[direction] = {'nodes':len(primary),'leading_edge_edges':primary.number_of_edges(),
            'leading_edge_isolates':len(list(nx.isolates(primary))), 'communities':len(groups),
            'non_singleton_communities':sum(len(g)>1 for g in groups),
            'period_support':subset.period_support.value_counts().to_dict(),
            'both_period_directions_agree':int(subset.both_period_NES_directions_agree.sum())}
    nodes['community_id'] = nodes.node_id.map(assignments)
    nodes['leading_edge_degree'] = nodes.node_id.map(degrees)
    save(nodes, 'network_nodes'); save(pd.DataFrame(communities), 'communities')
    save(pd.DataFrame(drivers), 'community_gene_drivers'); save(pd.DataFrame(sensitivity), 'overlap_threshold_sensitivity')
    # Preserve the original global correction family; no selected-node readjustment.
    source = combined.set_index('pathway').loc[nodes.pathway]
    assert np.array_equal(source.BH_across_collections.to_numpy(), nodes.BH_across_collections.to_numpy())
    assert np.array_equal(source.NES.to_numpy(), nodes.NES.to_numpy())
    assert before == {p.relative_to(ROOT).as_posix():sha(p) for p in paths}
    result = {'created_utc':datetime.now(timezone.utc).isoformat(), 'config':cfg,
        'input_sha256':before,'script_sha256':sha(__file__),'original_combined_tests':len(combined),
        'ranked_gene_universe':len(universe),'new_enrichment_tests':False,'summary':stats,
        'versions':{'networkx':nx.__version__,'pandas':pd.__version__,'numpy':np.__version__}}
    (PROV/'analysis.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(stats,indent=2))
    print(pd.DataFrame(communities).query('pathways >= 7')[['community_id','pathways','representative_pathway']].to_string(index=False))


if __name__ == '__main__':
    main()
