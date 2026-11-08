"""Freeze STRING v12.5 and saved DepMap correlations for 57 drugZ candidates."""
from pathlib import Path
from datetime import datetime, timezone
from urllib.request import Request, urlopen
from urllib.parse import urlencode
import sys, json, hashlib, time, itertools

THESIS = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(THESIS.parent / '.figure_python_deps'))
import pandas as pd
import numpy as np

OUT = THESIS / 'figures/chapter4'
DATA = OUT / 'source_data/ht1080_depmap_string_network_57'
RAW = DATA / 'raw'
RAW.mkdir(parents=True, exist_ok=True)
DZ = THESIS.parent / 'crispr_screen/ht1080_loy_screen_2/results/drugZ/LOY_vs_CTRL_unfiltered.txt'
DEP = OUT / 'source_data/ht1080_depmap_codependency'
BASE = 'https://version-12-5.string-db.org'
LOG = []


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def query(method, parameters, stem):
    path, meta = RAW / f'{stem}.json', RAW / f'{stem}_request.json'
    url = f'{BASE}/api/json/{method}'
    if not path.exists():
        req = Request(url, data=urlencode(parameters).encode(),
                      headers={'User-Agent': 'HT1080-thesis-network/1.0'})
        with urlopen(req, timeout=45) as response:
            payload = response.read()
        assert isinstance(json.loads(payload), list)
        path.write_bytes(payload)
        meta.write_text(json.dumps(dict(url=url, parameters=parameters,
            retrieved_utc=datetime.now(timezone.utc).isoformat()), indent=2), encoding='utf-8')
        time.sleep(1)
    metadata = json.loads(meta.read_text(encoding='utf-8'))
    assert metadata['parameters'] == parameters and metadata['url'] == url
    LOG.append({**metadata, 'file':path.name, 'sha256':sha(path)})
    return json.loads(path.read_text(encoding='utf-8'))


def main():
    scores = pd.read_csv(DZ, sep='\t')
    depleted = scores.loc[(scores.normZ < 0) & (scores.fdr_synth < .25)].copy()
    enriched = scores.loc[(scores.normZ > 0) & (scores.fdr_supp < .25)].copy()
    depleted['direction'], enriched['direction'] = 'depleted', 'enriched'
    nodes = pd.concat([depleted, enriched]).sort_values('normZ', key=abs, ascending=False).reset_index(drop=True)
    nodes['screen_fdr'] = np.where(nodes.direction.eq('depleted'), nodes.fdr_synth, nodes.fdr_supp)
    assert len(nodes) == 57 and len(depleted) == 48 and len(enriched) == 9 and nodes.GENE.is_unique
    genes = nodes.GENE.tolist()
    common = {'species':9606, 'caller_identity':'HT1080_LOY_thesis'}
    mapping = pd.DataFrame(query('get_string_ids', {**common,
        'identifiers':'\r'.join(genes), 'echo_query':1, 'limit':1}, 'mapping_57'))
    assert set(mapping.queryItem) == set(genes)
    assert mapping.queryItem.is_unique and mapping.stringId.is_unique and mapping.ncbiTaxonId.eq(9606).all()
    mapping.to_csv(DATA / 'identifier_mapping.csv', index=False)
    id_to_gene = mapping.set_index('stringId').queryItem.to_dict()
    string = pd.DataFrame(query('network', {**common, 'identifiers':'\r'.join(mapping.stringId),
        'required_score':400, 'add_nodes':0, 'network_type':'functional'}, 'functional_400_57'))
    assert (set(string.stringId_A) | set(string.stringId_B)) <= set(mapping.stringId)
    string['ncbiTaxonId'] = pd.to_numeric(string.ncbiTaxonId)
    assert string.ncbiTaxonId.eq(9606).all() and string.score.between(.4,1).all()
    string['gene_a'] = string.stringId_A.map(id_to_gene)
    string['gene_b'] = string.stringId_B.map(id_to_gene)
    string['pair'] = ['|'.join(sorted([a,b])) for a,b in zip(string.gene_a,string.gene_b)]
    assert string.pair.is_unique and string.gene_a.ne(string.gene_b).all()
    string.sort_values('pair').to_csv(DATA / 'string_functional_edges.csv', index=False)

    r = pd.read_csv(DEP / 'seed_pairwise_pearson.csv', index_col=0).loc[genes, genes]
    n = pd.read_csv(DEP / 'seed_pairwise_n.csv', index_col=0).loc[genes, genes]
    adjusted = pd.read_csv(DEP / 'seed_pairwise_lineage_sex_residual_pearson.csv', index_col=0).loc[genes, genes]
    assert np.allclose(r.to_numpy(), r.to_numpy().T, equal_nan=True)
    assert np.allclose(np.diag(r.drop(index='ZCCHC16', columns='ZCCHC16')), 1)
    r.to_csv(DATA / 'candidate_pearson.csv')
    n.to_csv(DATA / 'candidate_pairwise_n.csv')
    adjusted.to_csv(DATA / 'candidate_lineage_sex_residual_pearson.csv')
    string_score = string.set_index('pair').score.to_dict()
    pairs = []
    for a,b in itertools.combinations(genes,2):
        pair = '|'.join(sorted([a,b]))
        available = pd.notna(r.loc[a,b]) and n.loc[a,b] >= 200
        dep_edge = bool(available and abs(r.loc[a,b]) > .2)
        str_edge = pair in string_score
        cat = 'both' if dep_edge and str_edge else 'depmap_only' if dep_edge else 'string_only' if str_edge else 'neither'
        if not available:
            cat = 'string_with_depmap_unavailable' if str_edge else 'depmap_unavailable'
        pairs.append(dict(gene_a=a,gene_b=b,pair=pair,pearson_r=r.loc[a,b],n_models=n.loc[a,b],
            depmap_available=bool(available),depmap_edge=dep_edge,
            depmap_sign='negative' if dep_edge and r.loc[a,b]<0 else 'positive' if dep_edge else '',
            string_score=string_score.get(pair,np.nan),string_edge=str_edge,evidence_category=cat,
            adjusted_pearson_r=adjusted.loc[a,b]))
    pairs = pd.DataFrame(pairs)
    assert len(pairs)==1596 and pairs.pair.is_unique and pairs.depmap_edge.sum()==41
    pairs.to_csv(DATA / 'pairwise_evidence.csv',index=False)
    pairs.loc[pairs.depmap_edge].to_csv(DATA / 'depmap_edges.csv',index=False)
    nodes = nodes.merge(mapping, left_on='GENE', right_on='queryItem', validate='one_to_one')
    nodes['depmap_available'] = nodes.GENE.ne('ZCCHC16')
    for name, edges in [('depmap',pairs.loc[pairs.depmap_edge]),('string',string)]:
        deg=pd.concat([edges.gene_a,edges.gene_b]).value_counts()
        nodes[name+'_degree']=nodes.GENE.map(deg).fillna(0).astype(int)
    nodes.to_csv(DATA/'candidates.csv',index=False)
    # Check the old 48-candidate subnetwork is preserved by the extended query.
    old = pd.read_csv(OUT/'source_data/ht1080_all_sl_string/functional_edges.csv').set_index('pair')
    current48=string.loc[string.gene_a.isin(depleted.GENE) & string.gene_b.isin(depleted.GENE)].set_index('pair')
    assert set(old.index)==set(current48.index)
    max_diff=float((old.score-current48.score).abs().max())
    assert max_diff<1e-12
    manifest=dict(created_date='2026-10-19',input_count=57,depleted=48,enriched=9,
        input_rule='Original drugZ directional FDR <0.25 and matching normZ direction; no extra connectors',
        depmap_release='Public 25Q3',depmap_cohort_models=1186,depmap_minimum_pairwise_n=200,
        depmap_missing=['ZCCHC16'],depmap_edges=int(pairs.depmap_edge.sum()),
        depmap_positive=int((pairs.depmap_sign=='positive').sum()),depmap_negative=int((pairs.depmap_sign=='negative').sum()),
        depmap_connected=int(nodes.depmap_degree.gt(0).sum()),
        string_version='12.5',string_network_type='functional',string_threshold=.4,added_nodes=0,
        string_mapped=len(mapping),string_edges=len(string),string_connected=int(nodes.string_degree.gt(0).sum()),
        shared_edges=int((pairs.depmap_edge & pairs.string_edge).sum()),
        evidence_category_counts=pairs.evidence_category.value_counts().to_dict(),
        old_48_string_edges_preserved=True,old_48_score_max_absolute_difference=max_diff,
        requests=LOG,sources={str(p):sha(p) for p in [DZ,DEP/'seed_pairwise_pearson.csv',DEP/'seed_pairwise_n.csv',
            DEP/'seed_pairwise_lineage_sex_residual_pearson.csv',DEP/'analysis_manifest.json',Path(__file__)]})
    (DATA/'analysis_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in manifest.items() if k not in ['requests','sources']},indent=2))


if __name__=='__main__':
    main()
