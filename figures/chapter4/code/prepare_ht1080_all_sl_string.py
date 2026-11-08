"""Freeze STRING associations among the 48 primary drugZ depleted candidates.

Cached API responses make re-runs local and preserve the initial retrieval.
"""
from pathlib import Path
import sys, json, hashlib
from datetime import datetime, timezone
from urllib.request import Request, urlopen
from urllib.parse import urlencode

THESIS = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(THESIS.parent / '.figure_python_deps'))
import pandas as pd

OUT = THESIS / 'figures/chapter4/source_data/ht1080_all_sl_string'
RAW = OUT / 'raw'
RAW.mkdir(parents=True, exist_ok=True)
SOURCE = THESIS / 'figures/chapter4/source_data/ht1080_candidate_discovery/depleted_candidates.csv'
DZ = pd.read_csv(SOURCE)
assert len(DZ) == 48 and DZ.GENE.is_unique
assert (DZ.fdr_synth < .25).all() and (DZ.normZ < 0).all()
DZ.to_csv(OUT/'input_candidates.csv', index=False)
LOG = []

def query(base, method, params, stem):
    path = RAW / f'{stem}.json'
    meta_path = RAW / f'{stem}_request.json'
    if not path.exists():
        url = f'{base}/api/json/{method}'
        payload = urlencode(params).encode() if params else None
        req = Request(url, data=payload, headers={'User-Agent':'HT1080-thesis-figure/1.0'})
        with urlopen(req, timeout=55) as response:
            data = response.read()
        parsed = json.loads(data)
        if not isinstance(parsed, list):
            raise ValueError(f'Unexpected response for {method}: {parsed}')
        path.write_bytes(data)
        meta_path.write_text(json.dumps({'url':url, 'parameters':params,
            'retrieved_utc':datetime.now(timezone.utc).isoformat()},indent=2), encoding='utf-8')
    data = json.loads(path.read_text(encoding='utf-8'))
    meta = json.loads(meta_path.read_text(encoding='utf-8'))
    LOG.append({**meta, 'file':path.name, 'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
    return data

version = query('https://string-db.org', 'version', {}, 'version')[0]
base = version.get('stable_address', version.get('string_stable_address'))
assert base.startswith('https://version-') and base.endswith('.string-db.org')
common = {'species':9606, 'caller_identity':'HT1080_LOY_thesis'}
mapping = pd.DataFrame(query(base, 'get_string_ids', {**common,
    'identifiers':'\r'.join(DZ.GENE), 'echo_query':1, 'limit':1}, 'mapping'))
assert set(mapping.queryItem) == set(DZ.GENE), 'Unmapped input: review before proceeding'
assert mapping.queryItem.is_unique and mapping.stringId.is_unique
assert (mapping.ncbiTaxonId == 9606).all()
mapping.to_csv(OUT/'identifier_mapping.csv', index=False)
id_to_gene = mapping.set_index('stringId').queryItem.to_dict()
network_params = {**common, 'identifiers':'\r'.join(mapping.stringId),
    'required_score':400, 'add_nodes':0}
networks = {}
for kind in ['functional', 'physical']:
    d = pd.DataFrame(query(base, 'network', {**network_params,'network_type':kind}, f'{kind}_400'))
    assert len(d), f'No edges returned: {kind}'
    assert set(d.stringId_A) | set(d.stringId_B) <= set(mapping.stringId)
    d['ncbiTaxonId'] = pd.to_numeric(d.ncbiTaxonId)
    assert (d.ncbiTaxonId == 9606).all() and d.score.between(.4,1).all()
    d['gene_a'] = d.stringId_A.map(id_to_gene)
    d['gene_b'] = d.stringId_B.map(id_to_gene)
    d['pair'] = ['|'.join(sorted([a,b])) for a,b in zip(d.gene_a,d.gene_b)]
    assert d.pair.is_unique and (d.gene_a != d.gene_b).all()
    d.sort_values(['gene_a','gene_b']).to_csv(OUT/f'{kind}_edges.csv',index=False)
    networks[kind] = d

nodes = DZ.merge(mapping, left_on='GENE', right_on='queryItem',validate='one_to_one')
for kind, edges in networks.items():
    degree = pd.concat([edges.gene_a,edges.gene_b]).value_counts()
    nodes[f'{kind}_degree'] = nodes.GENE.map(degree).fillna(0).astype(int)
nodes.to_csv(OUT/'all_nodes.csv', index=False)
summary = {'string_version':version['string_version'], 'stable_address':base,
    'input_rule':'Primary drugZ depleted candidates: normZ < 0 and fdr_synth < 0.25',
    'input_count':len(nodes), 'mapped_count':len(mapping), 'added_nodes':0,
    'required_score':400, 'score_display_threshold':.4,
    'networks':{k:{'edges':len(e), 'connected_nodes':int((nodes[f'{k}_degree']>0).sum()),
        'isolated_nodes':nodes.loc[nodes[f'{k}_degree']==0,'GENE'].tolist()}
        for k,e in networks.items()},
    'source':str(SOURCE), 'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
    'requests':LOG,
    'interpretation':'Database associations among exploratory screen candidates; not validation of LOY specificity or physical interactions in HT1080. No enrichment test performed.',
    'api_documentation':'https://string-db.org/help/api/'}
(OUT/'provenance.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
print(json.dumps({k:v for k,v in summary.items() if k not in ['requests','source']},indent=2))
