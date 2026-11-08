"""A549 counterpart to HT1080 Figure 2, using the original saved P4 branch.

Screen scores/FDRs and MAGeCK LFCs are read unchanged. New calculations are
descriptive NTC density and public DepMap candidate correlations only.
"""
from pathlib import Path
from datetime import datetime, timezone
from urllib.request import Request, urlopen
from urllib.parse import urlencode
import hashlib
import itertools
import json
import re
import sys
import time

THESIS = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(THESIS.parent / '.figure_python_deps'))
import numpy as np
import pandas as pd

OUT = THESIS / 'figures/chapter4'
DATA = OUT / 'source_data/a549_candidate_discovery'
SCREEN = THESIS.parent / 'crispr_screen/a549_loy_screen1/CRISPR_screen'
DEPMAP = Path('A:/DEPMAP/DepMap')
BASE = 'https://version-12-5.string-db.org'
GAMETOLOGS = ['DDX3X', 'EIF1AX', 'ZFX', 'RPS4X']
ENDPOINTS = ['CTRL_P4_R2', 'LOY_P4_R1', 'LOY_P4_R2']
BASELINES = ['CTRL_D0_1', 'CTRL_D0_2', 'LOY_D0_1', 'LOY_D0_2']
XLIM = (-10.0, 6.25)


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as fh:
        for block in iter(lambda: fh.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def query(method, parameters, stem):
    raw = DATA / 'raw'
    raw.mkdir(parents=True, exist_ok=True)
    path, meta = raw / f'{stem}.json', raw / f'{stem}_request.json'
    url = f'{BASE}/api/json/{method}'
    if not path.exists():
        req = Request(url, data=urlencode(parameters).encode(),
                      headers={'User-Agent': 'A549-thesis-network/1.0'})
        with urlopen(req, timeout=45) as response:
            payload = response.read()
        assert isinstance(json.loads(payload), list)
        path.write_bytes(payload)
        meta.write_text(json.dumps(dict(url=url, parameters=parameters,
            retrieved_utc=datetime.now(timezone.utc).isoformat()), indent=2), encoding='utf-8')
        time.sleep(1)
    metadata = json.loads(meta.read_text(encoding='utf-8'))
    assert metadata['parameters'] == parameters and metadata['url'] == url
    return json.loads(path.read_text(encoding='utf-8')), {
        **metadata, 'file': str(path), 'sha256': sha(path)}


def main():
    DATA.mkdir(parents=True, exist_ok=True)
    paths = {
        'original_drugz': SCREEN / 'drugZ/A549_P4_LOY_vs_CTRL.txt',
        'rerun_drugz': SCREEN / 'results/rerun/drugz/A549_P4_LOY_vs_CTRL_rerun.txt',
        'counts': SCREEN / 'Screen1/Mageck/A549_full.count.txt',
        'mageck_guides': SCREEN / 'Screen1/Mageck/P4_LOY_vs_CTRL.sgrna_summary.txt',
        'mageck_log': SCREEN / 'Screen1/Mageck/P4_LOY_vs_CTRL.log',
        'brunello': SCREEN / 'Screen1/broadgpp-brunello-library-contents.txt',
        'depmap': DEPMAP / 'CRISPRGeneEffect.csv',
        'depmap_readme': DEPMAP / 'README.txt',
    }
    fingerprints = {k: {'path': str(p), 'sha256': sha(p)} for k, p in paths.items()}
    # Ensure the public matrix is exactly the release used for HT1080 Figure 2.
    ht = json.loads((OUT / 'source_data/ht1080_depmap_codependency/analysis_manifest.json').read_text())
    assert fingerprints['depmap']['sha256'] == ht['sources'][str(paths['depmap'])]['sha256']
    dz = pd.read_csv(paths['original_drugz'], sep='\t')
    assert len(dz) == 19112 and dz.GENE.is_unique
    assert not dz.GENE.eq('Non-Targeting-Control').any()
    dz['direction'] = np.where(dz.normZ.lt(0), 'depleted', 'enriched')
    dz['screen_fdr'] = np.where(dz.normZ.lt(0), dz.fdr_synth, dz.fdr_supp)
    dz['display_class'] = np.select([dz.screen_fdr.lt(.05), dz.screen_fdr.lt(.25)],
                                    ['FDR < 0.05', '0.05–0.25'], default='FDR ≥ 0.25')
    dz.to_csv(DATA / 'all_gene_scores.csv', index=False)
    named = dz.loc[dz.screen_fdr.lt(.05)].sort_values('rank_synth').copy()
    nodes = dz.loc[dz.screen_fdr.lt(.25)].sort_values('rank_synth').copy()
    assert len(named) == 19 and len(nodes) == 73 and nodes.direction.eq('depleted').all()
    named.to_csv(DATA / 'main_named_hits.csv', index=False)
    # Preserve the rerun as a comparison, never mix its FDR family into this figure.
    rerun = pd.read_csv(paths['rerun_drugz'], sep='\t')
    comparison = dz[['GENE', 'normZ', 'rank_synth', 'fdr_synth']].merge(
        rerun[['GENE', 'normZ', 'rank_synth', 'fdr_synth']], on='GENE',
        suffixes=('_original', '_rerun'), validate='one_to_one')
    comparison.to_csv(DATA / 'original_rerun_comparison.csv', index=False)

    counts = pd.read_csv(paths['counts'], sep='\t')
    mg = pd.read_csv(paths['mageck_guides'], sep='\t')
    joined = counts.merge(mg[['sgrna', 'Gene', 'LFC', 'FDR']],
        left_on=['sgRNA', 'Gene'], right_on=['sgrna', 'Gene'],
        how='left', validate='one_to_one').drop(columns='sgrna')
    assert len(joined) == 77441 and joined.LFC.notna().all()
    joined['endpoint_read_total'] = joined[ENDPOINTS].sum(axis=1)
    joined['all_retained_samples_zero'] = joined[ENDPOINTS + BASELINES].sum(axis=1).eq(0)
    joined['lfc_displayed'] = joined.endpoint_read_total.gt(0) & np.isfinite(joined.LFC)
    genes = list(dict.fromkeys(GAMETOLOGS + named.GENE.tolist()))
    selected = joined.loc[joined.Gene.isin(genes)].sort_values(['Gene', 'sgRNA']).copy()
    selected['drugz_directional_fdr'] = selected.Gene.map(dz.set_index('GENE').fdr_synth)
    selected['drugz_normZ'] = selected.Gene.map(dz.set_index('GENE').normZ)
    selected.rename(columns={'LFC': 'mageck_reported_log2FC',
                             'FDR': 'mageck_guide_FDR'}, inplace=True)
    assert selected.groupby('Gene').size().eq(4).all()
    selected.to_csv(DATA / 'selected_guides.csv', index=False)
    ntc = joined.loc[joined.Gene.eq('Non-Targeting-Control')].copy()
    assert len(ntc) == 1000 and ntc.lfc_displayed.sum() == 997
    ntc.to_csv(DATA / 'non_targeting_guides.csv', index=False)
    vals = ntc.loc[ntc.lfc_displayed, 'LFC'].to_numpy()
    bandwidth = .9 * min(vals.std(ddof=1), np.subtract(*np.percentile(vals, [75, 25])) / 1.34) * len(vals)**(-.2)
    grid = np.linspace(*XLIM, 3251)
    density = np.exp(-.5 * ((grid[:, None] - vals[None, :]) / bandwidth)**2).mean(axis=1) / (bandwidth * np.sqrt(2*np.pi))
    pd.DataFrame({'guide_log2FC': grid, 'kde_density': density,
                  'display_density_peak_scaled': density / density.max()}).to_csv(DATA / 'ntc_density.csv', index=False)
    assert selected.loc[selected.lfc_displayed, 'mageck_reported_log2FC'].between(*XLIM).all()

    genes = nodes.GENE.tolist()
    print('Retrieving STRING mapping and associations for 73 candidate symbols...', flush=True)
    common = {'species': 9606, 'caller_identity': 'A549_LOY_thesis'}
    response, request1 = query('get_string_ids', {
        **common, 'identifiers': '\r'.join(genes), 'echo_query': 1, 'limit': 1}, 'mapping_73')
    mapping = pd.DataFrame(response)
    assert mapping.queryItem.is_unique and mapping.stringId.is_unique
    assert set(mapping.queryItem) <= set(genes) and mapping.ncbiTaxonId.eq(9606).all()
    # NARR (Brunello Entrez 100861437) is not RAB34 (83871). Do not merge
    # this unresolved target into another gene on the basis of STRING's name hit.
    mapping['accepted'] = ~mapping.queryItem.eq('NARR')
    mapping['mapping_note'] = np.where(mapping.queryItem.eq('NARR'),
        'Rejected: STRING RAB34 differs from the Brunello NARR target Gene ID',
        'Exact symbol' )
    mapping.loc[mapping.queryItem.eq('TMEM136'), 'mapping_note'] = 'TLCD5 alias; same Entrez 219902'
    mapping.to_csv(DATA / 'identifier_mapping.csv', index=False)
    accepted = mapping.loc[mapping.accepted]
    response, request2 = query('network', {
        **common, 'identifiers': '\r'.join(accepted.stringId), 'required_score': 400,
        'add_nodes': 0, 'network_type': 'functional'}, 'functional_400_verified_73')
    st = pd.DataFrame(response)
    id_to_gene = accepted.set_index('stringId').queryItem.to_dict()
    assert (set(st.stringId_A) | set(st.stringId_B)) <= set(accepted.stringId)
    assert st.score.between(.4, 1).all()
    st['gene_a'] = st.stringId_A.map(id_to_gene)
    st['gene_b'] = st.stringId_B.map(id_to_gene)
    st['pair'] = ['|'.join(sorted([a, b])) for a, b in zip(st.gene_a, st.gene_b)]
    assert st.pair.is_unique and st.gene_a.ne(st.gene_b).all()
    st.to_csv(DATA / 'string_functional_edges.csv', index=False)

    print('Reading matching DepMap Public 25Q3 profiles...', flush=True)
    header = pd.read_csv(paths['depmap'], nrows=0).columns.tolist()
    symbol_to_column = {re.sub(r' \(\d+\)$', '', col): col for col in header[1:]}
    assert len(symbol_to_column) == len(header)-1
    library = pd.read_csv(paths['brunello'], sep='\t', dtype=str)
    lib_ids = library[['Target Gene Symbol', 'Target Gene ID']].drop_duplicates()
    assert lib_ids['Target Gene Symbol'].is_unique
    lib_ids = lib_ids.set_index('Target Gene Symbol')['Target Gene ID'].to_dict()
    entrez_to_column = {re.search(r'\((\d+)\)$', col).group(1): col for col in header[1:]}
    mapped_columns, mapping_rows = {}, []
    for g in genes:
        gene_id = lib_ids.get(g)
        col = symbol_to_column.get(g)
        method = 'Exact symbol'
        if col is None and gene_id in entrez_to_column:
            col, method = entrez_to_column[gene_id], 'Brunello target Entrez ID'
        if col is not None:
            if gene_id:
                assert re.search(r'\((\d+)\)$', col).group(1) == gene_id, (g, col, gene_id)
            mapped_columns[g] = col
        mapping_rows.append(dict(GENE=g, brunello_entrez=gene_id,
            depmap_column=col, mapping=method if col else 'Unavailable; no substitution'))
    assert len(set(mapped_columns.values())) == len(mapped_columns)
    pd.DataFrame(mapping_rows).to_csv(DATA / 'depmap_identifier_mapping.csv', index=False)
    present = list(mapped_columns)
    ge = pd.read_csv(paths['depmap'],
        usecols=[header[0]] + list(mapped_columns.values()), index_col=0)
    ge.rename(columns={v: k for k, v in mapped_columns.items()}, inplace=True)
    assert len(ge) == 1186 and ge.index.is_unique
    ge = ge.reindex(columns=genes)
    ge.to_csv(DATA / 'candidate_gene_effect.csv', index_label='ModelID')
    r = ge.corr(method='pearson', min_periods=200)
    valid = ge.notna().astype(int)
    ns = valid.T @ valid
    r.to_csv(DATA / 'candidate_pearson.csv', index_label='gene')
    ns.to_csv(DATA / 'candidate_pairwise_n.csv', index_label='gene')
    st_scores = st.set_index('pair').score
    rows, errors = [], []
    for a, b in itertools.combinations(genes, 2):
        pair = '|'.join(sorted([a, b]))
        available = bool(np.isfinite(r.loc[a, b]) and ns.loc[a, b] >= 200)
        de = bool(available and abs(r.loc[a, b]) > .2)
        se = pair in st_scores
        rows.append(dict(gene_a=a, gene_b=b, pair=pair, pearson_r=r.loc[a, b],
            n_models=int(ns.loc[a, b]), depmap_available=available, depmap_edge=de,
            depmap_sign='positive' if de and r.loc[a, b] > 0 else 'negative' if de else '',
            string_edge=se, string_score=st_scores.get(pair, np.nan)))
        if available:
            arr = ge[[a, b]].dropna().to_numpy()
            assert len(arr) == ns.loc[a, b]
            errors.append(abs(np.corrcoef(arr.T)[0, 1] - r.loc[a, b]))
    pairs = pd.DataFrame(rows)
    assert len(pairs) == len(genes)*(len(genes)-1)//2 and max(errors) < 1e-10
    pairs.to_csv(DATA / 'pairwise_evidence.csv', index=False)
    pairs.loc[pairs.depmap_edge].to_csv(DATA / 'depmap_edges.csv', index=False)
    nodes['depmap_available'] = nodes.GENE.isin(present)
    nodes['string_available'] = nodes.GENE.isin(accepted.queryItem)
    nodes['string_preferred_name'] = nodes.GENE.map(mapping.set_index('queryItem').preferredName)
    nodes['depmap_column'] = nodes.GENE.map(mapped_columns)
    for name, edges in [('depmap', pairs.loc[pairs.depmap_edge]), ('string', st)]:
        deg = pd.concat([edges.gene_a, edges.gene_b]).value_counts()
        nodes[name + '_degree'] = nodes.GENE.map(deg).fillna(0).astype(int)
    nodes.to_csv(DATA / 'candidates.csv', index=False)
    manifest = dict(
        created_utc=datetime.now(timezone.utc).isoformat(),
        status='Primary A549 P4 reporting source retained 2026-10-21; original invocation unresolved',
        primary_account_audit='research/a549_primary_account_20261021/README.md',
        original_invocation_recovered=False,
        primary_branch='original saved P4 drugZ; saved FDRs unchanged',
        samples={'CTRL': ['CTRL_P4_R2'], 'LOY': ['LOY_P4_R1', 'LOY_P4_R2']},
        P3_used=False, D0_used_in_endpoint_scores=False, screen_genes=len(dz),
        strong_depleted=len(named), exploratory_depleted=len(nodes), enriched_FDR25=0,
        named_selection='All directional FDR <0.05 genes; no added pathway-selected genes',
        prior_genes=GAMETOLOGS, guide_unique_genes=selected.Gene.nunique(),
        guide_unique_total=len(selected), guide_unique_displayed=int(selected.lfc_displayed.sum()),
        guide_exclusions=selected.loc[~selected.lfc_displayed, ['sgRNA', 'Gene']].to_dict('records'),
        guide_source='Original saved MAGeCK P4 LOY/CTRL LFC; control-guide normalization; not D0-adjusted',
        ntc_n=997, ntc_bandwidth=float(bandwidth), guide_xlim=XLIM,
        ntc_kde='Gaussian, bandwidth=0.9*min(sample_SD,IQR/1.34)*n^(-1/5); peak-scaled for display',
        depmap_release='Public 25Q3', depmap_models=len(ge),
        depmap_mapping='Exact symbol, then verified Brunello Entrez ID for renamed genes',
        depmap_renamed=[r for r in mapping_rows if r['mapping']=='Brunello target Entrez ID'],
        depmap_missing=nodes.loc[~nodes.depmap_available, 'GENE'].tolist(),
        depmap_edges=int(pairs.depmap_edge.sum()),
        depmap_positive=int(pairs.depmap_sign.eq('positive').sum()),
        depmap_negative=int(pairs.depmap_sign.eq('negative').sum()),
        depmap_rule='Pairwise-complete Pearson, n>=200, |r|>0.2; no significance test',
        string_version='12.5', string_type='functional', string_threshold=.4,
        string_edges=len(st), string_missing=nodes.loc[~nodes.string_available, 'GENE'].tolist(),
        string_renamed=nodes.loc[nodes.string_available & nodes.GENE.ne(nodes.string_preferred_name),
                                ['GENE', 'string_preferred_name']].to_dict('records'),
        string_mapping_rejected=['NARR -> RAB34: different Brunello target gene ID'],
        shared_edges=int((pairs.depmap_edge & pairs.string_edge).sum()),
        no_edge_candidates=int(nodes.depmap_degree.add(nodes.string_degree).eq(0).sum()),
        validation={'correlations_checked':len(errors), 'max_numpy_difference':float(max(errors))},
        caveats=['One CTRL endpoint; no control reproducibility estimate',
                 'Endpoint contrasts can retain D0 population differences',
                 'Public networks are descriptive context, not LOY-specific validation',
                 'Original/rerun distinction retained; no FDR recomputation'],
        requests=[request1, request2], sources=fingerprints, script_sha256=sha(Path(__file__)))
    (DATA / 'analysis_manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    print(json.dumps({k: v for k, v in manifest.items() if k not in ['sources', 'requests']}, indent=2))


if __name__ == '__main__':
    main()
