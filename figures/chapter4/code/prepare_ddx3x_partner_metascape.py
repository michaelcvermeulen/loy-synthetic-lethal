"""Import the 31-partner Metascape report; preserve original full-list statistics."""
from pathlib import Path, PurePosixPath
from html.parser import HTMLParser
import hashlib
import json
import shutil
import stat
import zipfile
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd

THESIS = Path(__file__).resolve().parents[3]
OUT = THESIS / 'figures/chapter4'
DATA = OUT / 'source_data/ddx3x_partner_metascape'
ROOT = THESIS / 'research/chapter4_metascape_ddx3x_partners_2026-10-20'
REPORT = ROOT / 'report'
DOWNLOAD = Path('C:/Users/micha/Downloads/all.tw0djeyyo.zip')
ORA = OUT / 'source_data/ddx3x_partner_ora'
PATTERN = OUT / 'source_data/ddx3x_genomewide_y_pattern'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


class ReportText(HTMLParser):
    def __init__(self):
        super().__init__()
        self.parts = []

    def handle_data(self, value):
        if value.strip():
            self.parts.append(value.strip())


def import_report():
    ROOT.mkdir(parents=True, exist_ok=True)
    archived = ROOT / DOWNLOAD.name
    if not archived.exists():
        shutil.copy2(DOWNLOAD, archived)
    elif DOWNLOAD.exists():
        assert sha(archived) == sha(DOWNLOAD)
    entries, seen = [], set()
    with zipfile.ZipFile(archived) as z:
        assert sum(i.file_size for i in z.infolist()) < 100_000_000
        for item in z.infolist():
            member = PurePosixPath(item.filename)
            assert not member.is_absolute() and '..' not in member.parts
            assert '\\' not in item.filename and ':' not in item.filename
            assert not stat.S_ISLNK(item.external_attr >> 16)
            assert item.filename.casefold() not in seen
            seen.add(item.filename.casefold())
            target = REPORT.joinpath(*member.parts).resolve()
            assert target.is_relative_to(REPORT.resolve())
            entries.append((item, target))
        for item, target in entries:
            if item.is_dir():
                target.mkdir(parents=True, exist_ok=True)
                continue
            payload = z.read(item)
            target.parent.mkdir(parents=True, exist_ok=True)
            if target.exists():
                assert target.read_bytes() == payload
            else:
                target.write_bytes(payload)
    record = dict(import_date='2026-10-20', original_download=str(DOWNLOAD),
                  archived_zip=str(archived), sha256=sha(archived),
                  entries=len(entries), paths_and_symlinks_validated=True,
                  files=[dict(path=str(p.relative_to(ROOT)), bytes=p.stat().st_size,
                              sha256=sha(p)) for i, p in entries if not i.is_dir()])
    (ROOT / 'import_manifest.json').write_text(json.dumps(record, indent=2),
                                              encoding='utf-8')
    return record


def parse_network(name):
    root = ET.parse(REPORT / 'Enrichment_PPI/xgmml' / name).getroot()
    nodes, edges, id_map = [], [], {}
    for node in root.findall('{*}node'):
        attrs = {a.attrib['name']: a.attrib.get('value') for a in node.findall('{*}att')}
        symbol = attrs['Symbol']
        id_map[node.attrib['id']] = symbol
        graphics = node.find('{*}graphics')
        nodes.append(dict(symbol=symbol, entrez_id=int(attrs['Gene']),
                          cluster=int(attrs.get('MCODE_CLUSTER_ID', 0)),
                          x=float(graphics.attrib['x']), y=float(graphics.attrib['y'])))
    for edge in root.findall('{*}edge'):
        attrs = {a.attrib['name']: a.attrib.get('value') for a in edge.findall('{*}att')}
        edges.append(dict(source=id_map[edge.attrib['source']],
                          target=id_map[edge.attrib['target']],
                          score=attrs.get('SCORE'), support=attrs.get('support')))
    return pd.DataFrame(nodes), pd.DataFrame(edges)


def main():
    DATA.mkdir(parents=True, exist_ok=True)
    imported = import_report()
    expected = set((ORA / 'query_31.txt').read_text().splitlines())
    annotation = pd.read_excel(REPORT / 'metascape_result.xlsx', sheet_name='Annotation')
    assert len(annotation) == len(expected) == 31
    assert annotation.MyList.is_unique and annotation['Gene ID'].is_unique
    assert set(annotation.MyList) == expected
    selected = pd.read_csv(PATTERN / 'selected_partners.csv')
    assert set(selected.gene) == expected
    annotation.merge(selected, left_on='MyList', right_on='gene',
                     validate='one_to_one').to_csv(DATA / 'input_31_verified.csv', index=False)
    original = pd.read_csv(REPORT / 'Enrichment_GO/GO_AllLists.csv')
    final = pd.read_csv(REPORT / 'Enrichment_GO/_FINAL_GO.csv')
    assert len(original) == 291 and original.GO.is_unique and final.GO.is_unique
    assert set(original.GO) == set(final.GO)
    assert original['#GeneInHitList'].eq(31).all()
    assert original['#TotalGeneInLibrary'].eq(17768).all()
    gene_map = dict(zip(annotation['Gene ID'].astype(str), annotation.MyList))
    original['source_row'] = np.arange(len(original))
    original['input_gene_symbols'] = original.GeneID.apply(
        lambda value: '|'.join(gene_map[x] for x in value.split('|')))
    assert original.apply(lambda r: set(r.input_gene_symbols.split('|')) ==
                           set(r.Hits.split('|')), axis=1).all()
    assert original.apply(lambda r: len(r.input_gene_symbols.split('|')) ==
                           r['#GeneInGOAndHitList'], axis=1).all()
    original['approx_q_value'] = 10. ** original['Log(q-value)']
    original['q_boundary_rounding_ambiguous'] = original['Log(q-value)'].eq(-1.3)
    original['reported_q_clearly_lt_005'] = original['Log(q-value)'].lt(-1.35)
    original = original.merge(final[['GO', 'GROUP_ID']], on='GO', validate='one_to_one')
    original.to_csv(DATA / 'original_full31_enrichment.csv', index=False)
    reps = original.sort_values(['LogP', 'source_row'], kind='stable').drop_duplicates(
        'GROUP_ID').copy()
    assert len(reps) == 14 and not reps.q_boundary_rounding_ambiguous.any()
    reps.to_csv(DATA / 'cluster_representatives.csv', index=False)
    audit = original[['GO', 'Description', 'LogP', 'Log(q-value)', '#GeneInHitList',
                       '#GeneInGOAndHitList']].merge(
        final[['GO', 'LogP', 'Log(q-value)', '#GeneInHitList',
               '#GeneInGOAndHitList', 'EvidenceCutoff']],
        on='GO', suffixes=('_original_full31', '_final_prioritised'),
        validate='one_to_one')
    audit.to_csv(DATA / 'original_vs_prioritised_audit.csv', index=False)
    modules = pd.read_csv(REPORT / 'Enrichment_PPI/MCODE.csv')
    assert set(modules.Symbol) <= expected and modules.Symbol.is_unique
    modules.to_csv(DATA / 'mcode_members.csv', index=False)
    network_counts = {}
    for name, stem in [('MyList.xgmml', 'ppi'), ('MyList_MCODE_ALL.xgmml', 'mcode')]:
        nodes, edges = parse_network(name)
        assert set(nodes.symbol) <= expected
        assert set(edges.source) | set(edges.target) <= set(nodes.symbol)
        nodes.to_csv(DATA / f'{stem}_nodes.csv', index=False)
        edges.to_csv(DATA / f'{stem}_edges.csv', index=False)
        network_counts[stem] = {'nodes': len(nodes), 'edges': len(edges)}
    nodes = pd.read_csv(DATA / 'ppi_nodes.csv')
    pd.DataFrame({'gene': sorted(expected - set(nodes.symbol))}).to_csv(
        DATA / 'genes_without_exported_ppi_connection.csv', index=False)
    module_summary = modules.groupby('Cluster').agg(
        genes=('Symbol', lambda s: '; '.join(sorted(s))),
        n_genes=('Symbol', 'size'), score=('Score', 'first')).reset_index()
    module_summary.to_csv(DATA / 'mcode_summary.csv', index=False)
    text = ReportText()
    text.feed((REPORT / 'AnalysisReport.html').read_text(encoding='utf-8'))
    report_text = '\n'.join(text.parts)
    assert '17786 genes' in report_text and '2026-07-01' in report_text
    (DATA / 'report_text.txt').write_text(report_text, encoding='utf-8')
    record = dict(
        analysis_date='2026-10-20', archive_sha256=imported['sha256'],
        archive=str(ROOT / DOWNLOAD.name), input_genes=31, exact_input_match=True,
        input_unique_entrez=31, report_custom_background_genes=17786,
        local_supplied_background_genes=len((ORA / 'background_17786.txt').read_text().splitlines()),
        effective_table_background_genes=17768,
        background_identity_limit='Report confirms submitted count, but does not export '
                                  'the full mapped background or identify the 18-gene difference.',
        annotation_database_date='2026-07-01',
        original_reported_terms=len(original), displayed_cluster_representatives=len(reps),
        original_terms_clearly_q_lt_005=int(original.reported_q_clearly_lt_005.sum()),
        original_terms_q_boundary_ambiguous=int(original.q_boundary_rounding_ambiguous.sum()),
        representative_terms_clearly_q_lt_005=int(reps.reported_q_clearly_lt_005.sum()),
        inference='Reported enrichment statistics are retained. Original full-list terms '
                  'are a filtered export, not all hypotheses tested; do not recompute BH '
                  'on this table. q values are approximate due to exported log rounding.',
        statistical_source='GO_AllLists.csv, all 31 input genes',
        final_table_use='GROUP_ID only; prioritised-subset statistics are not used in the figure.',
        display_selection='Smallest original full-list LogP in each final-table GROUP_ID; '
                          'source order resolves rounded ties; all 14 groups retained.',
        network_counts=network_counts, mcode=module_summary.to_dict('records'),
        sources={str(p): sha(p) for p in [
            Path(__file__), REPORT / 'metascape_result.xlsx',
            REPORT / 'Enrichment_GO/GO_AllLists.csv',
            REPORT / 'Enrichment_GO/_FINAL_GO.csv',
            REPORT / 'Enrichment_PPI/MCODE.csv', REPORT / 'AnalysisReport.html',
            ORA / 'query_31.txt', ORA / 'background_17786.txt',
            PATTERN / 'selected_partners.csv']})
    (DATA / 'manifest.json').write_text(json.dumps(record, indent=2), encoding='utf-8')
    print(json.dumps({k: v for k, v in record.items() if k != 'sources'}, indent=2))
    print(reps[['Description', 'approx_q_value', '#GeneInGOAndHitList']].to_string(index=False))


if __name__ == '__main__':
    main()
