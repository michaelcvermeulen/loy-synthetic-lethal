"""Check HDAC figure semantics and standalone report against complete source tables."""
import csv
import hashlib
import json
from pathlib import Path
from urllib.parse import unquote, urlparse

import pandas as pd
from bs4 import BeautifulSoup

root = Path(__file__).resolve().parents[1]
table = lambda name: pd.read_csv(root / 'outputs/tables' / (name + '.csv'))
src = lambda name: pd.read_csv(root / 'outputs/source_data' / (name + '.csv'))
labels = {
    'HDAC2_only_vs_HDAC1': 'HDAC2 only',
    'HDAC1_only_vs_HDAC2': 'HDAC1 only',
    'concordant_shared': 'Shared, concordant',
    'double_knockdown_only': 'Double KD only',
    'conflicting_direction': 'Conflicting direction',
}
panel = src('Fig08_HDAC12_gene_drivers_a')
assert panel['label'].notna().all()
assert panel['label'].equals(panel['partition'].map(labels)), 'Membership labels are misassigned'
members = table('HDAC12_membership_partitions')
for row in panel.itertuples():
    expected = members.represented & members.partition.eq(row.partition) & members.partition_arm.eq(row.source_arm)
    assert expected.sum() == row.genes
sample = src('SuppFig05_HDAC12_scores_processes_a')
assert len(sample) == 48 and sample.sample_key.nunique() == 12
assert sample.groupby(['source', 'batch', 'treatment']).size().eq(3).all()
score = table('HDAC12_percentile_score_contrasts')
matched = src('Fig08_HDAC12_gene_drivers_b')
assert len(matched) == 24 and matched.contrast.nunique() == 6
check = matched.merge(score, on=['source', 'score_arm', 'contrast'], suffixes=('_plot', '_table'), validate='one_to_one')
assert (check.percentile_difference_plot - check.percentile_difference_table).abs().max() < 1e-12
readouts = table('HDAC12_readout_shortlist')
heat = src('Fig08_HDAC12_gene_drivers_c')
assert len(heat) == len(readouts) * 6
for row in heat.itertuples():
    gene = readouts[readouts.gene_name.eq(row.gene_name)].iloc[0]
    assert abs(gene[row.block] - row.difference) < 1e-12
assert set(readouts.source_arm) == {'UP', 'DN'}
assert all(g not in readouts.gene_name.tolist() for g in ['HDAC1', 'HDAC2', 'HDAC3'])

report = root / 'outputs/report/HDAC12_gene_drivers.html'
soup = BeautifulSoup(report.read_text(encoding='utf-8'), 'html.parser')
broken = []
for a in soup.find_all('a', href=True):
    href = a['href']
    if href.startswith(('#', 'http:', 'https:', 'mailto:', 'data:', 'javascript:')):
        continue
    path = unquote(urlparse(href).path)
    if path and not (report.parent / path).exists():
        broken.append(href)
assert not broken, broken
images = soup.find_all('img')
assert len(images) == 2 and all(i.get('src', '').startswith('data:') for i in images)
assert len(soup.find_all('table')) == 4
assert '94 HDAC2-only DOWN genes' in soup.get_text(' ', strip=True)

rows = list(csv.DictReader((root / 'outputs/figures/figure_manifest.csv').open(encoding='utf-8')))
unique = {row['figure']: row for row in rows}
new = {'Fig08_HDAC12_gene_drivers', 'SuppFig05_HDAC12_scores_processes',
       'Fig09_HDAC_external_comparison', 'SuppFig06_HDAC_external_QC_sensitivity',
       'Fig10_PRC2_gene_programs', 'Fig11_epigenetic_context',
       'Fig12_GSEA_PRC2_H3K27me3', 'Fig13_GSEA_other_chromatin'}
original = {k: v for k, v in unique.items() if k not in new}
assert len(original) == 11
for row in original.values():
    for kind in ['pdf', 'svg', 'png']:
        actual = hashlib.sha256((root / row[kind]).read_bytes()).hexdigest()
        assert actual == row[kind + '_sha256'], 'An earlier figure changed: ' + row[kind]
result = dict(membership_labels_verified=True, library_points=len(sample), timing_match_points=len(matched),
    readout_heatmap_cells=len(heat), standalone_figures=len(images), standalone_tables=len(soup.find_all('table')),
    broken_local_links=broken, unchanged_earlier_figure_exports=33)
(root / 'provenance/HDAC12_step2_deliverable_checks.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
print('Passed HDAC figure labels and source values, all libraries/matches, standalone report, and 33 unchanged prior exports.')
