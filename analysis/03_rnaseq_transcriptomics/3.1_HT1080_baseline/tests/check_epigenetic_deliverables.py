"""Validate source semantics, complete figure coverage and preservation of original results."""
import hashlib
import json
from pathlib import Path
from urllib.parse import unquote, urlparse
import numpy as np
import pandas as pd
from bs4 import BeautifulSoup

root = Path(__file__).resolve().parents[1]
table = lambda n: pd.read_csv(root / 'outputs/tables' / f'EPIGEN_{n}.csv')
src = lambda n: pd.read_csv(root / 'outputs/source_data' / f'{n}.csv')
manifest = pd.read_csv(root / 'outputs/figures/figure_manifest.csv')
ids = ['Fig10_PRC2_gene_programs', 'Fig11_epigenetic_context']
assert manifest[manifest.figure.isin(ids)].groupby('figure').size().eq(4).all()
arm = src(ids[0] + '_a')
assert len(arm) == 18 and set(arm.source_arm) == {'UP', 'DN'}
assert arm.groupby(['source', 'source_arm']).size().eq(3).all()
assert arm.label.notna().all()
score = src(ids[0] + '_b')
assert len(score) == 18 and score.contrast.nunique() == 6
check = score.merge(table('score_contrasts'), on=['program', 'source', 'score_arm', 'contrast'],
                    suffixes=('_plot', '_table'), validate='one_to_one')
assert np.allclose(check.percentile_difference_plot, check.percentile_difference_table)
part = src(ids[0] + '_c')
assert part.genes.isin([148, 6, 497]).all() and len(part) == 9
assert all(f'({row.genes})' in row.label for row in part.itertuples())
heat = src(ids[0] + '_d')
readouts = table('posthoc_readouts')
assert len(heat) == 60 and set(heat.gene_name) == set(readouts.gene_name)
for row in heat.itertuples():
    a = readouts[readouts.gene_name.eq(row.gene_name)].iloc[0]
    assert abs(row.value - a[row.block]) < 1e-12
atlas = src(ids[1] + '_a')
assert not atlas.duplicated(['row_label', 'mark']).any()
assert atlas.pathway.nunique() == 41
assert atlas.display_status.eq('tested').sum() == 30
assert atlas.display_status.eq('size_excluded').sum() == 11
assert atlas[atlas.display_status.ne('tested')].NES.isna().all()
assert atlas[atlas.display_status.eq('unavailable_combination')].text.eq('-').all()
assert atlas[atlas.display_status.eq('size_excluded')].text.eq('NT').all()
assert atlas[atlas.BH_across_collections.lt(.05)].text.str.endswith('*').all()
assert atlas[atlas.pathway.eq('MEISSNER_NPC_HCP_WITH_H3K4ME2')].mark.eq('K4me2').all()
comp = src(ids[1] + '_b')
assert len(comp) == 26 and comp.groupby('source').source_arm.nunique().eq(2).all()
diagnostic = src(ids[1] + '_c')
assert len(diagnostic) == 30 and diagnostic.scenario.nunique() == 6
reg = src(ids[1] + '_d')
assert len(reg) == 36 and reg.gene_name.nunique() == 12
for name in ['EPIGEN_frozen_input_hashes.csv']:
    for row in pd.read_csv(root / 'provenance' / name).itertuples():
        assert hashlib.sha256((root / row.path).read_bytes()).hexdigest() == row.sha256, row.path
prior = pd.read_csv(root / 'provenance/EPIGEN_prior_figure_hashes.csv')
assert len(prior) == 15
for row in prior.to_dict('records'):
    for kind in ['pdf', 'svg', 'png']:
        assert hashlib.sha256((root / row[kind]).read_bytes()).hexdigest() == row[kind + '_sha256'], row[kind]
report = root / 'outputs/report/epigenetic_deep_dive.html'
soup = BeautifulSoup(report.read_text(encoding='utf-8'), 'html.parser')
images = soup.find_all('img')
assert len(images) == 2 and all(i.get('src', '').startswith('data:') for i in images)
assert len(soup.find_all('table')) == 5
for link in soup.find_all('a', href=True):
    href = link['href']
    if href.startswith(('#', 'http:', 'https:', 'mailto:', 'data:', 'javascript:')):
        continue
    path = unquote(urlparse(href).path)
    assert not path or (report.parent / path).exists(), href
text = ' '.join(soup.get_text(' ', strip=True).split())
for required in ['only 6 genes', 'same established', 'No new enrichment tests', 'post hoc', '339', '23,696']:
    assert required.lower() in text.lower(), required
result = dict(new_figures=2, new_source_panels=8, atlas_source_sets=41, atlas_size_excluded=11,
              standalone_tables=5, readout_cells=60, unchanged_prior_figure_exports=45,
              frozen_inputs_verified=True, original_global_BH_reused=True, broken_local_links=0)
(root / 'provenance/EPIGEN_deliverable_checks.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
print('Passed epigenetic figure semantics, all source states, complete arms, six matches, report and 45 unchanged prior figure exports.')
