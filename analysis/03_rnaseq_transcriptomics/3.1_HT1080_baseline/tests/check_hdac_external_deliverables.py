"""Verify figure coverage, readout missingness, report links and unchanged prior exports."""
import hashlib
import json
from pathlib import Path
from urllib.parse import unquote, urlparse
import numpy as np
import pandas as pd
from bs4 import BeautifulSoup

root = Path(__file__).resolve().parents[1]
table = lambda n: pd.read_csv(root / 'outputs/tables' / (n + '.csv'))
src = lambda n: pd.read_csv(root / 'outputs/source_data' / (n + '.csv'))
cc = table('HDACEXT_contrast_manifest')
selected = set(cc.loc[cc.display_tier.eq('primary') | cc.comparison_type.eq('active_transgene_control_vs_WT'), 'contrast'])
a, b, c, d = [src('Fig09_HDAC_external_comparison_' + p) for p in 'abcd']
assert len(selected) == 19 and [len(q) for q in [a, b, c, d]] == [57, 38, 54, 48]
assert set(a.external_contrast) == selected == set(b.external_contrast)
assert a.label.notna().all() and b.label.notna().all()
assert a.label.nunique() == b.label.nunique() == 19
assert a.set_index('external_contrast').label.to_dict() == b.set_index('external_contrast').label.to_dict()
joined = a.merge(table('HDACEXT_cross_study_concordance'), on=['scenario', 'contrast', 'external_contrast'],
                 suffixes=('_figure', '_table'), validate='one_to_one')
np.testing.assert_allclose(joined.spearman_figure, joined.spearman_table, rtol=0, atol=1e-12)
arms = table('HDACEXT_HT1080_source_arm_effects')
for row in b.itertuples():
    q = arms[arms.scenario.eq(row.scenario) & arms.contrast.eq(row.contrast) & arms.external_contrast.eq(row.external_contrast)].iloc[0]
    assert abs(row.mean_effect - q[row.source_arm + '_mean']) < 1e-12
    assert row.genes == q[row.source_arm + '_genes']
assert c.pathway.nunique() == 54 and c.family_size.eq(54).all()
assert c.scenario.eq('original').all() and c.contrast.eq('combined').all()
np.testing.assert_allclose(c.fgsea_support, -np.log10(c.fgsea_secondary_BH), atol=1e-12)
np.testing.assert_allclose(c.camera_support, -np.log10(c.camera_fixed_secondary_BH), atol=1e-12)
assert (c.camera_estimated_secondary_BH < .05).sum() == 0
shortlist = table('HDAC12_readout_shortlist')
assert set(d.gene) == set(shortlist.gene_name) and d.groupby('gene').size().eq(6).all()
assert d.value.isna().sum() == 4 and d.loc[d.value.isna(), 'gene'].eq('BEX1').all()
assert d.loc[d.scope.str.startswith('HT1080'), 'value'].notna().all()
evidence = table('HDACEXT_eight_readout_comparison')
for row in d.itertuples():
    if row.scope.startswith('HT1080'):
        value = shortlist.loc[shortlist.gene_name.eq(row.gene), row.contrast].iloc[0]
    else:
        value = evidence.loc[evidence.HT1080_symbol.eq(row.gene) & evidence.external_contrast.eq(row.contrast), 'external_logFC'].iloc[0]
    np.testing.assert_allclose(row.value, value, rtol=0, atol=1e-12, equal_nan=True)
sa, sb, sc = [src('SuppFig06_HDAC_external_QC_sensitivity_' + p) for p in 'abc']
assert len(sa) == 63 and sa.sample_key.nunique() == 63 and sa.groupby('model_group').size().eq(3).all()
assert sa.annotation.notna().sum() == 21
assert len(sb) == 9 and sb.source_worksheet.nunique() == 9
assert len(sc) == 108 and sc.external_contrast.nunique() == 27 and sc.scenario.nunique() == 4
assert sc.label.nunique() == 27 and sc.label.notna().all() and sc['mask'].notna().all()
report = root / 'outputs/report/HDAC12_external_comparison.html'
soup = BeautifulSoup(report.read_text(encoding='utf-8'), 'html.parser')
broken = []
for link in soup.find_all('a', href=True):
    href = link['href']
    if href.startswith(('#', 'https:', 'http:', 'data:', 'javascript:', 'mailto:')):
        continue
    path = unquote(urlparse(href).path)
    if path and not (report.parent / path).exists():
        broken.append(href)
assert not broken, broken
images = soup.find_all('img')
assert len(images) == 2 and all(i.get('src', '').startswith('data:') for i in images)
assert len(soup.find_all('table')) == 5
text = ' '.join(soup.get_text(' ', strip=True).lower().split())
for phrase in ['does not reproduce a specific HDAC2 catalytic-loss state', '14,881', 'independent clone',
               'does not verify its treatment identity', 'Optional LINCS, complex-membership analysis and addback integration remain unrun']:
    assert phrase.lower() in text, phrase
prior = pd.read_csv(root / 'provenance/HDACEXT_prior_figure_hashes.csv')
assert len(prior) == 13
for row in prior.to_dict('records'):
    for kind in ['pdf', 'svg', 'png']:
        assert hashlib.sha256((root / row[kind]).read_bytes()).hexdigest() == row[kind + '_sha256'], row[kind]
frozen = pd.read_csv(root / 'provenance/HDACEXT_frozen_input_hashes.csv')
for row in frozen.itertuples():
    assert hashlib.sha256((root / row.path).read_bytes()).hexdigest() == row.sha256, row.path
result = dict(primary_and_active_control_contrasts=19, complete_secondary_sets=54, readout_cells=48,
    explicit_missing_external_cells=4, external_libraries=63, sensitivity_contrast_universes=108,
    standalone_figures=len(images), standalone_tables=5, broken_local_links=broken,
    unchanged_prior_exports=39, frozen_inputs_verified=len(frozen))
(root / 'provenance/HDACEXT_deliverable_checks.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
print(json.dumps(result, indent=2))
