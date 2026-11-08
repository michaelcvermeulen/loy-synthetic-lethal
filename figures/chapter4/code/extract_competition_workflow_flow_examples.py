"""Extract two prespecified WT D0 wells using saved gates and classifier.

No classifier fitting or source mutation. First replicate of each orientation.
"""
from pathlib import Path
import csv
import hashlib
import json
import sys

THESIS = Path(__file__).resolve().parents[3]
PROJECTS = THESIS.parent
sys.path.insert(0, str(PROJECTS / '.figure_python_deps'))
sys.path.insert(0, str(PROJECTS / 'comp_assay'))
import numpy as np
from spectral.a8_pipeline.competition_report import _primary_features, load_cef_gates

REPO = PROJECTS / 'comp_assay'
SOURCE = REPO / 'analysis/spectral_growth_assay_D0_20261015'
EXPORT = REPO / '20261016_171604_Export/comp_growth_assay_D0_OCT_15'
OUT = THESIS / 'figures/chapter4/source_data/competition_workflow'
SAMPLES = ['A01-LOY_CHERRY_WT_1', 'A04-LOY_GFP_WT_1_001']
result = json.loads((SOURCE / 'results.json').read_text(encoding='utf-8'))
rows = {row['sample']: row for row in result['wells']}
cef = next(EXPORT.glob('*.cef'))
gates = load_cef_gates(cef)
weights = np.asarray(result['classifier']['weights'])
intercept = result['classifier']['intercept']
gfp_limit = result['negative_rules']['broad']['gfp_below']
mcherry_limit = result['negative_rules']['broad']['mcherry_below']

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

manifest = {
    'acquisition_date': '2026-10-15',
    'selection': 'First replicate of each standard WT dye orientation (A01 and A04); no outcome-based selection.',
    'source_report': str(SOURCE / 'report.md'),
    'classifier': result['classifier'],
    'negative_rule': {'name': 'broad', 'limits': [gfp_limit, mcherry_limit]},
    'transform': 'arcsinh(unmixed intensity / 1000)',
    'gate_hierarchy': ['Scatter', 'SSC Singlets', 'FSC Singlets'],
    'sources': [{'path': str(p), 'sha256': digest(p)} for p in (cef, SOURCE / 'results.json', SOURCE / 'wells.csv')],
    'wells': [],
}
for mixture, name in enumerate(SAMPLES, 1):
    path = EXPORT / 'Plate_004' / (name + '.fcs')
    before = digest(path)
    values, count, _ = _primary_features(path, gates)
    negative = (values[:, 0] < gfp_limit) & (values[:, 1] < mcherry_limit)
    gfp = values @ weights + intercept > 0
    row = rows[name]
    assert count == row['primary_singlets']
    assert int((~negative).sum()) == row['reporter_positive_singlets']
    assert int((gfp & ~negative).sum()) == row['gfp_count']
    loy = gfp if mixture == 2 else ~gfp
    assert int((loy & ~negative).sum()) == row['loy_count']
    assert digest(path) == before
    csv_path = OUT / f'flow_D0_WT_mixture_{mixture}_events.csv'
    with csv_path.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.writer(handle)
        writer.writerow(['singlet_index', 'eGFP_arcsinh', 'mCherry_arcsinh', 'reporter_class', 'excluded'])
        for idx, (point, cls, excluded) in enumerate(zip(values, gfp, negative)):
            writer.writerow([idx, *point, 'GFP' if cls else 'mCherry', int(excluded)])
    entry = {key: row[key] for key in ('sample', 'well', 'orientation', 'replicate', 'primary_singlets', 'reporter_positive_singlets', 'gfp_count', 'gfp_fraction', 'loy_count', 'loy_fraction', 'negative_broad_fraction')}
    entry.update({'mixture': mixture, 'events_csv': csv_path.name, 'source_fcs': str(path), 'source_sha256': before, 'data_min': values.min(axis=0).tolist(), 'data_max': values.max(axis=0).tolist()})
    manifest['wells'].append(entry)
(OUT / 'flow_examples_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n', encoding='utf-8')
print(json.dumps(manifest['wells'], indent=2))
