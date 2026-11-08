"""Export prespecified reporter controls and a historical unstained example.

Keep each acquisition's own CEF singlet gates. The September unstained sample
is descriptive only: do not apply October exclusion or reporter classification.
"""
from pathlib import Path
import csv
import hashlib
import json
import sys

THESIS = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(THESIS.parent / '.figure_python_deps'))
sys.path.insert(0, str(THESIS.parent / 'comp_assay'))
import numpy as np
from spectral.a8_pipeline.competition_report import _primary_features, load_cef_gates
from spectral.a8_pipeline.fcs_io import load_fcs

REPO = THESIS.parent / 'comp_assay'
OUT = THESIS / 'figures/chapter4/source_data/competition_workflow'
OCT = REPO / '20261016_171604_Export/comp_growth_assay_D0_OCT_15'
SEP = REPO / 'spectral/20260901_125230_Mike_Run1/Mike_Run1'
RESULTS = REPO / 'analysis/spectral_growth_assay_D0_20261015/results.json'
saved = json.loads(RESULTS.read_text(encoding='utf-8'))
rows = {row['sample']: row for row in saved['wells']}
weights = np.asarray(saved['classifier']['weights'])
intercept = saved['classifier']['intercept']
limits = [saved['negative_rules']['broad']['gfp_below'],
          saved['negative_rules']['broad']['mcherry_below']]

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

specifications = [
    ('GFP only', 'A10-CTRL_GFP_SOLO_WT_003', OCT, 'Plate_004', '2026-10-15'),
    ('mCherry only', 'A09-CTRL_CHERRY_SOLO_WT_002', OCT, 'Plate_004', '2026-10-15'),
    ('Unstained', 'A01-unstained_CTRL', SEP, '1', '2026-09-01'),
]
manifest = {
    'selection': 'CTRL WT single-reporter wells and CTRL unstained, selected by sample role before viewing; no outcome-based selection.',
    'transform': 'arcsinh(unmixed intensity / 1000)',
    'gate_hierarchy': ['Scatter', 'SSC Singlets', 'FSC Singlets'],
    'gate_source': 'Each sample uses its own acquisition CEF.',
    'historical_control_limit': 'No unstained sample in the October plate. September unmixed background is descriptive, not a transferable October gate or a same-day QC control.',
    'source_report': str(RESULTS), 'source_report_sha256': digest(RESULTS),
    'controls': [],
}
for label, name, root, plate, date in specifications:
    path = root / plate / (name + '.fcs')
    cef = next(root.glob('*.cef'))
    before = digest(path)
    values, count, _ = _primary_features(path, load_cef_gates(cef))
    data = load_fcs(path)
    historical = label == 'Unstained'
    if not historical:
        excluded = (values[:, 0] < limits[0]) & (values[:, 1] < limits[1])
        gfp = values @ weights + intercept > 0
        row = rows[name]
        assert count == row['primary_singlets']
        assert int((~excluded).sum()) == row['reporter_positive_singlets']
        assert int((gfp & ~excluded).sum()) == row['gfp_count']
    assert digest(path) == before
    csv_path = OUT / ('flow_control_' + name + '_events.csv')
    with csv_path.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.writer(handle)
        writer.writerow(['singlet_index', 'eGFP_arcsinh', 'mCherry_arcsinh', 'reporter_class', 'excluded'])
        for idx, point in enumerate(values):
            writer.writerow([idx, *point,
                             'unclassified' if historical else ('GFP' if gfp[idx] else 'mCherry'),
                             '' if historical else int(excluded[idx])])
    record = {
        'label': label, 'sample': name, 'genotype': 'CTRL', 'experiment_date': date,
        'acquisition_date': '2026-09-01' if historical else '2026-10-16',
        'historical_unstained': historical, 'primary_singlets': count,
        'events_csv': csv_path.name, 'source_fcs': str(path), 'source_sha256': before,
        'source_cef': str(cef), 'source_cef_sha256': digest(cef),
        'unmixing_matrix': data.matrix_name,
        'fcs_acquisition_date': data.metadata.get('date'),
        'date_note': 'October experiment is labelled 15 October; its FCS acquisition date is 16 October.' if not historical else '',
        'data_min': values.min(axis=0).tolist(), 'data_max': values.max(axis=0).tolist(),
        'october_classifier_applied': not historical, 'october_negative_rule_applied': not historical,
    }
    if not historical:
        record.update({key: row[key] for key in ('reporter_positive_singlets', 'gfp_count', 'gfp_fraction', 'negative_broad_fraction')})
    manifest['controls'].append(record)
(OUT / 'flow_controls_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n', encoding='utf-8')
print(json.dumps(manifest['controls'], indent=2))
