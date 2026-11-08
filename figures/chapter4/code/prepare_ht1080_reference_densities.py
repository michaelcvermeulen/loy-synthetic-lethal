"""Prepare a descriptive guide-LFC density overlay from saved MAGeCK estimates."""
from pathlib import Path
import hashlib
import json
import sys

THESIS = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(THESIS.parent / '.figure_python_deps'))
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

DATA = THESIS / 'figures/chapter4/source_data/ht1080_figure1_qc'


def main():
    sources = json.loads((DATA / 'source_manifest.json').read_text())
    keys = ['counts', 'mageck_guides', 'essential_reference', 'nonessential_reference', 'ntc_ids']
    for key in keys:
        entry = sources[key]
        assert hashlib.sha256(Path(entry['path']).read_bytes()).hexdigest() == entry['sha256'], key
    guides = pd.read_csv(sources['mageck_guides']['path'], sep='\t').set_index('sgrna')
    counts = pd.read_csv(sources['counts']['path'], sep='\t').set_index('sgRNA')
    endpoints = ['POST_CTRL_REP1', 'POST_CTRL_REP2', 'POST_LOY_REP1', 'POST_LOY_REP2']
    total = counts[endpoints].sum(axis=1)
    ntc_ids = set(Path(sources['ntc_ids']['path']).read_text().splitlines())
    categories = {'Non-targeting': guides.index.isin(ntc_ids)}
    refs = {}
    for category, key in [('Essential', 'essential_reference'), ('Non-essential', 'nonessential_reference')]:
        refs[category] = set(pd.read_csv(sources[key]['path']).Gene.str.extract(r'^([^\s(]+)', expand=False))
        categories[category] = guides.Gene.isin(refs[category])
    assert not refs['Essential'] & refs['Non-essential']
    parts = []
    for category, mask in categories.items():
        part = guides.loc[mask, ['Gene', 'LFC']].copy()
        part['category'] = category
        part['endpoint_raw_reads'] = total.reindex(part.index)
        assert part.endpoint_raw_reads.notna().all() and np.isfinite(part.LFC).all()
        part['included'] = part.endpoint_raw_reads.gt(0)
        parts.append(part)
    table = pd.concat(parts).reset_index()
    assert table.sgrna.is_unique and len(table) == 8820
    table.to_csv(DATA / 'endpoint_reference_density_guides.csv', index=False)
    eligible = table.loc[table.included]
    ntc = eligible.loc[eligible.category.eq('Non-targeting'), 'LFC'].to_numpy()
    # Use one absolute bandwidth for all categories so smoothing is comparable.
    bandwidth = float(np.std(ntc, ddof=1) * len(ntc) ** (-1 / 5))
    lo, hi = float(eligible.LFC.min()), float(eligible.LFC.max())
    grid = np.linspace(lo - 4 * bandwidth, hi + 4 * bandwidth, 2401)
    curves = {'LFC': grid}
    groups = {}
    for category in categories:
        values = eligible.loc[eligible.category.eq(category), 'LFC'].to_numpy()
        density = gaussian_kde(values, bw_method=bandwidth / np.std(values, ddof=1))(grid)
        curves[category] = density
        integral = float(np.trapezoid(density, grid))
        assert abs(integral - 1) < 0.0001
        q025, q975 = np.quantile(values, [.025, .975])
        all_group = table.loc[table.category.eq(category)]
        groups[category] = {'total_guides': len(all_group), 'included_guides': len(values),
                            'excluded_zero_endpoint_guides': int((~all_group.included).sum()),
                            'genes': int(eligible.loc[eligible.category.eq(category), 'Gene'].nunique()),
                            'median': float(np.median(values)), 'sd': float(np.std(values, ddof=1)),
                            'q025': float(q025), 'q975': float(q975),
                            'minimum': float(values.min()), 'maximum': float(values.max()),
                            'density_integral': integral}
    pd.DataFrame(curves).to_csv(DATA / 'endpoint_reference_density_curves.csv', index=False)
    record = {'contrast': 'Saved NTC-normalised MAGeCK guide LFC: LOY endpoints / CTRL endpoints',
              'unit': 'sgRNA, not gene or independent culture',
              'eligibility': 'At least one raw read across the four endpoints; no D0 filter.',
              'bandwidth_log2_units': bandwidth, 'bandwidth_rule': 'NTC sample SD times n^(-1/5), shared across categories',
              'kernel': 'Gaussian', 'density_normalisation': 'Each curve integrates to one independently',
              'histogram_bin_width_log2_units': 0.1, 'groups': groups,
              'sources': {key: sources[key] for key in keys},
              'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (DATA / 'endpoint_reference_density_manifest.json').write_text(json.dumps(record, indent=2), encoding='utf-8')
    print(json.dumps({'bandwidth': bandwidth, 'groups': groups}, indent=2))


if __name__ == '__main__':
    main()
