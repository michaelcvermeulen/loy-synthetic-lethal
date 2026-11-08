"""A549 D0/P4 QC counterparts to the HT1080 panels; original inputs are read-only."""
from pathlib import Path
import hashlib
import json
import shutil
import sys
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT.parent / '.figure_python_deps'))
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde
from prepare_ht1080_overview_qc import roc_curve

INITIAL = ROOT / 'figures/chapter4/source_data/a549_initial_qc'
OUT = ROOT / 'figures/chapter4/source_data/a549_matched_qc'
SCREEN = ROOT.parent / 'crispr_screen/a549_loy_screen1/CRISPR_screen'

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    original = json.loads((INITIAL / 'analysis_manifest.json').read_text())
    sources = dict(original['inputs'])
    for key, path in {
        'retained_counts': INITIAL / 'retained_counts.tsv',
        'baseline_gene_qc': INITIAL / 'baseline_gene_qc.csv',
        'sample_qc': INITIAL / 'sample_qc.csv',
        'sample_metadata': INITIAL / 'sample_metadata.csv',
        'mageck_guides': SCREEN / 'Screen1/Mageck/P4_LOY_vs_CTRL.sgrna_summary.txt',
    }.items():
        sources[key] = {'path': str(path), 'sha256': sha(path)}
    for key, item in sources.items():
        assert sha(item['path']) == item['sha256'], key
    samples = original['retained_samples']
    assert len(samples) == 7 and all('_P3_' not in x for x in samples)
    counts = pd.read_csv(sources['retained_counts']['path'], sep='\t').set_index('sgRNA')
    assert len(counts) == 77441 and counts.index.is_unique
    raw = counts[samples].astype(float)
    log = np.log2(raw.div(raw.sum(), axis=1) * 1e6 + 1)
    log.to_csv(OUT / 'guide_log2_cpm.csv')
    raw.corr(method='spearman').to_csv(OUT / 'spearman_guide_counts.csv')
    for name in ['sample_metadata.csv', 'sample_qc.csv', 'baseline_gene_qc.csv']:
        shutil.copy2(INITIAL / name, OUT / name)

    # Match the HT1080 scaled PCA (R prcomp scale.=TRUE); old unscaled QC is preserved.
    sd = log.std(axis=1, ddof=1)
    use = sd.gt(0)
    z = log.loc[use].sub(log.loc[use].mean(axis=1), axis=0).div(sd.loc[use], axis=0)
    u, singular, vt = np.linalg.svd(z.to_numpy().T, full_matrices=False)
    scores = u * singular
    for component in range(3):
        anchor = np.argmax(np.abs(vt[component]))
        if vt[component, anchor] < 0:
            scores[:, component] *= -1
    explained = singular ** 2 / np.sum(singular ** 2)
    eigen = np.linalg.eigvalsh(z.to_numpy().T @ z.to_numpy())[::-1]
    np.testing.assert_allclose(explained[:3], eigen[:3] / eigen.sum(), atol=1e-10)
    meta = pd.read_csv(OUT / 'sample_metadata.csv')
    meta[['PC1', 'PC2', 'PC3']] = scores[:, :3]
    meta.to_csv(OUT / 'pca_sample_scores.csv', index=False)
    pd.DataFrame({'component': np.arange(1, 8), 'variance_percent': 100 * explained}).to_csv(OUT / 'pca_variance.csv', index=False)
    gini = pd.read_csv(OUT / 'sample_qc.csv').set_index('sample')[['gini_mageck_export']]
    for sample in samples:
        vals = np.sort(np.log2(raw[sample].to_numpy() + 1))
        n = len(vals)
        gini.loc[sample, 'recomputed_log_count_gini'] = 2 * np.dot(np.arange(1, n + 1), vals) / (n * vals.sum()) - (n + 1) / n
    gini['absolute_difference'] = abs(gini.gini_mageck_export - gini.recomputed_log_count_gini)
    # The original report prints four significant digits (P4 values have 4 decimals).
    assert gini.absolute_difference.max() < 5e-5
    gini.to_csv(OUT / 'gini_per_sample.csv')

    pooled = pd.DataFrame({g: raw[[s for s in samples if s.startswith(g + '_D0')]].sum(axis=1) for g in ['CTRL', 'LOY']})
    pooled_log = np.log2(pooled.div(pooled.sum(), axis=1) * 1e6 + 1)
    pooled_log.to_csv(OUT / 'pooled_d0_log2_cpm.csv')
    pairs = {'D0_pooled': (pooled_log.CTRL, pooled_log.LOY),
             'CTRL_D0_technical': (log.CTRL_D0_1, log.CTRL_D0_2),
             'LOY_D0_technical': (log.LOY_D0_1, log.LOY_D0_2),
             'LOY_P4': (log.LOY_P4_R1, log.LOY_P4_R2)}
    regressions = {}
    for key, (x, y) in pairs.items():
        slope, intercept = np.polyfit(x, y, 1)
        regressions[key] = {'n_guides': len(x), 'pearson_r': float(np.corrcoef(x, y)[0, 1]),
                            'slope': float(slope), 'intercept': float(intercept),
                            'maximum_abundance': float(max(x.max(), y.max()))}
    genes = pd.read_csv(OUT / 'baseline_gene_qc.csv')
    roc = {}
    old_summary = json.loads((INITIAL / 'summary.json').read_text())
    for group in ['CTRL', 'LOY']:
        d = genes.loc[genes.genotype.eq(group) & genes.category.isin(['Essential', 'Non-essential'])]
        curve, auc = roc_curve(d.category.eq('Essential'), -d.lfc)
        assert abs(auc - next(x['roc_auc'] for x in old_summary['biological_qc'] if x['genotype'] == group)) < 1e-12
        curve.to_csv(OUT / f'roc_curve_{group}.csv', index=False)
        roc[group] = {'auc': auc, 'n_essential': int(d.category.eq('Essential').sum()), 'n_nonessential': int(d.category.eq('Non-essential').sum())}

    # Use saved NTC-normalised endpoint LFCs, exactly as for the HT1080 density panel.
    guides = pd.read_csv(sources['mageck_guides']['path'], sep='\t').set_index('sgrna')
    assert guides.index.is_unique and set(guides.index) == set(counts.index)
    assert (guides.Gene == counts.Gene.reindex(guides.index)).all()
    refs = {cat: set(pd.read_csv(sources[key]['path']).iloc[:, 0].str.replace(r'\s+\(.*$', '', regex=True))
            for cat, key in [('Essential', 'essential_reference'), ('Non-essential', 'nonessential_reference')]}
    assert not refs['Essential'] & refs['Non-essential']
    ntc_ids = set(Path(sources['ntc_reference']['path']).read_text().splitlines())
    assert len(ntc_ids) == 1000
    masks = {'Non-targeting': guides.index.isin(ntc_ids), **{c: guides.Gene.isin(r) for c, r in refs.items()}}
    endpoint_total = raw[[s for s in samples if '_P4_' in s]].sum(axis=1)
    parts = []
    for cat, mask in masks.items():
        d = guides.loc[mask, ['Gene', 'LFC']].copy()
        d['category'] = cat
        d['endpoint_raw_reads'] = endpoint_total.reindex(d.index)
        d['included'] = d.endpoint_raw_reads.gt(0)
        parts.append(d)
    table = pd.concat(parts).reset_index()
    assert table.sgrna.is_unique and np.isfinite(table.LFC).all()
    table.to_csv(OUT / 'endpoint_reference_density_guides.csv', index=False)
    eligible = table.loc[table.included]
    ntc = eligible.loc[eligible.category.eq('Non-targeting'), 'LFC'].to_numpy()
    bandwidth = float(np.std(ntc, ddof=1) * len(ntc) ** (-1 / 5))
    grid = np.linspace(eligible.LFC.min() - 4 * bandwidth, eligible.LFC.max() + 4 * bandwidth, 2401)
    curves, groups = {'LFC': grid}, {}
    for cat in masks:
        values = eligible.loc[eligible.category.eq(cat), 'LFC'].to_numpy()
        curves[cat] = gaussian_kde(values, bw_method=bandwidth / np.std(values, ddof=1))(grid)
        integral = float(np.trapezoid(curves[cat], grid))
        assert abs(integral - 1) < .0001
        groups[cat] = {'included_guides': len(values), 'excluded_zero_endpoint_guides': int((table.category.eq(cat) & ~table.included).sum()),
                       'median': float(np.median(values)), 'sd': float(np.std(values, ddof=1)),
                       'q025': float(np.quantile(values, .025)), 'q975': float(np.quantile(values, .975)), 'integral': integral}
    pd.DataFrame(curves).to_csv(OUT / 'endpoint_reference_density_curves.csv', index=False)
    summary = {'samples': samples, 'excluded_samples': original['excluded_samples'], 'exclusion_reason': original['exclusion_reason'],
               'n_guides': len(raw), 'regression': regressions, 'roc': roc,
               'pca': {'scaled': True, 'ddof': 1, 'n_variable_guides': int(use.sum()), 'n_constant_guides_excluded': int((~use).sum()), 'variance_percent': (100 * explained[:3]).tolist()},
               'density': {'bandwidth': bandwidth, 'histogram_bin_width': .1, 'groups': groups,
                           'contrast': 'Saved NTC-normalised MAGeCK LOY P4 / CTRL P4 guide LFC', 'eligibility': 'At least one raw read across the three P4 samples; no D0 filter'},
               'biological_qc': 'Unchanged A549 scores: pooled D0 technical counts, mean P4 CPM, CPM pseudocount 1, >=50 reads in each D0 library; no D0 outlier filter.',
               'unavailable_counterparts': ['CTRL endpoint replicate comparison: one retained CTRL P4 sample', 'No verified A549 screen-passage doubling or ddPCR series supplied to this figure']}
    (OUT / 'summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    for key, item in sources.items():
        assert sha(item['path']) == item['sha256'], key
    manifest = {'inputs': sources, 'source_data_unchanged': True, 'script_sha256': sha(__file__),
                'ht1080_template': str(ROOT / 'figures/chapter4/code/make_ht1080_overview_qc.py'),
                'adaptation': 'Equivalent plot types, labels adapted to technical D0 libraries and one CTRL endpoint. Scaled PCA newly computed; baseline depletion preserved.'}
    (OUT / 'analysis_manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    print(json.dumps(summary, indent=2))

if __name__ == '__main__':
    main()
