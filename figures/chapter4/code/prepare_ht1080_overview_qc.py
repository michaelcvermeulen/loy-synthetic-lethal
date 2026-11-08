"""Verify and export source data for the Chapter 4 overview/QC drafts.

Reads original screen and Chapter 3 sources without modifying them.
No differential hit calling, FDR recalculation, or reference-set updates.
"""
from pathlib import Path
import hashlib
import json
import sys
import numpy as np
import pandas as pd

THESIS = Path(__file__).resolve().parents[3]
SCREEN = THESIS.parent / 'crispr_screen/ht1080_loy_screen_2'
OUT = THESIS / 'figures/chapter4/source_data/ht1080_figure1_qc'
SAMPLES = ['D0_CTRL', 'D0_LOY', 'POST_CTRL_REP1', 'POST_CTRL_REP2',
           'POST_LOY_REP1', 'POST_LOY_REP2']
SOURCES = {
    'counts': SCREEN / 'count/ht1080_screen2_R1.count.txt',
    'sequencing_qc': SCREEN / 'count/ht1080_screen2_R1.countsummary.txt',
    'saved_gene_lfc': SCREEN / 'results/module8_essential_depletion_per_condition.csv',
    'essential_reference': SCREEN / 'resources/AchillesCommonEssentialControls.csv',
    'nonessential_reference': SCREEN / 'resources/AchillesNonessentialControls.csv',
    'original_qc_code': SCREEN / 'code/ht1080_screen2_analysis.Rmd',
    'mageck_guides': SCREEN / 'test/POST_LOY_vs_CTRL.sgrna_summary.txt',
    'mageck_command': SCREEN / 'code/ht1080_screen2_07_mageck_test.sh',
    'ntc_ids': SCREEN / 'brunello/NTC_control_library.txt',
    'sample_metadata': SCREEN / 'sample_metadata.txt',
    'passages': SCREEN / 'results/screen_assay_passages/HT1080_Brunello_Passage_Tracking_with_4_Passages.xlsx',
    'dna': THESIS / 'figures/chapter3/source_data/ht1080_lpwgs_confirmation/regional_medians_verified.csv',
    'rna': THESIS / 'figures/chapter3/source_data/ht1080_y_expression_confirmation/relative_expression_group_means.csv',
    'y_cytobands': THESIS / 'figures/chapter3/source_data/ht1080_y_expression_confirmation/T2T_chrY_cytobands.csv',
    'y_gene_coordinates': THESIS / 'figures/chapter3/source_data/ht1080_y_expression_confirmation/gene_mapping_plotted.csv',
}


def roc_curve(labels, scores):
    """Empirical ROC with tied scores admitted together; essentials are positive."""
    frame = pd.DataFrame({'score': scores, 'positive': np.asarray(labels, dtype=int)})
    grouped = frame.groupby('score').positive.agg(['sum', 'count']).sort_index(ascending=False)
    tp = np.r_[0, grouped['sum'].cumsum().to_numpy()]
    fp = np.r_[0, (grouped['count'] - grouped['sum']).cumsum().to_numpy()]
    result = pd.DataFrame({'false_positive_rate': fp / fp[-1], 'true_positive_rate': tp / tp[-1]})
    auc = float(np.trapezoid(result.true_positive_rate, result.false_positive_rate))
    return result, auc


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    counts = pd.read_csv(SOURCES['counts'], sep='\t')
    assert len(counts) == 77441 and counts.sgRNA.is_unique
    raw = counts[SAMPLES].astype(float)
    assert raw.ge(0).all().all() and raw.notna().all().all()
    cpm = raw.div(raw.sum(), axis=1) * 1e6
    log = np.log2(cpm + 1)
    pd.concat([counts[['sgRNA', 'Gene']], log], axis=1).to_csv(OUT / 'guide_log2_cpm.csv', index=False)
    pearson = log.corr(method='pearson')
    spearman = raw.corr(method='spearman')
    pearson.to_csv(OUT / 'pearson_log2_cpm.csv')
    spearman.to_csv(OUT / 'spearman_guide_counts.csv')
    qc = pd.read_csv(SOURCES['sequencing_qc'], sep='\t').set_index('Label').loc[SAMPLES].copy()
    qc['percent_mapped'] = 100 * qc.Mapped / qc.Reads
    qc['percent_detected'] = 100 * (1 - qc.Zerocounts / qc.TotalsgRNAs)
    assert np.array_equal(qc.Zerocounts.to_numpy(), raw.eq(0).sum().to_numpy())
    qc.to_csv(OUT / 'sequencing_qc.csv')
    # Preserve MAGeCK's saved index; independently check its log-count scale.
    gini = qc[['GiniIndex']].copy()
    n_guides = len(raw)
    for sample in SAMPLES:
        values = np.sort(np.log2(raw[sample].to_numpy() + 1))
        gini.loc[sample, 'recomputed_log_count_gini'] = (
            2 * np.dot(np.arange(1, n_guides + 1), values) / (n_guides * values.sum())
            - (n_guides + 1) / n_guides)
    gini['absolute_difference'] = (gini.GiniIndex - gini.recomputed_log_count_gini).abs()
    assert gini.absolute_difference.max() < 1e-5
    gini.to_csv(OUT / 'gini_per_sample.csv')

    # Independently reproduce the saved R QC scores (R sd uses ddof=1).
    d0 = log[['D0_CTRL', 'D0_LOY']]
    z = (d0 - d0.mean()) / d0.std(ddof=1)
    excluded = z.abs().gt(3).any(axis=1)
    assert excluded.sum() == 173
    counts.loc[excluded, ['sgRNA', 'Gene']].to_csv(OUT / 'd0_outlier_guides.csv', index=False)
    filtered = counts.loc[~excluded].copy()
    fraw = filtered[SAMPLES].astype(float)
    fcpm = fraw.div(fraw.sum(), axis=1) * 1e6
    saved = pd.read_csv(SOURCES['saved_gene_lfc'])
    saved_real = saved.loc[~saved.Gene.str.startswith('RandomGene_')].set_index('Gene')
    references = {}
    for cat, key in [('Essential', 'essential_reference'), ('Non-Essential', 'nonessential_reference')]:
        references[cat] = set(pd.read_csv(SOURCES[key]).Gene.str.extract(r'^([^\s(]+)', expand=False))
        assert set(saved_real.index[saved_real.category.eq(cat)]) == set(saved_real.index) & references[cat]
    assert not references['Essential'] & references['Non-Essential']
    reproduction = {}
    roc_summary = {}
    for group in ['CTRL', 'LOY']:
        endpoint = fcpm[[f'POST_{group}_REP1', f'POST_{group}_REP2']].mean(axis=1)
        lfc = np.log2((endpoint + 1) / (fcpm[f'D0_{group}'] + 1))
        retained = filtered[f'D0_{group}'].ge(50) & filtered.Gene.ne('Non-Targeting-Control')
        reproduced = lfc[retained].groupby(filtered.loc[retained, 'Gene']).mean()
        target = saved_real[f'LFC_{group}'].dropna()
        assert set(reproduced.index) == set(target.index)
        err = float((reproduced.loc[target.index] - target).abs().max())
        assert err < 1e-10, (group, err)
        reproduction[group] = {'max_absolute_difference': err, 'guides_passing_d0_filter': int(retained.sum())}
        use = saved_real.loc[saved_real.category.isin(references) & saved_real[f'LFC_{group}'].notna()].copy()
        use['score'] = -use[f'LFC_{group}']
        use[['category', f'LFC_{group}', 'score']].to_csv(OUT / f'roc_input_{group}.csv')
        curve, auc = roc_curve(use.category.eq('Essential'), use.score)
        curve.to_csv(OUT / f'roc_curve_{group}.csv', index=False)
        roc_summary[group] = {'auc': auc,
            'n_essential': int(use.category.eq('Essential').sum()),
            'n_nonessential': int(use.category.eq('Non-Essential').sum()),
            'median_essential_lfc': float(use.loc[use.category.eq('Essential'), f'LFC_{group}'].median()),
            'median_nonessential_lfc': float(use.loc[use.category.eq('Non-Essential'), f'LFC_{group}'].median())}
    saved_real.loc[saved_real.category.isin(references)].to_csv(OUT / 'reference_gene_depletion.csv')

    guides = pd.read_csv(SOURCES['mageck_guides'], sep='\t')
    ntc_ids = set(SOURCES['ntc_ids'].read_text(encoding='utf-8').splitlines())
    ntc = guides.loc[guides.sgrna.isin(ntc_ids)].copy()
    assert len(ntc) == len(ntc_ids) == 1000 and ntc.Gene.eq('Non-Targeting-Control').all()
    # Confirm LFC direction from normalised means; saved precision is limited.
    assert ((ntc.LFC > 0) == (ntc.treat_mean > ntc.control_mean)).all()
    ntc.to_csv(OUT / 'ntc_guide_endpoint_lfc.csv', index=False)
    q025, q975 = ntc.LFC.quantile([.025, .975])
    passages = pd.read_excel(SOURCES['passages'])
    assert len(passages) == 16 and set(passages.Passage) == {'P1', 'P2', 'P3', 'P4'}
    for _, rows in passages.groupby(['Group', 'Replicate']):
        rows = rows.sort_values('Passage')
        np.testing.assert_allclose(rows['Population Doublings'], np.log2(rows['Cells Harvested (M)'] / rows['Cells Plated (M)']), atol=1e-10)
        np.testing.assert_allclose(rows['Cumulative Doublings'], rows['Population Doublings'].cumsum(), atol=1e-10)
    passages.to_csv(OUT / 'population_doublings.csv', index=False)

    dna = pd.read_csv(SOURCES['dna'])
    dna = dna.loc[dna.region.isin(['MSY Yp', 'MSY Yq (euchromatic)'])].copy()
    dna['analysis_floor'] = dna.median_relative_cn.le(.001000001)
    dna.to_csv(OUT / 'dna_inset.csv', index=False)
    genes = ['RPS4Y1', 'ZFY', 'USP9Y', 'DDX3Y', 'UTY', 'NLGN4Y', 'KDM5D', 'EIF1AY']
    rna = pd.read_csv(SOURCES['rna']).set_index('gene_name').loc[genes]
    rna.to_csv(OUT / 'rna_inset.csv')
    regressions = {}
    pairs = {'CTRL': ('POST_CTRL_REP1', 'POST_CTRL_REP2'),
             'LOY': ('POST_LOY_REP1', 'POST_LOY_REP2'), 'D0': ('D0_CTRL', 'D0_LOY')}
    for group, (x_name, y_name) in pairs.items():
        x, y = log[x_name].to_numpy(), log[y_name].to_numpy()
        slope, intercept = np.polyfit(x, y, 1)
        regressions[group] = {'n_guides': len(x), 'pearson_r': float(np.corrcoef(x, y)[0, 1]),
                              'slope': float(slope), 'intercept': float(intercept),
                              'x_sample': x_name, 'y_sample': y_name}
    summary = {'status': 'Draft figure sources verified; differential hit analyses unchanged',
        'n_guides': len(counts), 'd0_outliers_excluded_for_essential_qc': int(excluded.sum()),
        'lfc_reproduction': reproduction, 'roc': roc_summary, 'regression': regressions,
        'gini': {'source': 'Saved MAGeCK GiniIndex; plotted without substitution',
                 'scale_check': 'Gini of log2(raw guide count + 1)',
                 'max_absolute_difference_from_recomputation': float(gini.absolute_difference.max())},
        'spearman_range': [float(spearman.to_numpy().min()), float(spearman.to_numpy().max())],
        'ntc': {'n': len(ntc), 'mean': float(ntc.LFC.mean()), 'median': float(ntc.LFC.median()),
                'sd': float(ntc.LFC.std(ddof=1)), 'q025': float(q025), 'q975': float(q975)},
        'population_doublings_check': 'Per-passage log2(harvested/plated) and cumulative sums verified',
        'dna_rna': 'Reused Chapter 3 population summaries with separate assay scales',
        'reference_set_provenance': 'Local Achilles control lists used by the original R analysis; release not inferred'}
    (OUT / 'summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    manifest = {k: {'path': str(v), 'sha256': hashlib.sha256(v.read_bytes()).hexdigest()} for k, v in SOURCES.items()}
    manifest['preparation_script'] = {'path': str(Path(__file__).resolve()), 'sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (OUT / 'source_manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
