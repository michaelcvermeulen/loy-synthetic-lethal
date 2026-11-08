"""Assemble complete gene-level evidence and independently check custom RNA enrichment."""
from pathlib import Path
import argparse
import hashlib
import json
import itertools
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'outputs/tables/chromatin_rna'
ATLAS = ROOT / 'outputs/tables/chromatin_atlas'
PROV = ROOT / 'provenance/chromatin_rna'


def bh(p, n=None, by=False):
    p = np.asarray(p, float)
    n = len(p) if n is None else n
    good = np.flatnonzero(np.isfinite(p))
    order = good[np.argsort(p[good], kind='stable')]
    result = np.full(len(p), np.nan)
    values = p[order] * n / np.arange(1, len(order) + 1)
    if by:
        values *= np.sum(1 / np.arange(1, n + 1))
    result[order] = np.minimum(1, np.minimum.accumulate(values[::-1])[::-1])
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--primary-only', action='store_true')
    args = parser.parse_args()
    cfg = json.loads((ROOT / 'config/chromatin_rna.json').read_text())
    if not args.primary_only:
        complete = json.loads((PROV / 'statistics_completion.json').read_text())
        assert complete['status'] == 'statistics_complete'
    manifest = pd.read_csv(PROV / 'input_manifest.csv')
    for row in manifest.itertuples():
        assert hashlib.sha256((ROOT / row.path).read_bytes()).hexdigest() == row.sha256, row.path
    if args.primary_only:
        fg = pd.concat([pd.read_csv(p) for p in (OUT/'partitions').glob('part_*/fgsea_results.csv')], ignore_index=True)
        camera = pd.concat([pd.read_csv(p) for p in (OUT/'partitions').glob('part_*/camera_results.csv')], ignore_index=True)
        fg = fg[fg.scenario.eq('window_2000')].copy()
        camera = camera[camera.scenario.eq('window_2000')].copy()
        assert len(fg) == 60 and len(camera) == 120
        for frame, p in [(fg, 'pval'), (camera, 'PValue')]:
            frame['analysis_role'] = np.where(frame.contrast.eq('combined'), 'primary', 'sensitivity')
            frame['BH_pooled_sensitivities'] = np.nan
            if frame is camera:
                for _, d in frame[frame.analysis_role.eq('sensitivity')].groupby('correlation_mode'):
                    frame.loc[d.index, 'BH_pooled_sensitivities'] = bh(d[p])
            else:
                d = frame[frame.analysis_role.eq('sensitivity')]
                frame.loc[d.index, 'BH_pooled_sensitivities'] = bh(d[p])
    else:
        fg = pd.read_csv(OUT / 'fgsea_results.csv')
        camera = pd.read_csv(OUT / 'camera_results.csv')
    # These snapshots supply complete primary/period figures without waiting on sensitivity jobs.
    # Pooled sensitivity q values are deliberately omitted from the primary-only snapshots.
    for frame, name in [(fg, 'fgsea'), (camera, 'camera')]:
        frame[frame.scenario.eq('window_2000')].drop(columns='BH_pooled_sensitivities').to_csv(
            OUT / f'{name}_primary_complete.csv', index=False)
    effects = pd.read_csv(OUT / 'gene_effect_evidence.csv').set_index('gene_name', drop=False)
    ranks = pd.read_csv(ROOT / 'outputs/tables/pathway_TF_effect_ranking.csv')
    members = pd.read_csv(ATLAS / 'gene_set_membership_RNA_universe.csv')
    catalogue = pd.read_csv(ATLAS / 'gene_set_catalogue.csv')
    sets = {n: set(members.loc[members.set_name.eq(n), 'gene_name']) for n in catalogue.set_name}
    reference = pd.read_csv(OUT / 'reference_features.csv').set_index('gene_name')
    assert len(effects) == 25461 and len(sets) == 20
    assert all(len(sets[r.set_name]) == r.n_RNA_universe for r in catalogue.itertuples())
    assert np.allclose(ranks.combined, effects.loc[ranks.gene_name, 'mean_difference'], atol=1e-12)
    for frame, p, group in [(fg, 'pval', ['scenario', 'contrast']),
                             (camera, 'PValue', ['scenario', 'contrast', 'correlation_mode'])]:
        for key, part in frame.groupby(group):
            assert len(part) == 20, key
            assert np.allclose(bh(part[p], 20), part.BH_custom_20, equal_nan=True), key
            assert np.allclose(bh(part[p], 20, True), part.BY_custom_20, equal_nan=True), key
        sens = frame[frame.analysis_role.eq('sensitivity')]
        groups = sens.groupby('correlation_mode') if p == 'PValue' else [('fgsea', sens)]
        for _, part in groups:
            assert np.allclose(bh(part[p]), part.BH_pooled_sensitivities, equal_nan=True)
    assert (camera.residual_df == 4).all()
    assert (camera.effective_correlation.dropna() >= 0).all()
    assert (fg.loc[fg.scenario.eq('evaluable_promoters_only'), 'universe_size'] == 24994).all()
    assert set(fg[fg.analysis_role.eq('primary')].pathway) == set(sets)
    bad = fg.test_status.eq('numerical_failure')
    assert not (bad & fg.analysis_role.eq('primary')).any()
    assert fg.loc[bad, 'pval'].isna().all()
    assert (fg.loc[bad, 'nPermSimple_used'] == cfg['fgsea_retry_nPermSimple']).all()
    assert not (camera.test_status == 'numerical_failure').any()

    summaries, passage_rows, influence, curve_rows, leading_rows = [], [], [], [], []
    passage_columns = [f'{period}__passage_{i}' for period in ['AUG_2025', 'MAY_2026'] for i in range(1, 4)]
    sorted_ranks = ranks.sort_values(['combined', 'gene_name'], ascending=[False, True]).reset_index(drop=True)
    values = sorted_ranks.combined.to_numpy()
    for name, genes in sets.items():
        d = effects.loc[sorted(genes)]
        rest = effects.loc[~effects.index.isin(genes)]
        row = {'pathway': name, 'size': len(d)}
        for contrast, column in [('combined', 'mean_difference'), ('AUG_2025', 'AUG_2025'), ('MAY_2026', 'MAY_2026')]:
            v = d[column]
            row.update({f'{contrast}_mean': v.mean(), f'{contrast}_median': v.median(),
                        f'{contrast}_q25': v.quantile(.25), f'{contrast}_q75': v.quantile(.75),
                        f'{contrast}_fraction_positive': (v > 0).mean(),
                        f'{contrast}_complement_mean': rest[column].mean()})
        row['fraction_genes_period_direction_agrees'] = d.period_direction_agrees.mean()
        row['fraction_genes_all_six_same_sign'] = d.matches_agreeing_with_mean.eq(6).mean()
        row['median_CPM'] = d.mean_CPM.median()
        row['median_feature_length_bp'] = d.feature_length_bp.median()
        means = d[passage_columns].mean()
        row['set_mean_same_sign_all_six'] = bool((np.sign(means) == np.sign(d.mean_difference.mean())).all())
        for column in passage_columns:
            passage_rows.append({'pathway': name, 'passage': column,
                                 'set_mean_difference': means[column],
                                 'complement_mean_difference': rest[column].mean(), 'n_genes': len(d)})
        main = fg[fg.pathway.eq(name) & fg.analysis_role.eq('primary')].iloc[0]
        hit = sorted_ranks.gene_name.isin(genes).to_numpy()
        hit_weights = np.abs(values[hit])
        inc = np.full(len(values), -1 / (~hit).sum())
        inc[hit] = hit_weights / hit_weights.sum()
        curve = np.cumsum(inc)
        es = curve.max() if curve.max() > -curve.min() else curve.min()
        assert abs(es - main.ES) < 1e-10, (name, es, main.ES)
        assert abs(curve[-1]) < 1e-10
        peak = np.argmax(curve) if es > 0 else np.argmin(curve)
        lead_mask = hit & ((np.arange(len(hit)) <= peak) if es > 0 else (np.arange(len(hit)) > peak))
        lead = set(sorted_ranks.loc[lead_mask, 'gene_name'])
        assert lead == set(str(main.leadingEdge).split(';')), name
        row['leading_edge_n'] = len(lead)
        row['leading_edge_fraction'] = len(lead) / len(genes)
        w = d.mean_difference.abs().sort_values(ascending=False)
        row['largest_gene_fraction_total_hit_weight'] = w.iloc[0] / w.sum()
        row['top5_genes_fraction_total_hit_weight'] = w.iloc[:5].sum() / w.sum()
        row['top5_absolute_effect_genes'] = ';'.join(w.index[:5])
        row['leading_edge_top_arm'] = effects.loc[sorted(lead), 'arm'].value_counts().index[0]
        row['leading_edge_top_arm_fraction'] = effects.loc[sorted(lead), 'arm'].value_counts().iloc[0] / len(lead)
        # Influence diagnostics change neither gene sets nor the reported testing family.
        for remove in [1, 5]:
            removed = set(w.index[:remove])
            eligible = ~sorted_ranks.gene_name.isin(removed).to_numpy()
            h, v = hit[eligible], values[eligible]
            increments = np.full(len(v), -1 / (~h).sum())
            increments[h] = np.abs(v[h]) / np.abs(v[h]).sum()
            cs = increments.cumsum()
            reduced_es = cs.max() if cs.max() > -cs.min() else cs.min()
            influence.append({'pathway': name, 'removed_n': remove, 'removed_genes': ';'.join(sorted(removed)),
                              'original_ES': es, 'reduced_ES': reduced_es,
                              'analysis_role': 'post_hoc_influence_diagnostic_no_P_value',
                              'same_ES_direction': np.sign(es) == np.sign(reduced_es)})
        for gene in sorted(lead):
            leading_rows.append({'pathway': name, 'gene_name': gene})
        if name in cfg['main_curve_panels']:
            curve_rows.append(pd.DataFrame({'pathway': name, 'rank': np.arange(1, len(values) + 1),
                                          'gene_name': sorted_ranks.gene_name, 'rank_metric': values,
                                          'running_ES': curve, 'is_hit': hit, 'is_leading_edge': lead_mask}))
        summaries.append(row)
    summary = pd.DataFrame(summaries)
    main_fg = fg[fg.analysis_role.eq('primary')].drop(columns='size')
    summary = summary.merge(main_fg, on='pathway', validate='1:1')
    for mode, prefix in [('fixed_0.01', 'camera_fixed'), ('estimated_nonnegative', 'camera_estimated')]:
        z = camera[camera.analysis_role.eq('primary') & camera.correlation_mode.eq(mode)]
        columns = ['Direction', 'PValue', 'BH_custom_20', 'BY_custom_20', 'effective_correlation']
        z = z[['pathway'] + columns].rename(columns={c: f'{prefix}_{c}' for c in columns})
        summary = summary.merge(z, on='pathway', validate='1:1')
    summary.to_csv(OUT / 'primary_results_with_effects.csv', index=False)
    pd.DataFrame(passage_rows).to_csv(OUT / 'set_passage_effects.csv', index=False)
    pd.DataFrame(influence).to_csv(OUT / 'top_gene_ES_influence.csv', index=False)
    pd.concat(curve_rows).to_csv(OUT / 'selected_running_curves.csv.gz', index=False, compression={'method':'gzip', 'mtime':0})
    lead_table = pd.DataFrame(leading_rows).merge(effects.reset_index(drop=True), on='gene_name', validate='m:1')
    lead_table.to_csv(OUT / 'leading_edge_gene_evidence.csv', index=False)

    overlaps = []
    lead_sets = {n: set(lead_table.loc[lead_table.pathway.eq(n), 'gene_name']) for n in sets}
    for a, b in itertools.combinations(sets, 2):
        overlaps.append({'set_a': a, 'set_b': b, 'shared_genes': len(sets[a] & sets[b]),
                         'membership_jaccard': len(sets[a] & sets[b]) / len(sets[a] | sets[b]),
                         'shared_leading_edge': len(lead_sets[a] & lead_sets[b]),
                         'leading_edge_jaccard': len(lead_sets[a] & lead_sets[b]) / len(lead_sets[a] | lead_sets[b])})
    pd.DataFrame(overlaps).to_csv(OUT / 'set_and_leading_edge_overlap.csv', index=False)

    old = pd.read_csv(ROOT / 'outputs/tables/fgsea_all_collections.csv')
    selected = old[old.contrast.eq('combined') & old.pathway.str.match(
        r'^(MIKKELSEN_(MEF|MCV6)_HCP_WITH_H3K27ME3$|PRC2_EED_DN.V1_|SENESE_HDAC)')]
    legacy_rows, legacy_genes = [], []
    for r in selected.itertuples():
        lead = set(r.leadingEdge.split(';'))
        for name, genes in sets.items():
            shared = sorted(lead & genes)
            legacy_rows.append({'source_pathway': r.pathway, 'source_NES': r.NES,
                                'source_original_global_BH': r.BH_across_collections,
                                'custom_set': name, 'source_leading_edge_n': len(lead),
                                'source_LE_unevaluable_promoters': sum(reference.loc[sorted(lead), 'reference_status'].ne('evaluable_under_frozen_rules')),
                                'shared_n': len(shared), 'fraction_source_LE': len(shared) / len(lead),
                                'custom_set_fraction_primary_background': len(genes) / len(effects),
                                'descriptive_overlap_rate_ratio': (len(shared) / len(lead)) / (len(genes) / len(effects)),
                                'shared_genes': ';'.join(shared)})
        for gene in sorted(lead):
            legacy_genes.append({'source_pathway': r.pathway, 'gene_name': gene})
    pd.DataFrame(legacy_rows).to_csv(OUT / 'original_epigenetic_leading_edge_reference_overlap.csv', index=False)
    pd.DataFrame(legacy_genes).merge(reference.reset_index(), on='gene_name', validate='m:1').to_csv(
        OUT / 'original_epigenetic_leading_edge_annotations.csv', index=False)
    # Covariate imbalance is descriptive; it does not identify an artifact by itself.
    cov = []
    for name, genes in sets.items():
        for flag, d in [('in_set', effects.loc[effects.index.isin(genes)]),
                        ('complement', effects.loc[~effects.index.isin(genes)])]:
            cov.append({'pathway': name, 'group': flag, 'n': len(d), 'mean_logCPM_median': d.mean_logCPM.median(),
                        'feature_length_bp_median': d.feature_length_bp.median(),
                        'fraction_CPM_below_2': d.mean_CPM.lt(2).mean()})
    pd.DataFrame(cov).to_csv(OUT / 'covariate_imbalance.csv', index=False)
    validation = {'status': 'passed', 'input_hashes_checked': len(manifest),
                  'primary_ES_and_leading_edge_checks': len(sets), 'all_BH_BY_and_pooled_BH_reproduced': True,
                  'fgsea_rows': len(fg), 'camera_rows': len(camera), 'gene_identity_and_size_checks': True,
                  'unestimable_fgsea_sensitivity_tests': int(bad.sum()),
                  'original_statistics_unchanged': True}
    if args.primary_only:
        validation['status'] = 'primary_checks_passed'
    (PROV / ('primary_validation.json' if args.primary_only else 'validation.json')).write_text(json.dumps(validation, indent=2))
    print(summary[['pathway', 'size', 'NES', 'BH_custom_20', 'camera_fixed_BH_custom_20',
                   'camera_estimated_BH_custom_20', 'combined_mean', 'combined_fraction_positive']].to_string(index=False))
    print(json.dumps(validation))


if __name__ == '__main__':
    main()
