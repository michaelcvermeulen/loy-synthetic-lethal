"""Pearson co-dependencies in LOY and Y-retaining males for frozen partner pairs.

Retain the Public25Q3 gene effects used in the companion figure. Match the
current Chapter 3 (26Q1) primary cohort/calls by exact ModelID, recording this
cross-release annotation step rather than silently changing dependency data.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

THESIS = Path(__file__).resolve().parents[3]
OUT = THESIS / 'figures/chapter4'
BASE = OUT / 'source_data/xy_paralog_codependency'
DATA = OUT / 'source_data/xy_paralog_ranked_partners_y_status'
SOURCE = Path('A:/DEPMAP/DepMap/CRISPRGeneEffect.csv')
EXPORT = THESIS / 'research/depmap_loy_reanalysis_v1/loy_reanalysis_v1_export_2026-10-07'
COHORT = EXPORT / 'derived/loy_reanalysis_v1/cohort.csv'
CALLS = EXPORT / 'derived/loy_reanalysis_v1/y_status_calls.csv'


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def main():
    DATA.mkdir(parents=True, exist_ok=True)
    original = json.loads((BASE / 'analysis_manifest.json').read_text())
    expected = next(v['sha256'] for k, v in original['sources'].items()
                    if k.replace('\\', '/').endswith('/CRISPRGeneEffect.csv'))
    digest = sha256(SOURCE)
    assert digest == expected, 'Gene-effect source differs from frozen partner-selection analysis.'
    top = pd.read_csv(BASE / 'x_member_top10_long.csv')
    needed = set(top.seed) | set(top.partner)
    ge = pd.read_csv(SOURCE, usecols=lambda c: c.startswith('Unnamed:') or
                     c.rsplit(' (', 1)[0] in needed, index_col=0)
    ge.columns = ge.columns.str.replace(r' \(\d+\)$', '', regex=True)
    models = pd.read_csv(BASE / 'model_cohort.csv').set_index('ModelID')
    cohort = pd.read_csv(COHORT).set_index('ModelID')
    calls = pd.read_csv(CALLS).set_index('ModelID')
    assert all(df.index.is_unique for df in [ge, models, cohort, calls])
    assert ge.columns.is_unique and set(ge.index) == set(models.index)
    assert set(ge.columns) == needed
    assert cohort.Sex.eq('Male').all()
    assert cohort.OncotreePrimaryDisease.ne('Non-Cancerous').all()
    assert cohort.EngineeredModel.fillna('').eq('').all()
    assert cohort.in_primary.eq(cohort.y_call.notna()).all()
    pd.testing.assert_series_equal(cohort.y_call, calls.y_call.reindex(cohort.index), check_names=False)
    cols = ['Sex', 'PatientID', 'y_call', 'y_tier', 'in_primary', 'rna_call',
            'dna_call', 'dna_source', 'dna_partial_type', 'possible_female', 'call_stable']
    joined = models.reindex(ge.index).join(cohort[cols], rsuffix='_calls', validate='one_to_one')
    matched = joined.Sex_calls.notna()
    assert joined.loc[matched, 'Sex'].eq(joined.loc[matched, 'Sex_calls']).all()
    joined['group'] = joined.y_call.map({'Y_positive': 'XY', 'Y_negative': 'LOY'})
    included = joined.Sex.eq('Male') & joined.in_primary.eq(True) & joined.group.notna()
    joined['included'] = included
    joined['exclusion_reason'] = np.select(
        [included, joined.Sex.ne('Male'), ~matched],
        ['', 'Not recorded male in dependency cohort', 'Outside Chapter 3 eligible male cohort'],
        default='Unresolved primary Y-status call')
    joined.to_csv(DATA / 'model_cohort_and_exclusions.csv', index_label='ModelID')
    selected = joined.loc[included]
    rows, pooled_errors, independent_errors = [], [], []
    for pair in top.itertuples():
        pooled = ge[[pair.seed, pair.partner]].dropna()
        assert len(pooled) == pair.n_models
        pooled_errors.append(abs(pooled.corr().iloc[0, 1] - pair.pearson_r))
        for group in ['XY', 'LOY']:
            ids = selected.index[selected.group.eq(group)]
            values = ge.loc[ids, [pair.seed, pair.partner]].dropna()
            n = len(values)
            reason = ''
            if n < 3:
                r, reason = np.nan, 'Fewer than three paired models'
            elif not values.std().gt(0).all():
                r, reason = np.nan, 'Constant gene-effect profile'
            else:
                r = float(values.corr().iloc[0, 1])
                centered = values.to_numpy() - values.to_numpy().mean(axis=0)
                independent = float(centered[:, 0] @ centered[:, 1] /
                                    np.sqrt(np.sum(centered[:, 0]**2) * np.sum(centered[:, 1]**2)))
                independent_errors.append(abs(r - independent))
            rows.append(dict(seed=pair.seed, partner=pair.partner,
                             rank_positive=pair.rank_positive,
                             pooled_pearson_r=pair.pearson_r, pooled_n_models=pair.n_models,
                             group=group, pearson_r=r, n_models=n,
                             n_patients=int(selected.loc[values.index, 'PatientID'].nunique()),
                             missing_reason=reason))
    assert max(pooled_errors) < 1e-9 and max(independent_errors) < 1e-12
    result = pd.DataFrame(rows)
    assert len(result) == 320 and not result[['seed', 'partner', 'group']].duplicated().any()
    assert result.pearson_r.dropna().between(-1, 1).all()
    result.to_csv(DATA / 'ranked_partners_by_y_status_long.csv', index=False)
    selected.groupby(['group','OncotreeLineage']).size().rename('n_models').reset_index().to_csv(
        DATA / 'lineage_composition.csv', index=False)
    selected.groupby(['group','y_tier']).size().rename('n_models').reset_index().to_csv(
        DATA / 'call_tier_counts.csv', index=False)
    summary = {
        'dependency_release': original['release'], 'y_status_release': 'DepMap 26Q1; current Chapter 3 primary calls, export 2026-10-07',
        'cross_release_join': 'Exact ModelID; dependency values and pooled partner ranks remain Public25Q3. No Y-status reclassification was run.',
        'selection': 'Original 160 pooled top-ten pairs; no within-group reranking or new hit selection.',
        'correlation': 'Unadjusted pairwise-complete Pearson within each eligible male Y-status group.',
        'group_definition': {'XY': 'Primary Y_positive (Y-retaining); not a confirmed diploid XY karyotype.',
                             'LOY': 'Primary Y_negative (Y-depleted); includes original nonconcordant evidence tiers.'},
        'eligibility': 'Recorded male in both sources, cancer, not engineered, in the Chapter 3 primary cohort, with definite primary Y status and available Public25Q3 gene effects.',
        'model_counts': {str(k): int(v) for k, v in selected.group.value_counts().items()},
        'patient_counts': {str(k): int(v) for k, v in selected.groupby('group').PatientID.nunique().items()},
        'excluded_males': {str(k): int(v) for k, v in joined.loc[joined.Sex.eq('Male') & ~included, 'exclusion_reason'].value_counts().items()},
        'pooled_selection_min_models': 200,
        'within_group_selection': 'No reranking or 200-model re-filter. Estimability requires at least three pairs and nonconstant profiles; all counts retained.',
        'displayed_pairs': len(top), 'stratified_values': len(result),
        'by_group': {g: {'n_range': [int(t.n_models.min()), int(t.n_models.max())],
                          'r_range': [float(t.pearson_r.min()), float(t.pearson_r.max())],
                          'not_estimable': int(t.pearson_r.isna().sum()),
                          'pairs_below_200_models': int(t.n_models.lt(200).sum())}
                     for g, t in result.groupby('group')},
        'validation': {'max_pooled_r_error': max(pooled_errors), 'all_pooled_n_match': True,
                       'independent_centered_dot_checks': len(independent_errors),
                       'max_independent_r_error': max(independent_errors),
                       'no_matched_sex_disagreements': True, 'primary_calls_match_call_table': True},
        'interpretation': 'Descriptive model-level estimates; lineage mix and donor relationships are not adjusted. Differences in r are not an interaction test or evidence of an LOY-specific dependency.',
        'sources': {str(SOURCE): digest, **{str(p): sha256(p) for p in
                    [BASE/'x_member_top10_long.csv', BASE/'model_cohort.csv', BASE/'analysis_manifest.json', COHORT, CALLS]}}
    }
    (DATA / 'analysis_manifest.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    print(json.dumps({k: v for k, v in summary.items() if k != 'sources'}, indent=2))


if __name__ == '__main__':
    main()
