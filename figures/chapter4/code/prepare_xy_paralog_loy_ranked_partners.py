"""Select genome-wide top-ten positive partners within LOY male models."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from prepare_xy_ranked_partners_by_y_status import sha256

THESIS = Path(__file__).resolve().parents[3]
OUT = THESIS / 'figures/chapter4'
BASE = OUT / 'source_data/xy_paralog_codependency'
PREVIOUS = OUT / 'source_data/xy_paralog_ranked_partners_y_status'
DATA = OUT / 'source_data/xy_paralog_loy_ranked_partners'
SOURCE = Path('A:/DEPMAP/DepMap/CRISPRGeneEffect.csv')
MIN_LOY_N = 100


def cross_cor(x, y):
    """Pairwise-complete Pearson r; mask unestimable values, not low-coverage ranks."""
    mx = np.isfinite(x).astype(float)
    my = np.isfinite(y).astype(float)
    x, y = np.nan_to_num(x), np.nan_to_num(y)
    n = mx.T @ my
    sx, sy = x.T @ my, mx.T @ y
    sxx, syy = (x*x).T @ my, mx.T @ (y*y)
    with np.errstate(invalid='ignore', divide='ignore'):
        numerator = x.T @ y - sx*sy/n
        denominator = np.sqrt(np.maximum(sxx-sx*sx/n, 0) * np.maximum(syy-sy*sy/n, 0))
        r = numerator / denominator
    r[(n < 3) | (denominator <= 1e-12)] = np.nan
    return np.clip(r, -1, 1), n.astype(int)


def main():
    DATA.mkdir(parents=True, exist_ok=True)
    previous = json.loads((PREVIOUS/'analysis_manifest.json').read_text())
    digest = sha256(SOURCE)
    assert digest == previous['sources'][str(SOURCE)]
    old = pd.read_csv(BASE/'x_member_top10_long.csv')
    core = ['DDX3X', 'EIF1AX', 'RPS4X', 'ZFX']
    seeds = core + sorted(set(old.seed)-set(core))
    assert len(seeds) == 16
    print('Reading frozen 25Q3 genome-wide gene effects...', flush=True)
    ge = pd.read_csv(SOURCE, index_col=0)
    ge.columns = ge.columns.str.replace(r' \(\d+\)$', '', regex=True)
    assert ge.columns.is_unique and ge.index.is_unique
    cohort = pd.read_csv(PREVIOUS/'model_cohort_and_exclusions.csv').set_index('ModelID')
    assert ge.index.tolist() == cohort.index.tolist()
    selected = cohort.loc[cohort.included]
    assert selected.Sex.eq('Male').all() and selected.Sex_calls.eq('Male').all()
    cohort.to_csv(DATA/'model_cohort_and_exclusions.csv', index_label='ModelID')
    arrays = {}
    for group in ['XY', 'LOY']:
        values = ge.loc[selected.index[selected.group.eq(group)]]
        assert len(values) == previous['model_counts'][group]
        r, n = cross_cor(values[seeds].to_numpy(), values.to_numpy())
        arrays[f'{group.lower()}_r'], arrays[f'{group.lower()}_n'] = r, n
        print(f'{group}: {len(values)} models, {len(seeds)} seeds x {ge.shape[1]} partners', flush=True)
    np.savez_compressed(DATA/'genomewide_correlations.npz', seeds=np.array(seeds),
                        genes=ge.columns.to_numpy(dtype=str), **arrays)
    tables, coverage, comparison_rows, check_errors, old_errors = [], [], [], [], []
    gene_position = {g: i for i, g in enumerate(ge.columns)}
    old_stratified = pd.read_csv(PREVIOUS/'ranked_partners_by_y_status_long.csv')
    for pair in old_stratified.itertuples():
        i, j = seeds.index(pair.seed), gene_position[pair.partner]
        assert arrays[f'{pair.group.lower()}_n'][i, j] == pair.n_models
        old_errors.append(abs(arrays[f'{pair.group.lower()}_r'][i, j]-pair.pearson_r))
    for i, seed in enumerate(seeds):
        table = pd.DataFrame(dict(seed=seed, partner=ge.columns,
                                 loy_pearson_r=arrays['loy_r'][i], loy_n_models=arrays['loy_n'][i],
                                 xy_pearson_r=arrays['xy_r'][i], xy_n_models=arrays['xy_n'][i]))
        eligible = table.partner.ne(seed) & table.loy_n_models.ge(MIN_LOY_N) & table.loy_pearson_r.notna()
        ranked = table.loc[eligible & table.loy_pearson_r.gt(0)].sort_values(
            ['loy_pearson_r', 'partner'], ascending=[False, True])
        top = ranked.head(10).copy()
        assert len(top) == 10 and top.partner.is_unique
        top.insert(2, 'rank_positive', np.arange(1, 11))
        tables.append(top)
        previous_set = set(old.loc[old.seed.eq(seed), 'partner'])
        coverage.append(dict(seed=seed, loy_seed_models=int(arrays['loy_n'][i, gene_position[seed]]),
                             eligible_nonself_partners=int(eligible.sum()),
                             eligible_positive_partners=len(ranked),
                             overlap_with_pooled_top10=len(set(top.partner) & previous_set),
                             selected_loy_n_min=int(top.loy_n_models.min()),
                             selected_loy_n_max=int(top.loy_n_models.max())))
        for pair in top.itertuples():
            for group in ['XY', 'LOY']:
                values = ge.loc[selected.index[selected.group.eq(group)], [seed, pair.partner]].dropna()
                expected_r = getattr(pair, f'{group.lower()}_pearson_r')
                expected_n = getattr(pair, f'{group.lower()}_n_models')
                assert len(values) == expected_n
                direct = float(np.corrcoef(values.to_numpy().T)[0, 1])
                assert np.isfinite(expected_r) and np.isfinite(direct)
                check_errors.append(abs(direct-expected_r))
                comparison_rows.append(dict(seed=seed, partner=pair.partner,
                                            rank_positive=pair.rank_positive, group=group,
                                            pearson_r=expected_r, n_models=expected_n,
                                            n_patients=int(selected.loc[values.index, 'PatientID'].nunique())))
    top = pd.concat(tables, ignore_index=True)
    assert len(top) == 160 and max(check_errors) < 1e-9 and max(old_errors) < 1e-9
    top.to_csv(DATA/'loy_top10_partners.csv', index=False)
    pd.DataFrame(coverage).to_csv(DATA/'selection_coverage.csv', index=False)
    pd.DataFrame(comparison_rows).to_csv(DATA/'selected_partners_by_group_long.csv', index=False)
    record = {
        'dependency_release': previous['dependency_release'],
        'y_status_release': previous['y_status_release'],
        'cross_release_join': 'Exact ModelID; unchanged 25Q3 dependency profiles and current Chapter 3 26Q1 primary cohort/calls.',
        'cohort_source': str(PREVIOUS/'model_cohort_and_exclusions.csv'),
        'model_counts': previous['model_counts'], 'patient_counts': previous['patient_counts'],
        'eligibility': previous['eligibility'], 'excluded_males': previous['excluded_males'],
        'gene_features': ge.shape[1], 'seed_count': len(seeds), 'displayed_pairs': len(top),
        'partner_selection': 'Top ten positive within-LOY Pearson correlations per X-linked seed, genome-wide, excluding self; r descending with gene-symbol ascending tie break.',
        'min_paired_loy_models': MIN_LOY_N,
        'coverage_rationale': 'Exploratory subset ranking requires 100 paired LOY models. The former 200-model pooled rule is not reused because two seeds have only 199 measured LOY models.',
        'xy_comparison': 'Same LOY-selected partners; XY r does not affect selection or order.',
        'estimation': 'Unadjusted, pairwise-complete Pearson r; at least three paired observations and nonconstant profiles for a finite estimate.',
        'by_group': {group: {'n_range': [int(top[f'{group.lower()}_n_models'].min()), int(top[f'{group.lower()}_n_models'].max())],
                              'r_range': [float(top[f'{group.lower()}_pearson_r'].min()), float(top[f'{group.lower()}_pearson_r'].max())]}
                     for group in ['XY','LOY']},
        'validation': {'old_stratified_values_reproduced': len(old_errors), 'max_old_r_error': max(old_errors),
                       'direct_selected_pair_checks': len(check_errors), 'max_direct_error': max(check_errors)},
        'interpretation': 'Descriptive selection on LOY r, not a differential-correlation or LOY-specificity test. LOY r is subject to selection inflation; lineage mix and donor relationships are not adjusted.',
        'sources': {str(SOURCE): digest, **{str(p): sha256(p) for p in [
            BASE/'x_member_top10_long.csv', PREVIOUS/'model_cohort_and_exclusions.csv',
            PREVIOUS/'analysis_manifest.json', PREVIOUS/'ranked_partners_by_y_status_long.csv']}}
    }
    (DATA/'analysis_manifest.json').write_text(json.dumps(record, indent=2), encoding='utf-8')
    print(json.dumps({k:v for k,v in record.items() if k!='sources'}, indent=2))
    print(pd.DataFrame(coverage).to_string(index=False))


if __name__ == '__main__':
    main()
