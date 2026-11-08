"""Recompute sex-stratified Pearson r for the frozen pooled top-ten pairs."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

THESIS = Path(__file__).resolve().parents[3]
DATA = THESIS / 'figures/chapter4/source_data/xy_paralog_codependency'
SOURCE = Path('A:/DEPMAP/DepMap/CRISPRGeneEffect.csv')


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def main():
    original = json.loads((DATA / 'analysis_manifest.json').read_text())
    expected = next(v['sha256'] for k, v in original['sources'].items()
                    if k.replace('\\', '/').endswith('/CRISPRGeneEffect.csv'))
    digest = sha256(SOURCE)
    assert digest == expected, 'Gene-effect source differs from the frozen pooled analysis.'
    top = pd.read_csv(DATA / 'x_member_top10_long.csv')
    needed = set(top.seed) | set(top.partner)
    ge = pd.read_csv(SOURCE, usecols=lambda c: c.startswith('Unnamed:') or
                     c.rsplit(' (', 1)[0] in needed, index_col=0)
    ge.columns = ge.columns.str.replace(r' \(\d+\)$', '', regex=True)
    cohort = pd.read_csv(DATA / 'model_cohort.csv').set_index('ModelID')
    assert ge.index.is_unique and cohort.index.is_unique and ge.columns.is_unique
    assert set(ge.index) == set(cohort.index)
    cohort = cohort.reindex(ge.index)
    assert cohort.Sex.notna().all() and set(ge.columns) == needed
    rows, pooled_errors, independent_errors = [], [], []
    for pair in top.itertuples():
        pooled = ge[[pair.seed, pair.partner]].dropna()
        assert len(pooled) == pair.n_models
        pooled_errors.append(abs(pooled.corr().iloc[0, 1] - pair.pearson_r))
        for sex in ['Male', 'Female']:
            values = ge.loc[cohort.Sex.eq(sex), [pair.seed, pair.partner]].dropna()
            n = len(values)
            assert n >= 3 and values.std().gt(0).all()
            r = float(values.corr().iloc[0, 1])
            centered = values.to_numpy() - values.to_numpy().mean(axis=0)
            independent = float(centered[:, 0] @ centered[:, 1] /
                                np.sqrt(np.sum(centered[:, 0]**2) * np.sum(centered[:, 1]**2)))
            independent_errors.append(abs(r - independent))
            rows.append(dict(seed=pair.seed, partner=pair.partner,
                             rank_positive=pair.rank_positive,
                             pooled_pearson_r=pair.pearson_r, pooled_n_models=pair.n_models,
                             sex=sex, pearson_r=r, n_models=n))
    assert max(pooled_errors) < 1e-9
    assert max(independent_errors) < 1e-12
    result = pd.DataFrame(rows)
    assert len(result) == 320 and not result[['seed', 'partner', 'sex']].duplicated().any()
    assert result.pearson_r.between(-1, 1).all()
    result.to_csv(DATA / 'ranked_partners_by_sex_long.csv', index=False)
    summary = {
        'release': original['release'],
        'selection': 'The original 160 pooled top-ten pairs; no sex-specific reranking or new selection.',
        'correlation': 'Unadjusted pairwise-complete Pearson within recorded Male or Female models.',
        'pooled_selection_min_models': 200,
        'within_sex_selection': 'None. All frozen pairs are estimable and retained; pairwise n is displayed.',
        'unknown_sex': 'Excluded from stratified values only; retained in the original pooled selection.',
        'model_counts': {str(k): int(v) for k, v in cohort.Sex.value_counts().items()},
        'displayed_pairs': len(top), 'stratified_values': len(result),
        'by_sex': {sex: {'n_range': [int(g.n_models.min()), int(g.n_models.max())],
                          'r_range': [float(g.pearson_r.min()), float(g.pearson_r.max())],
                          'pairs_below_200_models': int(g.n_models.lt(200).sum())}
                   for sex, g in result.groupby('sex')},
        'validation': {'max_pooled_r_error': max(pooled_errors),
                       'all_pooled_n_match': True,
                       'independent_centered_dot_checks': len(result),
                       'max_independent_r_error': max(independent_errors)},
        'interpretation': 'Descriptive correlations; lineage composition is not adjusted. Differences in r are not a sex-interaction test. Recorded sex is not measured Y status.',
        'sources': {str(SOURCE): digest,
                    **{str(DATA / f): sha256(DATA / f) for f in
                       ['x_member_top10_long.csv', 'model_cohort.csv', 'analysis_manifest.json']}}
    }
    (DATA / 'ranked_partners_by_sex_manifest.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    print(json.dumps({k: v for k, v in summary.items() if k != 'sources'}, indent=2))


if __name__ == '__main__':
    main()
