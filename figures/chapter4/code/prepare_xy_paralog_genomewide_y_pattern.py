"""Apply the frozen DDX3X descriptive co-dependency pattern to 16 X gametologs.

The selection is an effect-size rule, not a differential-correlation test.
Missing anchor and partner measurements are handled jointly, without imputation.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from prepare_ddx3x_partner_adjustment import nuisance_design
from prepare_xy_ranked_partners_by_y_status import sha256
from prepare_ddx3x_genomewide_y_pattern import (
    GROUPS, METHODS, MIN_MODELS, MIN_POSITIVE, MAX_ABS_XY, MIN_DELTA,
    select_pattern,
)

THESIS = Path(__file__).resolve().parents[3]
OUT = THESIS / 'figures/chapter4'
DATA = OUT / 'source_data/xy_paralog_genomewide_y_pattern'
PREVIOUS = OUT / 'source_data/ddx3x_genomewide_y_pattern'
COHORT = OUT / 'source_data/ddx3x_partner_adjustment'
ANCHORS = OUT / 'source_data/xy_paralog_loy_ranked_partners/selection_coverage.csv'
SOURCE = Path('A:/DEPMAP/DepMap/CRISPRGeneEffect.csv')


def correlations_for_group(ge, meta, seed):
    """Fit both genes on identical complete cases, grouped by missingness mask."""
    assert ge.index.equals(meta.index)
    anchor_finite = np.isfinite(ge[seed].to_numpy(float))
    ge, meta = ge.loc[anchor_finite], meta.loc[anchor_finite]
    anchor = ge[seed].to_numpy(float)
    y = ge.to_numpy(float)
    finite = np.isfinite(y)
    _, patterns = np.unique(np.packbits(finite, axis=0).T, axis=0,
                            return_inverse=True)
    arrays = {a: np.full(y.shape[1], np.nan) for a in METHODS}
    counts = finite.sum(axis=0)
    donor_counts = np.zeros(y.shape[1], int)
    rank = np.zeros(y.shape[1], int)
    max_error = 0.
    for pattern in np.unique(patterns):
        columns = np.flatnonzero(patterns == pattern)
        mask = finite[:, columns[0]]
        assert (finite[:, columns] == mask[:, None]).all()
        frame = meta.iloc[np.flatnonzero(mask)]
        donor_counts[columns] = frame.PatientID.nunique()
        if len(frame) < MIN_MODELS:
            continue
        pair = np.column_stack([anchor[mask], y[np.ix_(mask, columns)]])
        for method in METHODS:
            x = nuisance_design(frame, method).to_numpy(float)
            beta, _, q, _ = np.linalg.lstsq(x, pair, rcond=None)
            residual = pair - x @ beta
            assert len(frame) - q > 3
            if method == 'lineage_efficacy':
                rank[columns] = q
            max_error = max(max_error, float(np.max(np.abs(x.T @ residual))))
            num = residual[:, 0] @ residual[:, 1:]
            den = np.sqrt(np.sum(residual[:, 0]**2)
                          * np.sum(residual[:, 1:]**2, axis=0))
            r = np.divide(num, den, out=np.full(len(columns), np.nan),
                          where=den > 1e-10)
            arrays[method][columns] = np.clip(r, -1, 1)
    assert max_error < 1e-7
    return arrays, counts, donor_counts, rank, {
        'anchor_models': len(meta),
        'complete_case_patterns': int(len(np.unique(patterns))),
        'max_normal_equation_error': max_error,
    }


def main():
    DATA.mkdir(parents=True, exist_ok=True)
    prior_manifest = json.loads((PREVIOUS / 'manifest.json').read_text())
    digest = sha256(SOURCE)
    assert digest == prior_manifest['sources'][str(SOURCE)]
    seeds = pd.read_csv(ANCHORS).seed.tolist()
    assert len(seeds) == 16 and len(set(seeds)) == 16
    meta = pd.read_csv(COHORT / 'model_values.csv').set_index('ModelID')
    ge = pd.read_csv(SOURCE, index_col=0)
    original_names = ge.columns.tolist()
    ge.columns = ge.columns.str.replace(r' \(\d+\)$', '', regex=True)
    assert ge.columns.is_unique and ge.index.is_unique
    assert set(seeds) <= set(ge.columns)
    frames, summary, sensitivity, checks = [], [], [], {}
    independent_errors, previous_errors = [], []
    for seed in seeds:
        result = pd.DataFrame({'anchor': seed, 'gene': ge.columns,
                               'source_feature': original_names})
        checks[seed] = {}
        for group in GROUPS:
            frame = meta.loc[meta.group.eq(group)]
            values = ge.loc[frame.index]
            arrays, n, nd, q, check = correlations_for_group(values, frame, seed)
            prefix = group.lower()
            result[f'{prefix}_r'] = arrays['lineage_efficacy']
            result[f'{prefix}_raw_r'] = arrays['raw']
            result[f'{prefix}_lineage_r'] = arrays['lineage']
            result[f'{prefix}_n_models'] = n
            result[f'{prefix}_n_patients'] = nd
            result[f'{prefix}_design_rank'] = q
            checks[seed][group] = check
        result['eligible'] = (
            result.gene.ne(seed)
            & result[[f'{g.lower()}_n_models' for g in GROUPS]].ge(
                MIN_MODELS).all(axis=1)
            & result[['xy_r', 'loy_r', 'female_r']].notna().all(axis=1))
        result['loy_minus_xy'] = result.loy_r - result.xy_r
        result['female_minus_xy'] = result.female_r - result.xy_r
        result['min_y_absent_r'] = result[['loy_r', 'female_r']].min(axis=1)
        result['min_delta'] = result[['loy_minus_xy', 'female_minus_xy']].min(axis=1)
        result['selected'] = select_pattern(result)
        selected = result.loc[result.selected].sort_values(
            ['min_y_absent_r', 'min_delta', 'gene'],
            ascending=[False, False, True]).copy()
        selected.insert(0, 'within_anchor_rank', np.arange(1, len(selected) + 1))
        selected.to_csv(DATA / f'{seed}_selected_partners.csv', index=False)
        result.to_csv(DATA / f'{seed}_all_gene_correlations.csv.gz', index=False,
                      compression={'method': 'gzip', 'mtime': 0})
        frames.append(result)
        row = {
            'anchor': seed,
            'eligible_nonself_partners': int(result.eligible.sum()),
            'selected_partners': len(selected),
            'selected_percent': 100 * len(selected) / result.eligible.sum(),
            'top_five_partners': '; '.join(selected.gene.head(5)),
        }
        for group in GROUPS:
            prefix = group.lower()
            row[f'{prefix}_anchor_models'] = checks[seed][group]['anchor_models']
            row[f'{prefix}_selected_n_min'] = (
                int(selected[f'{prefix}_n_models'].min()) if len(selected) else None)
            row[f'{prefix}_selected_n_max'] = (
                int(selected[f'{prefix}_n_models'].max()) if len(selected) else None)
        summary.append(row)
        for positive in [.20, .25, .30]:
            for abs_xy in [.10, .15, .20]:
                for delta in [.15, .20, .25]:
                    keep = select_pattern(result, positive, abs_xy, delta)
                    sensitivity.append(dict(
                        anchor=seed, min_positive_r=positive,
                        max_abs_xy_r=abs_xy, min_delta=delta,
                        selected_partners=int(keep.sum()),
                        primary=(positive == MIN_POSITIVE and abs_xy == MAX_ABS_XY
                                 and delta == MIN_DELTA)))
        # Single-pair refits independently check batching and missing-anchor handling.
        for gene in selected.gene:
            for group in GROUPS:
                frame = meta.loc[meta.group.eq(group)]
                pair = ge.loc[frame.index, [seed, gene]].dropna()
                x = nuisance_design(frame.loc[pair.index],
                                    'lineage_efficacy').to_numpy(float)
                residual = pair.to_numpy() - x @ np.linalg.lstsq(
                    x, pair.to_numpy(), rcond=None)[0]
                independent = np.corrcoef(residual.T)[0, 1]
                saved = selected.loc[selected.gene.eq(gene),
                                     f'{group.lower()}_r'].iloc[0]
                independent_errors.append(abs(independent - saved))
        if seed == 'DDX3X':
            old = pd.read_csv(PREVIOUS / 'all_gene_correlations.csv')
            assert result.gene.tolist() == old.gene.tolist()
            assert result.selected.tolist() == old.selected.tolist()
            assert result.eligible.tolist() == old.eligible.tolist()
            for col in [c for c in result if c.endswith('_r')]:
                a, b = result[col].to_numpy(), old[col].to_numpy()
                np.testing.assert_allclose(a, b, rtol=0, atol=1e-10, equal_nan=True)
                previous_errors.extend(np.abs(a[np.isfinite(a)] - b[np.isfinite(a)]))
            assert len(selected) == 31
        print(f'{seed}: {len(selected)} / {int(result.eligible.sum())}; '
              f'top: {row["top_five_partners"]}', flush=True)
    assert max(independent_errors, default=0) < 1e-10
    summary = pd.DataFrame(summary)
    summary.to_csv(DATA / 'anchor_summary.csv', index=False)
    pd.DataFrame(sensitivity).to_csv(DATA / 'threshold_sensitivity_counts.csv',
                                     index=False)
    all_results = pd.concat(frames, ignore_index=True)
    selected = all_results.loc[all_results.selected].copy()
    selected = selected.sort_values(
        ['anchor', 'min_y_absent_r', 'min_delta', 'gene'],
        ascending=[True, False, False, True])
    selected['within_anchor_rank'] = selected.groupby('anchor').cumcount() + 1
    selected.to_csv(DATA / 'all_selected_partners.csv', index=False)
    # Retain partner recurrence and DDX3X overlap without implying independence.
    partner_anchors = selected.groupby('gene').anchor.agg(
        lambda x: '; '.join(sorted(x))).rename('anchors').to_frame()
    partner_anchors['n_anchors'] = selected.groupby('gene').anchor.nunique()
    partner_anchors.sort_values(['n_anchors', 'gene'], ascending=[False, True]
                                ).to_csv(DATA / 'partner_recurrence.csv')
    ddx = set(selected.loc[selected.anchor.eq('DDX3X'), 'gene'])
    overlaps = []
    for seed in seeds:
        partners = set(selected.loc[selected.anchor.eq(seed), 'gene'])
        overlaps.append(dict(anchor=seed, n_overlap_with_ddx3x=len(partners & ddx),
                             shared_genes='; '.join(sorted(partners & ddx))))
    pd.DataFrame(overlaps).to_csv(DATA / 'ddx3x_partner_overlap.csv', index=False)
    manifest = dict(
        analysis_date='2026-10-20',
        anchors=seeds,
        dependency_release=prior_manifest['dependency_release'],
        annotation_release=prior_manifest['annotation_release'],
        groups=prior_manifest['groups'],
        cohort='Identical saved model IDs and covariates to the DDX3X analysis.',
        adjustment=prior_manifest['adjustment'],
        min_paired_models_each_group=MIN_MODELS,
        selection=prior_manifest['selection'],
        missingness='Subset each group to finite anchor measurements first, then batch '
                    'partners by identical complete-case masks. Fit both genes on '
                    'the same pair-complete models. No imputation.',
        summary=summary.to_dict('records'),
        display_labels=prior_manifest['display_labels'],
        karyotype_note=prior_manifest['karyotype_note'],
        inference=prior_manifest['inference'],
        sensitivity='27 combinations per anchor; all combinations retained, '
                    'including the unchanged primary DDX3X rule.',
        validation={
            'ddx3x_selected_31_exactly_reproduced': True,
            'max_ddx3x_previous_coefficient_error': max(previous_errors, default=0),
            'independent_selected_pair_checks': len(independent_errors),
            'max_selected_pair_error': max(independent_errors, default=0),
            'group_checks': checks,
        },
        sources={str(SOURCE): digest, **{str(p): sha256(p) for p in [
            COHORT / 'model_values.csv', PREVIOUS / 'manifest.json',
            PREVIOUS / 'all_gene_correlations.csv', ANCHORS,
            Path(__file__), Path(__file__).with_name('prepare_ddx3x_partner_adjustment.py')
        ]}})
    (DATA / 'manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    print(summary[['anchor', 'selected_partners', 'eligible_nonself_partners',
                   'top_five_partners']].to_string(index=False))
    print('Validation:', {k: v for k, v in manifest['validation'].items()
                           if k != 'group_checks'}, flush=True)


if __name__ == '__main__':
    main()
