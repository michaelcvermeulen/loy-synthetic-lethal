"""Adjusted mean-dependency contrasts for DDX3X and its 31 selected partners.

Primary correction spans all eligible genes and all three pairwise contrasts.
FWL batch fits use donor-clustered CR1 uncertainty. No correlation re-selection.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import scipy
from scipy import linalg, stats
from prepare_ddx3x_partner_adjustment import nuisance_design
from prepare_xy_ranked_partners_by_y_status import sha256

THESIS = Path(__file__).resolve().parents[3]
OUT = THESIS / 'figures/chapter4'
DATA = OUT / 'source_data/ddx3x_partner_differential_dependency'
COHORT = OUT / 'source_data/ddx3x_partner_adjustment'
PATTERN = OUT / 'source_data/ddx3x_genomewide_y_pattern'
SOURCE = Path('A:/DEPMAP/DepMap/CRISPRGeneEffect.csv')
CONTRASTS = [('LOY_minus_XY', 'LOY', 'XY'),
             ('Female_minus_XY', 'Female', 'XY'),
             ('LOY_minus_Female', 'LOY', 'Female')]
PRIMARY = 'lineage_efficacy'
MIN_MODELS = 100


def fit_pair(y, frame, target, method):
    """Batch responses with a shared complete-case design; return CR1 estimates."""
    z = nuisance_design(frame, method).to_numpy(float)
    g = frame.group.eq(target).to_numpy(float)
    bz, _, rank_z, _ = np.linalg.lstsq(z, np.column_stack([g, y]), rcond=None)
    h = g - z @ bz[:, 0]
    ss = h @ h
    assert ss > 1e-8, 'Group coefficient is not estimable'
    effect = (h @ y) / ss
    residual = y - z @ bz[:, 1:] - h[:, None] * effect
    n, k = len(frame), int(rank_z + 1)
    codes, donors = pd.factorize(frame.PatientID, sort=True)
    ng = len(donors)
    score = np.zeros((ng, y.shape[1]))
    np.add.at(score, codes, h[:, None] * residual)
    correction = ng / (ng - 1) * (n - 1) / (n - k)
    se = np.sqrt(correction * np.sum(score**2, axis=0)) / ss
    assert np.isfinite(se).all() and (se > 0).all()
    t = effect / se
    p = 2 * stats.t.sf(np.abs(t), df=ng - 1)
    margin = stats.t.ppf(.975, df=ng - 1) * se
    error = float(max(np.max(np.abs(z.T @ residual)),
                      np.max(np.abs(g @ residual))))
    assert error < 1e-7
    return dict(effect=effect, se_cluster=se, ci_low=effect - margin,
                ci_high=effect + margin, t_stat=t, p_value=p,
                n_models=n, n_patients=ng, model_rank=k,
                inferential_df=ng - 1, cr1_factor=correction,
                max_normal_equation_error=error)


def fit_contrast(ge, meta, contrast, target, reference, method, shared_only=False):
    frame = meta.loc[meta.group.isin([target, reference])].copy()
    if shared_only:
        counts = pd.crosstab(frame.OncotreeLineage, frame.group)
        keep = counts.index[counts[[target, reference]].ge(5).all(axis=1)]
        frame = frame.loc[frame.OncotreeLineage.isin(keep)]
    values = ge.loc[frame.index]
    finite = np.isfinite(values.to_numpy(float))
    _, patterns = np.unique(np.packbits(finite, axis=0).T, axis=0,
                            return_inverse=True)
    records = []
    for pattern in np.unique(patterns):
        columns = np.flatnonzero(patterns == pattern)
        mask = finite[:, columns[0]]
        sub = frame.iloc[np.flatnonzero(mask)]
        n_target = int(sub.group.eq(target).sum())
        n_ref = int(sub.group.eq(reference).sum())
        assert n_target >= MIN_MODELS and n_ref >= MIN_MODELS
        y = values.to_numpy(float)[np.ix_(mask, columns)]
        fitted = fit_pair(y, sub, target, method)
        target_values = y[sub.group.eq(target)]
        ref_values = y[sub.group.eq(reference)]
        result = pd.DataFrame({
            'gene': values.columns[columns], 'contrast': contrast,
            'target_group': target, 'reference_group': reference,
            'analysis': 'shared_lineages' if shared_only else method,
            'n_target': n_target, 'n_reference': n_ref,
            'raw_mean_target': target_values.mean(axis=0),
            'raw_mean_reference': ref_values.mean(axis=0),
            'raw_sd_target': target_values.std(axis=0, ddof=1),
            'raw_sd_reference': ref_values.std(axis=0, ddof=1),
            **fitted})
        records.append(result)
    return pd.concat(records, ignore_index=True)


def independent_pair_check(y, frame, target, method):
    """Full-design QR/OLS and explicit cluster sandwich, independent of FWL."""
    z = nuisance_design(frame, method).to_numpy(float)
    _, _, piv = linalg.qr(z, mode='economic', pivoting=True)
    rz = np.linalg.matrix_rank(z)
    x = np.column_stack([z[:, piv[:rz]], frame.group.eq(target).to_numpy(float)])
    q, r = np.linalg.qr(x, mode='reduced')
    beta = linalg.solve_triangular(r, q.T @ y)
    residual = y - x @ beta
    inv = linalg.solve_triangular(r, np.eye(r.shape[0]))
    bread = inv @ inv.T
    meat = np.zeros((x.shape[1], x.shape[1]))
    for donor in frame.PatientID.unique():
        mask = frame.PatientID.eq(donor).to_numpy()
        score = x[mask].T @ residual[mask]
        meat += np.outer(score, score)
    n, k, ng = len(frame), x.shape[1], frame.PatientID.nunique()
    vcov = bread @ meat @ bread * ng / (ng - 1) * (n - 1) / (n - k)
    return float(beta[-1]), float(np.sqrt(vcov[-1, -1]))


def main():
    DATA.mkdir(parents=True, exist_ok=True)
    prior = json.loads((PATTERN / 'manifest.json').read_text())
    protocol = dict(
        analysis_date='2026-10-20',
        question='Do the 31 correlation-selected DDX3X partners differ in mean dependency?',
        contrasts=CONTRASTS, reference_gene='DDX3X',
        primary='Separate pooled two-group OLS fits: gene effect ~ Oncotree lineage '
                '+ screen efficacy + library/platform composition + group.',
        group_effect='Target minus reference in original Chronos gene-effect units; '
                     'negative means stronger target-group dependency.',
        uncertainty='PatientID-clustered sandwich CR1; finite sample factor '
                    'G/(G-1)*(n-1)/(n-k); two-sided t with G-1 degrees of freedom.',
        multiplicity='Primary BH jointly across all eligible genes and three contrasts; '
                     'BY sensitivity for the same family. BH across the 32 x 3 displayed '
                     'tests is exploratory and does not replace the genome-wide result.',
        eligibility='Existing 17,786 eligible nonself DDX3X partners plus DDX3X; '
                    'at least 100 measured models in each original group.',
        sensitivities=['Unadjusted group-only', 'Lineage only',
                       'Primary covariates restricted to lineages with >=5 models '
                       'in each group of a contrast'],
        excluded_adjustment='Do not adjust for DDX3X dependency itself. Do not use '
                            'within-group-centred residuals to estimate group means.',
        interpretation='Observational post hoc analysis; selected subset is not an '
                       'independent replication. Nonsignificance is not equivalence.',
        packages=dict(numpy=np.__version__, pandas=pd.__version__, scipy=scipy.__version__),
    )
    protocol_path = DATA / 'protocol.json'
    if protocol_path.exists():
        assert json.loads(protocol_path.read_text()) == protocol
    else:
        protocol_path.write_text(json.dumps(protocol, indent=2), encoding='utf-8')
    digest = sha256(SOURCE)
    assert digest == prior['sources'][str(SOURCE)]
    meta = pd.read_csv(COHORT / 'model_values.csv').set_index('ModelID')
    assert meta.index.is_unique and meta.PatientID.notna().all()
    old = pd.read_csv(PATTERN / 'all_gene_correlations.csv')
    universe = old.loc[old.eligible | old.gene.eq('DDX3X'), 'gene'].tolist()
    partners = pd.read_csv(PATTERN / 'selected_partners.csv').gene.tolist()
    panel = ['DDX3X'] + partners
    assert len(universe) == 17787 and len(panel) == 32
    ge = pd.read_csv(SOURCE, index_col=0)
    ge.columns = ge.columns.str.replace(r' \(\d+\)$', '', regex=True)
    assert ge.columns.is_unique
    ge = ge.loc[meta.index, universe]
    meta.to_csv(DATA / 'cohort_and_covariates.csv')
    ge[panel].to_csv(DATA / 'panel_model_gene_effects.csv')
    pd.crosstab(meta.OncotreeLineage, meta.group).to_csv(DATA / 'lineage_composition.csv')
    cross_ids = meta.groupby('PatientID').group.nunique()
    meta.loc[meta.PatientID.isin(cross_ids.index[cross_ids.gt(1)])].to_csv(
        DATA / 'donors_spanning_groups.csv')
    primary_frames, sensitivity_frames = [], []
    for contrast, target, reference in CONTRASTS:
        primary_frames.append(fit_contrast(ge, meta, contrast, target, reference, PRIMARY))
        print('Completed genome-wide', contrast, len(primary_frames[-1]), flush=True)
        for method, shared in [('raw', False), ('lineage', False), (PRIMARY, True)]:
            sensitivity_frames.append(fit_contrast(
                ge[panel], meta, contrast, target, reference, method, shared_only=shared))
    primary = pd.concat(primary_frames, ignore_index=True)
    assert len(primary) == 53361 and primary.p_value.notna().all()
    primary['fdr_bh_genomewide_joint'] = stats.false_discovery_control(
        primary.p_value.to_numpy(), method='bh')
    primary['fdr_by_genomewide_joint'] = stats.false_discovery_control(
        primary.p_value.to_numpy(), method='by')
    primary['is_partner'] = primary.gene.isin(partners)
    primary['is_reference_ddx3x'] = primary.gene.eq('DDX3X')
    primary.to_csv(DATA / 'genomewide_primary_results.csv.gz', index=False,
                   compression={'method': 'gzip', 'mtime': 0})
    display = primary.loc[primary.gene.isin(panel)].copy()
    display['fdr_bh_panel_joint_exploratory'] = stats.false_discovery_control(
        display.p_value.to_numpy(), method='bh')
    display['display_rank'] = display.gene.map({g: i for i, g in enumerate(panel)})
    display = display.sort_values(['contrast', 'display_rank'])
    display.to_csv(DATA / 'panel_primary_results.csv', index=False)
    sensitivity = pd.concat(sensitivity_frames, ignore_index=True)
    sensitivity['fdr_bh_panel_joint_within_analysis'] = sensitivity.groupby(
        'analysis').p_value.transform(lambda p: stats.false_discovery_control(p.to_numpy(), method='bh'))
    sensitivity.to_csv(DATA / 'panel_sensitivity_results.csv', index=False)
    # Independent validation for every displayed primary fit.
    errors = []
    for row in display.itertuples():
        frame = meta.loc[meta.group.isin([row.target_group, row.reference_group])]
        vals = ge.loc[frame.index, row.gene]
        frame, vals = frame.loc[vals.notna()], vals.loc[vals.notna()]
        b, se = independent_pair_check(vals.to_numpy(float), frame, row.target_group, PRIMARY)
        errors.append(dict(gene=row.gene, contrast=row.contrast,
                           effect_error=abs(b - row.effect), se_error=abs(se - row.se_cluster)))
    checks = pd.DataFrame(errors)
    assert checks.effect_error.max() < 1e-10 and checks.se_error.max() < 1e-10
    checks.to_csv(DATA / 'independent_validation.csv', index=False)
    raw = sensitivity.loc[sensitivity.analysis.eq('raw')]
    assert np.max(np.abs(raw.effect - (raw.raw_mean_target - raw.raw_mean_reference))) < 1e-12
    # Independently verify BH ordering and reverse cumulative minima.
    p = primary.p_value.to_numpy()
    order = np.argsort(p, kind='stable')
    expected = np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    np.testing.assert_allclose(primary.fdr_bh_genomewide_joint.to_numpy()[order],
                               np.clip(expected, 0, 1), rtol=1e-12, atol=1e-15)
    summaries = []
    for contrast, target, reference in CONTRASTS:
        rows = display.loc[display.contrast.eq(contrast) & display.is_partner]
        passing = rows.loc[rows.fdr_bh_genomewide_joint.lt(.05)]
        summaries.append(dict(
            contrast=contrast, partners=31,
            genomewide_bh_lt_005=len(passing),
            stronger_in_target=int(passing.effect.lt(0).sum()),
            weaker_in_target=int(passing.effect.gt(0).sum()),
            genes='; '.join(passing.sort_values('p_value').gene),
            genomewide_by_lt_005=int(rows.fdr_by_genomewide_joint.lt(.05).sum()),
            panel_bh_lt_005=int(rows.fdr_bh_panel_joint_exploratory.lt(.05).sum())))
    pd.DataFrame(summaries).to_csv(DATA / 'contrast_summary.csv', index=False)
    manifest = dict(
        **protocol, dependency_release=prior['dependency_release'],
        annotation_release=prior['annotation_release'], groups=prior['groups'],
        primary_genes=len(universe), primary_tests=len(primary),
        panel_genes=panel, panel_tests=len(display), summaries=summaries,
        missingness='Complete cases per gene and contrast, no imputation; identical '
                    'measurement subsets across primary, raw and lineage-only fits.',
        validation=dict(
            independent_full_design_checks=len(checks),
            max_effect_error=float(checks.effect_error.max()),
            max_cluster_se_error=float(checks.se_error.max()),
            max_normal_equation_error=float(primary.max_normal_equation_error.max()),
            raw_effects_equal_mean_differences=True, independent_bh_check_passed=True),
        sources={str(SOURCE): digest, **{str(p): sha256(p) for p in [
            Path(__file__), COHORT / 'model_values.csv',
            PATTERN / 'all_gene_correlations.csv', PATTERN / 'selected_partners.csv',
            protocol_path]}},
    )
    (DATA / 'manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    print(pd.DataFrame(summaries).to_string(index=False))
    fields = ['gene', 'contrast', 'effect', 'ci_low', 'ci_high', 'p_value',
              'fdr_bh_genomewide_joint', 'fdr_bh_panel_joint_exploratory']
    print(display.loc[display.fdr_bh_genomewide_joint.lt(.05)
                       | display.gene.isin(['DDX3X', 'IGBP1', 'GBF1', 'IK']), fields].to_string(index=False))
    print('Validation', manifest['validation'])


if __name__ == '__main__':
    main()
