"""DDX3Y expression, DDX3X dependency and the fixed 31-partner pattern.

Frozen DepMap 25Q3 expression/effects; existing Chapter 3 model eligibility.
No re-selection of partners, expression thresholds or screen hits.
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
DATA = OUT / 'source_data/ddx3y_expression_dependency'
DEP = Path('A:/DEPMAP/DepMap')
COHORT = OUT / 'source_data/ddx3x_partner_adjustment/model_values.csv'
EFFECTS = OUT / 'source_data/ddx3x_partner_differential_dependency/panel_model_gene_effects.csv'
PARTNERS = OUT / 'source_data/ddx3x_genomewide_y_pattern/selected_partners.csv'
EXPRESSION = DEP / 'OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv'
PROFILES = DEP / 'OmicsProfiles.csv'


def test_term(y, term, z, patients):
    """OLS test of one continuous term; FWL and donor-clustered CR1."""
    y, term, z = np.asarray(y, float), np.asarray(term, float), np.asarray(z, float)
    if y.ndim == 1:
        y = y[:, None]
    coef, _, rank, _ = np.linalg.lstsq(z, np.column_stack([term, y]), rcond=None)
    h = term - z @ coef[:, 0]
    residual_y = y - z @ coef[:, 1:]
    ss = h @ h
    assert ss > 1e-9
    beta = h @ residual_y / ss
    residual = residual_y - h[:, None] * beta
    codes, donors = pd.factorize(patients, sort=True)
    n, k, ng = len(y), int(rank + 1), len(donors)
    scores = np.zeros((ng, y.shape[1]))
    np.add.at(scores, codes, h[:, None] * residual)
    factor = ng / (ng - 1) * (n - 1) / (n - k)
    se = np.sqrt(factor * (scores**2).sum(axis=0)) / ss
    p = 2 * stats.t.sf(np.abs(beta / se), ng - 1)
    margin = stats.t.ppf(.975, ng - 1) * se
    r = h @ residual_y / np.sqrt(ss * (residual_y**2).sum(axis=0))
    assert np.max(np.abs(z.T @ residual)) < 1e-7
    return dict(effect=beta, se_cluster=se, ci_low=beta - margin,
                ci_high=beta + margin, p_value=p, partial_r=r,
                n_models=n, n_patients=ng, design_rank=k), h, residual_y


def independent_check(y, term, z, patients):
    """Validate FWL with a full rank QR fit and explicit sandwich covariance."""
    _, _, piv = linalg.qr(z, mode='economic', pivoting=True)
    rank = np.linalg.matrix_rank(z)
    x = np.column_stack([z[:, piv[:rank]], term])
    q, r = np.linalg.qr(x, mode='reduced')
    beta = linalg.solve_triangular(r, q.T @ y)
    residual = y - x @ beta
    inv = linalg.solve_triangular(r, np.eye(r.shape[0]))
    bread = inv @ inv.T
    meat = np.zeros((x.shape[1], x.shape[1]))
    for donor in np.unique(patients):
        ix = np.asarray(patients) == donor
        score = x[ix].T @ residual[ix]
        meat += np.outer(score, score)
    n, k, ng = len(y), x.shape[1], len(np.unique(patients))
    vcov = bread @ meat @ bread * ng / (ng - 1) * (n - 1) / (n - k)
    return beta[-1], np.sqrt(vcov[-1, -1])


def bootstrap_r(y, expression, design, patients, seed=20261020, n_boot=1000):
    """Patient-level resampling; refit both nuisance regressions in each draw."""
    rng = np.random.default_rng(seed)
    codes, donors = pd.factorize(patients, sort=True)
    values = np.column_stack([expression, y])
    results = []
    for _ in range(n_boot):
        counts = rng.multinomial(len(donors), np.full(len(donors), 1 / len(donors)))
        w = np.sqrt(counts[codes])[:, None]
        zw, yw = design * w, values * w
        residual = yw - zw @ np.linalg.lstsq(zw, yw, rcond=None)[0]
        cross = residual.T @ residual
        results.append(cross[0, 1] / np.sqrt(cross[0, 0] * cross[1, 1]))
    return np.asarray(results)


def main():
    DATA.mkdir(parents=True, exist_ok=True)
    protocol = dict(
        analysis_date='2026-10-20',
        primary_question='Within Y-retaining male models, does DDX3Y RNA expression '
                         'predict DDX3X Chronos gene effect?',
        release='DepMap Public 25Q3 for both expression and gene effects.',
        expression='Use the unique IsDefaultEntryForModel=Yes row, joined by exact '
                   'ModelID. Retain supplied log2(TPM+1); no averaging, imputation '
                   'or RNA detection cut-off.',
        primary='DDX3X gene effect ~ lineage + screen efficacy + library/platform '
                'fractions + DDX3Y expression. Report partial Pearson r and '
                'expression coefficient; positive means weaker DDX3X dependency '
                'at higher DDX3Y expression.',
        uncertainty='PatientID-clustered CR1 regression SE and two-sided t(G-1) '
                    'P/95% CI. Pearson-r CI: 1000 donor bootstraps, refit '
                    'nuisance regressions, percentile interval, seed 20261020.',
        sensitivities=['Raw', 'Lineage only', 'Add DDX3X RNA expression',
                       'Add RNA strandedness and sequencing year',
                       'Restrict XY to dna_call=retained'],
        partner_expression='Test DDX3Y expression against each of the fixed '
                           '31 partner dependencies using the primary covariates; '
                           'BH over these 31 exploratory tests.',
        strata='XY expression <=XY median vs >XY median, before examining '
               'co-dependencies. Correlate adjusted DDX3X and partner residuals '
               'within each stratum; show expression-covered LOY and female '
               'models alongside. Descriptive; no significance by visual contrast.',
        moderation='Within XY: DDX3X effect ~ primary covariates + centred '
                   'DDX3Y expression + centred partner effect + their product. '
                   'Test product with donor-clustered SE, BH across 31 partners. '
                   'Negative product coefficient means a more positive partner '
                   'slope at lower DDX3Y expression. This tests slope moderation, '
                   'not equality of correlations. Sensitivity adds quadratic '
                   'DDX3Y expression, then separately DDX3X RNA expression; '
                   'each has its own 31-test BH family.',
        limitations='Post hoc same-data partners; RNA is not protein or induced '
                    'compensatory capacity. Y calls used RNA/DNA and can restrict '
                    'expression range. Stratum correlations have unequal sample '
                    'sizes and residual variance. Association is not causation.',
    )
    path = DATA / 'protocol.json'
    if path.exists():
        assert json.loads(path.read_text()) == protocol
    else:
        path.write_text(json.dumps(protocol, indent=2), encoding='utf-8')
    partners = pd.read_csv(PARTNERS).gene.tolist()
    assert len(partners) == 31 and len(set(partners)) == 31
    meta = pd.read_csv(COHORT).set_index('ModelID')
    effects = pd.read_csv(EFFECTS, index_col=0)
    np.testing.assert_allclose(effects.loc[meta.index, 'DDX3X'], meta.DDX3X)
    for gene in partners:
        meta[gene] = effects.loc[meta.index, gene]
    assert meta.index.is_unique
    columns = ['SequencingID', 'ModelID', 'IsDefaultEntryForModel',
               'ModelConditionID', 'IsDefaultEntryForMC',
               'DDX3Y (8653)', 'DDX3X (1654)']
    expression = pd.read_csv(EXPRESSION, usecols=columns)
    default = expression.loc[expression.IsDefaultEntryForModel.eq('Yes')].copy()
    assert not default.ModelID.duplicated().any()
    omics = pd.read_csv(PROFILES)
    assert omics.SequencingID.is_unique
    default = default.merge(omics[['SequencingID', 'ModelID', 'DataType', 'Stranded',
                                  'SequencingDate', 'ProfileID', 'IsDefaultEntryForModel']],
                            on=['SequencingID', 'ModelID', 'IsDefaultEntryForModel'],
                            how='left', validate='one_to_one')
    assert default.DataType.eq('rna').all()
    default = default.rename(columns={'DDX3Y (8653)': 'DDX3Y_expression',
                                      'DDX3X (1654)': 'DDX3X_expression'})
    meta = meta.join(default.set_index('ModelID'), validate='one_to_one')
    meta['expression_available'] = meta[['DDX3X_expression', 'DDX3Y_expression']].notna().all(axis=1)
    meta['rna_stranded'] = meta.Stranded.map({True: 1., False: 0.})
    meta['rna_year'] = pd.to_datetime(meta.SequencingDate).dt.year
    meta.to_csv(DATA / 'all_models_expression_join.csv')
    coverage = meta.groupby('group').agg(
        eligible_models=('PatientID', 'size'),
        expression_covered=('expression_available', 'sum'))
    coverage['missing_expression'] = coverage.eligible_models - coverage.expression_covered
    coverage.to_csv(DATA / 'coverage.csv')
    frame = meta.loc[meta.expression_available].copy()
    xy = frame.loc[frame.group.eq('XY')].copy()
    assert xy[['rna_stranded', 'rna_year']].notna().all().all()
    median = float(xy.DDX3Y_expression.median())
    frame['display_group'] = frame.group
    frame.loc[xy.index, 'display_group'] = np.where(
        xy.DDX3Y_expression.le(median), 'XY_low_DDX3Y', 'XY_high_DDX3Y')
    frame.to_csv(DATA / 'analysis_model_values.csv')
    frame.groupby('display_group')[['DDX3Y_expression', 'DDX3X']].describe().to_csv(
        DATA / 'stratum_descriptives.csv')
    sensitivity_rows, residual_rows, validation, bootstrap_rows = [], [], [], []
    for analysis in ['raw', 'lineage', 'primary', 'add_DDX3X_expression',
                     'add_RNA_covariates', 'DNA_retained_only']:
        sub = xy.loc[xy.dna_call.eq('retained')].copy() if analysis == 'DNA_retained_only' else xy
        method = analysis if analysis in ['raw', 'lineage'] else 'lineage_efficacy'
        z = nuisance_design(sub, method).to_numpy(float)
        if analysis == 'add_DDX3X_expression':
            z = np.column_stack([z, sub.DDX3X_expression])
        if analysis == 'add_RNA_covariates':
            z = np.column_stack([z, sub.rna_stranded, sub.rna_year - sub.rna_year.mean()])
        fitted, h, ry = test_term(sub.DDX3X, sub.DDX3Y_expression, z, sub.PatientID)
        row = {k: v.item() if isinstance(v, np.ndarray) else v for k, v in fitted.items()}
        row['analysis'] = analysis
        row['raw_spearman_r'] = stats.spearmanr(sub.DDX3Y_expression, sub.DDX3X).statistic
        b, se = independent_check(sub.DDX3X.to_numpy(float), sub.DDX3Y_expression.to_numpy(float),
                                  z, sub.PatientID.to_numpy())
        validation.append(dict(analysis=analysis, gene='DDX3X',
                               beta_error=abs(b - row['effect']),
                               se_error=abs(se - row['se_cluster'])))
        if analysis in ['raw', 'primary']:
            boots = bootstrap_r(sub.DDX3X, sub.DDX3Y_expression, z, sub.PatientID)
            row['r_ci_low'], row['r_ci_high'] = np.quantile(boots, [.025, .975])
            bootstrap_rows.extend(dict(analysis=analysis, iteration=i, pearson_r=r)
                                  for i, r in enumerate(boots))
        sensitivity_rows.append(row)
        residual_rows.append(pd.DataFrame(dict(
            ModelID=sub.index, analysis=analysis, expression_residual=h,
            dependency_residual=ry[:, 0])))
        print('Expression test', analysis, 'n=', len(sub),
              'r=', row['partial_r'], 'p=', row['p_value'], flush=True)
    sensitivity = pd.DataFrame(sensitivity_rows)
    sensitivity.to_csv(DATA / 'ddx3x_expression_associations.csv', index=False)
    pd.concat(residual_rows).to_csv(DATA / 'expression_dependency_residuals.csv', index=False)
    pd.DataFrame(bootstrap_rows).to_csv(DATA / 'correlation_bootstrap.csv.gz', index=False,
                                      compression={'method': 'gzip', 'mtime': 0})
    z = nuisance_design(xy, 'lineage_efficacy').to_numpy(float)
    fits, _, _ = test_term(xy[partners], xy.DDX3Y_expression, z, xy.PatientID)
    partner_expression = pd.DataFrame(dict(gene=partners, **fits))
    partner_expression['fdr_bh_31'] = stats.false_discovery_control(partner_expression.p_value)
    partner_expression.to_csv(DATA / 'partner_expression_associations.csv', index=False)
    correlations = []
    for name, sub in frame.groupby('display_group', sort=False):
        vals = sub[['DDX3X'] + partners].to_numpy(float)
        z = nuisance_design(sub, 'lineage_efficacy').to_numpy(float)
        residual = vals - z @ np.linalg.lstsq(z, vals, rcond=None)[0]
        corr = np.corrcoef(residual, rowvar=False)[0, 1:]
        for gene, r in zip(partners, corr):
            correlations.append(dict(stratum=name, gene=gene, pearson_r=r,
                                     n_models=len(sub), n_patients=sub.PatientID.nunique()))
    pd.DataFrame(correlations).to_csv(DATA / 'stratum_correlations.csv', index=False)
    interactions = []
    expression_c = xy.DDX3Y_expression.to_numpy(float) - median
    for method in ['primary', 'quadratic_expression', 'add_DDX3X_expression']:
        for gene in partners:
            partner_c = xy[gene].to_numpy(float) - xy[gene].mean()
            z = np.column_stack([nuisance_design(xy, 'lineage_efficacy').to_numpy(float),
                                 expression_c, partner_c])
            if method == 'quadratic_expression':
                z = np.column_stack([z, expression_c**2])
            elif method == 'add_DDX3X_expression':
                z = np.column_stack([z, xy.DDX3X_expression])
            term = expression_c * partner_c
            fitted, _, _ = test_term(xy.DDX3X, term, z, xy.PatientID)
            row = {k: v.item() if isinstance(v, np.ndarray) else v for k, v in fitted.items()}
            interactions.append(dict(analysis=method, gene=gene, **row))
            b, se = independent_check(xy.DDX3X.to_numpy(float), term, z, xy.PatientID.to_numpy())
            validation.append(dict(analysis='interaction_' + method, gene=gene,
                                   beta_error=abs(b - row['effect']),
                                   se_error=abs(se - row['se_cluster'])))
    interactions = pd.DataFrame(interactions)
    interactions['fdr_bh_31_within_analysis'] = interactions.groupby('analysis').p_value.transform(
        lambda p: stats.false_discovery_control(p.to_numpy()))
    interactions.to_csv(DATA / 'expression_partner_interactions.csv', index=False)
    validation = pd.DataFrame(validation)
    assert validation.beta_error.max() < 1e-10 and validation.se_error.max() < 1e-10
    validation.to_csv(DATA / 'independent_validation.csv', index=False)
    summary = dict(
        xy_models=len(xy), xy_patients=xy.PatientID.nunique(), xy_expression_median=median,
        strata=frame.display_group.value_counts().to_dict(),
        partner_expression_fdr_lt_005=int(partner_expression.fdr_bh_31.lt(.05).sum()),
        interactions={a: dict(fdr_lt_005=int(s.fdr_bh_31_within_analysis.lt(.05).sum()),
                               genes=s.loc[s.fdr_bh_31_within_analysis.lt(.05), 'gene'].tolist())
                      for a, s in interactions.groupby('analysis')},
        validation_max_beta_error=float(validation.beta_error.max()),
        validation_max_se_error=float(validation.se_error.max()))
    manifest = dict(**protocol, summary=summary, partner_order=partners,
                    packages=dict(numpy=np.__version__, pandas=pd.__version__, scipy=scipy.__version__),
                    sources={str(p): sha256(p) for p in [
                        EXPRESSION, PROFILES, COHORT, EFFECTS, PARTNERS, path, Path(__file__)]},
                    expression_metadata_note='Actual table and metadata contain both stranded '
                                             'and unstranded default RNA profiles, despite the '
                                             'generic unstranded description in the release README.',
                    expression_reference='https://forum.depmap.org/t/clarifications-in-gene-expression-level-csv-file/2603/2',
                    batch_reference='https://forum.depmap.org/t/batch-correction-best-practices-for-depmap-omics-expression-data-25q3/4580')
    (DATA / 'manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    print(json.dumps(summary, indent=2))
    print(interactions.loc[interactions.analysis.eq('primary') &
                           interactions.gene.isin(['IGBP1', 'GBF1', 'IK'])].to_string(index=False))


if __name__ == '__main__':
    main()
