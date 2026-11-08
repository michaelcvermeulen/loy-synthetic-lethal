"""Exploratory DDX3X partner correlations, with female and nuisance checks.

Frozen 25Q3 effects; current Chapter 3 eligibility/Y annotations by ModelID.
No partner re-selection, screen-hit changes or differential-correlation tests.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from prepare_xy_ranked_partners_by_y_status import sha256

THESIS = Path(__file__).resolve().parents[3]
OUT = THESIS/'figures/chapter4'
DATA = OUT/'source_data/ddx3x_partner_adjustment'
DEP = Path('A:/DEPMAP/DepMap')
PREVIOUS = OUT/'source_data/ddx3x_partner_y_status_scatter'
COHORT = OUT/'source_data/xy_paralog_ranked_partners_y_status/model_cohort_and_exclusions.csv'
CALLS = THESIS/'research/depmap_loy_reanalysis_v1/loy_reanalysis_v1_export_2026-10-07/derived/loy_reanalysis_v1/y_status_calls.csv'
GENES = ['DDX3X', 'IGBP1', 'GBF1', 'IK']
GROUPS = ['XY', 'LOY', 'Female']
ANALYSES = ['raw', 'lineage', 'lineage_efficacy', 'lineage_essential_mean']


def nuisance_design(frame, analysis):
    if analysis == 'raw':
        design = pd.DataFrame({'intercept': np.ones(len(frame))}, index=frame.index)
    else:
        design = pd.get_dummies(frame.OncotreeLineage, prefix='lineage', dtype=float)
        if analysis == 'lineage_efficacy':
            # Composition captures systematic shifts between library/platforms.
            # One fraction is redundant with the full lineage intercept design.
            cols = ['screen_efficacy'] + [c for c in frame if c.startswith('library_fraction_')][1:]
        elif analysis == 'lineage_essential_mean':
            cols = ['essential_mean']
        else:
            cols = []
        for col in cols:
            design[col] = frame[col].to_numpy(float)
    return design


def fit_residuals(frame, analysis):
    """Fit both genes on the same observations and nuisance design, per group."""
    y = frame[GENES].to_numpy(float)
    design = nuisance_design(frame, analysis)
    x = design.to_numpy(float)
    beta, _, rank, _ = np.linalg.lstsq(x, y, rcond=None)
    residual = y-x@beta
    error = float(np.max(np.abs(x.T@residual)))
    assert error < 1e-8, (analysis, error)
    assert len(frame)-rank > 3
    coefs = pd.DataFrame(beta, index=design.columns, columns=GENES)
    return residual, int(rank), coefs, error


def bootstrap_correlations(selected, n_boot=1000):
    """Resample donors within each group and refit nuisance models each time."""
    rng = np.random.default_rng(20261019)
    records = []
    for group in GROUPS:
        frame = selected.loc[selected.group.eq(group)]
        donor_codes, donors = pd.factorize(frame.PatientID, sort=True)
        y = frame[GENES].to_numpy(float)
        designs = {a:nuisance_design(frame,a).to_numpy(float) for a in ANALYSES}
        for iteration in range(n_boot):
            donor_counts = rng.multinomial(len(donors), np.full(len(donors),1/len(donors)))
            root_weights = np.sqrt(donor_counts[donor_codes])[:,None]
            yw = y*root_weights
            for analysis,x in designs.items():
                xw = x*root_weights
                beta = np.linalg.lstsq(xw,yw,rcond=None)[0]
                rw = yw-xw@beta
                cross = rw.T@rw
                rs = cross[0,1:]/np.sqrt(cross[0,0]*np.diag(cross)[1:])
                assert np.isfinite(rs).all()
                for partner,r in zip(GENES[1:],rs):
                    records.append(dict(group=group,analysis=analysis,partner=partner,
                                        iteration=iteration,pearson_r=float(r)))
    return pd.DataFrame(records)


def main():
    DATA.mkdir(parents=True, exist_ok=True)
    source = DEP/'CRISPRGeneEffect.csv'
    old_manifest = json.loads((PREVIOUS/'manifest.json').read_text())
    digest = sha256(source)
    assert digest == old_manifest['source_sha256'][str(source)]
    cohort = pd.read_csv(COHORT).set_index('ModelID')
    calls = pd.read_csv(CALLS).set_index('ModelID')
    assert cohort.index.is_unique and calls.index.is_unique
    meta = calls.reindex(cohort.index)
    female = (cohort.Sex.eq('Female') & meta.Sex.eq('Female') &
              meta.OncotreePrimaryDisease.notna() &
              meta.OncotreePrimaryDisease.ne('Non-Cancerous') &
              meta.EngineeredModel.fillna('').eq(''))
    cohort['included_previous_male'] = cohort.included
    cohort['EngineeredModel_current'] = meta.EngineeredModel
    cohort['disease_current'] = meta.OncotreePrimaryDisease
    cohort.loc[cohort.Sex.eq('Female'), 'PatientID'] = meta.loc[cohort.Sex.eq('Female'), 'PatientID']
    cohort.loc[female, 'group'] = 'Female'
    cohort.loc[female, 'included'] = True
    cohort.loc[female, 'exclusion_reason'] = ''
    excluded_female = cohort.Sex.eq('Female') & ~female
    cohort.loc[excluded_female, 'exclusion_reason'] = np.where(
        cohort.loc[excluded_female, 'disease_current'].eq('Non-Cancerous'),
        'Female: non-cancerous', 'Female: engineered or drug-adapted model')
    cohort.to_csv(DATA/'cohort_and_exclusions.csv', index_label='ModelID')
    selected = cohort.loc[cohort.included].copy()
    assert selected.group.isin(GROUPS).all()
    assert selected.PatientID.notna().all() and selected.OncotreeLineage.notna().all()
    assert selected[selected.group.ne('Female')].index.tolist() == cohort[cohort.included_previous_male].index.tolist()

    controls = pd.read_csv(DEP/'AchillesCommonEssentialControls.csv').Gene
    control_symbols = controls.str.replace(r' \(\d+\)$', '', regex=True)
    wanted = set(control_symbols) | set(GENES)
    ge = pd.read_csv(source, usecols=lambda c: c.startswith('Unnamed:') or
                     c.rsplit(' (', 1)[0] in wanted, index_col=0)
    ge.columns = ge.columns.str.replace(r' \(\d+\)$', '', regex=True)
    assert ge.columns.is_unique and ge.index.tolist() == cohort.index.tolist()
    values = ge.reindex(selected.index)
    assert values[GENES].notna().all().all()
    control_audit = pd.DataFrame({'gene':control_symbols})
    control_audit['present'] = control_audit.gene.isin(ge.columns)
    control_audit['is_plotted_gene'] = control_audit.gene.isin(GENES)
    control_audit['n_models_measured'] = control_audit.gene.map(values.notna().sum()).fillna(0).astype(int)
    control_audit['used'] = (control_audit.present & ~control_audit.is_plotted_gene &
                             control_audit.n_models_measured.eq(len(selected)))
    control_audit.to_csv(DATA/'common_essential_control_selection.csv', index=False)
    used_controls = control_audit.loc[control_audit.used, 'gene'].tolist()
    assert len(used_controls) >= 500
    selected = selected.join(values[GENES])
    selected['essential_mean'] = values[used_controls].mean(axis=1)
    selected['essential_median'] = values[used_controls].median(axis=1)

    efficacy = pd.read_csv(DEP/'CRISPRInferredModelEfficacy.csv').rename(columns={'ModelID':'source_row_id'})
    libraries = [c for c in efficacy if c != 'source_row_id']
    assert efficacy[libraries].notna().sum(axis=1).eq(1).all()
    efficacy['ModelID'] = efficacy.source_row_id.str.extract(r'^(ACH-\d{6})(?:_|$)', expand=False)
    assert efficacy.ModelID.notna().all()
    assert set(efficacy.ModelID) == set(ge.index)
    long = efficacy.melt(id_vars=['source_row_id','ModelID'], value_vars=libraries,
                         var_name='library', value_name='efficacy').dropna(subset=['efficacy'])
    assert long.efficacy.between(0,1).all()
    long.to_csv(DATA/'efficacy_input_rows.csv', index=False)
    bylib = long.groupby(['ModelID','library']).efficacy.mean().unstack('library').reindex(columns=libraries)
    # Equal library/platform weighting, avoiding arbitrary replicate-number weights.
    model_efficacy = bylib.mean(axis=1)
    fractions = bylib.notna().div(bylib.notna().sum(axis=1), axis=0)
    fractions.columns = ['library_fraction_'+c for c in fractions.columns]
    selected['screen_efficacy'] = model_efficacy.reindex(selected.index)
    selected['n_efficacy_source_rows'] = long.groupby('ModelID').size().reindex(selected.index)
    selected['n_libraries'] = bylib.notna().sum(axis=1).reindex(selected.index)
    selected = selected.join(fractions)
    assert selected.screen_efficacy.notna().all()
    selected.to_csv(DATA/'model_values.csv', index_label='ModelID')
    selected.groupby(['group','OncotreeLineage']).size().rename('n_models').reset_index().to_csv(
        DATA/'lineage_composition.csv', index=False)

    stats, residual_rows, coefficient_rows, checks = [], [], [], []
    old = pd.read_csv(PREVIOUS/'correlations.csv').set_index(['partner','group'])
    # Within-group fits let lineage effects and sensitivity slopes differ by group.
    for group in GROUPS:
        frame = selected.loc[selected.group.eq(group)]
        for analysis in ANALYSES:
            residual, rank, coefs, error = fit_residuals(frame, analysis)
            checks.append(error)
            coefs = coefs.reset_index(names='term').melt(id_vars='term', var_name='gene', value_name='coefficient')
            coefs['group'], coefs['analysis'] = group, analysis
            coefficient_rows.append(coefs)
            residual_frame = pd.DataFrame(residual, index=frame.index, columns=GENES)
            residual_frame['group'], residual_frame['analysis'] = group, analysis
            residual_rows.append(residual_frame)
            for j,partner in enumerate(GENES[1:],1):
                r = float(np.corrcoef(residual[:,0],residual[:,j])[0,1])
                direct = float(residual[:,0]@residual[:,j]/np.sqrt(
                    (residual[:,0]@residual[:,0])*(residual[:,j]@residual[:,j])))
                assert abs(r-direct) < 1e-12
                rho = float(pd.DataFrame(residual[:,[0,j]]).corr(method='spearman').iloc[0,1])
                if analysis == 'raw' and group != 'Female':
                    assert abs(r-old.loc[(partner,group),'pearson_r']) < 1e-12
                    assert len(frame) == old.loc[(partner,group),'n_models']
                slope, intercept = np.polyfit(residual[:,0], residual[:,j], 1)
                stats.append(dict(partner=partner,group=group,analysis=analysis,
                                  n_models=len(frame),n_patients=frame.PatientID.nunique(),
                                  n_lineages=frame.OncotreeLineage.nunique(),design_rank=rank,
                                  residual_df=len(frame)-rank,pearson_r=r,residual_spearman_r=rho,
                                  residual_slope=float(slope),residual_intercept=float(intercept)))
    stats = pd.DataFrame(stats)
    assert len(stats) == 36 and stats.pearson_r.between(-1,1).all()
    boots = bootstrap_correlations(selected)
    boots.to_csv(DATA/'donor_bootstrap_correlations.csv',index=False)
    intervals = boots.groupby(['partner','group','analysis']).pearson_r.quantile([.025,.975]).unstack()
    intervals.columns = ['bootstrap_ci_low','bootstrap_ci_high']
    stats = stats.join(intervals,on=['partner','group','analysis'],validate='one_to_one')
    stats.to_csv(DATA/'correlations.csv', index=False)
    pd.concat(residual_rows).to_csv(DATA/'residuals.csv', index_label='ModelID')
    pd.concat(coefficient_rows).to_csv(DATA/'nuisance_coefficients.csv', index=False)
    selected.groupby('group')[GENES+['screen_efficacy','essential_mean','essential_median']].agg(['mean','std','min','max']).to_csv(DATA/'group_descriptives.csv')
    sources = [source,DEP/'README.txt',DEP/'CRISPRInferredModelEfficacy.csv',
               DEP/'AchillesCommonEssentialControls.csv',COHORT,CALLS,PREVIOUS/'correlations.csv']
    manifest = dict(
        dependency_release='DepMap Public 25Q3',annotation_release='Chapter 3 26Q1 primary calls/eligibility, export 2026-10-07',
        cross_release_join='Exact ModelID; no updated dependency values or Y-status reclassification.',
        groups={g:{'models':int(len(f)),'patients':int(f.PatientID.nunique()),'lineages':int(f.OncotreeLineage.nunique())}
                for g,f in selected.groupby('group')},
        female_eligibility='Recorded female in both sources, cancer and no engineered/drug-adapted annotation, using current Chapter 3 metadata. Female does not assert an XX karyotype.',
        excluded_females=cohort.loc[excluded_female,'exclusion_reason'].value_counts().to_dict(),
        male_eligibility='Exact saved 324 XY and 229 LOY cohorts, unchanged.',
        partners='IGBP1, GBF1 and IK, previously selected as the top three DDX3X partners within LOY; no reranking.',
        primary_adjustment='Within each group, separate OLS for each gene on Oncotree lineage indicators + mean Chronos inferred efficacy + library/platform composition fractions; correlate the paired residuals.',
        efficacy_aggregation='The ModelID-labelled input contains replicate-like IDs. Extract the anchored ACH-###### prefix, validate every ID against the frozen matrix, average efficacy rows within model/library-platform, then equally average measured library-platform means. Composition is equal fractions of available library-platforms; omit one redundant fraction in the design.',
        sensitivity_adjustment='Separate OLS on lineage + mean common-essential gene effect; same fixed controls in every model, excluding all four plotted genes and controls with any missing selected-model value. This proxy can include biology as well as technical sensitivity.',
        common_essential_controls_used=len(used_controls),
        control_exclusion_counts={'listed':len(control_audit),'plotted_genes':int(control_audit.is_plotted_gene.sum()),
                                  'absent':int((~control_audit.present).sum()),
                                  'measured_incompletely_nonplotted':int((control_audit.present & ~control_audit.is_plotted_gene & ~control_audit.used).sum())},
        raw_comparison='All four analyses use identical complete observations within group; no outlier exclusions or imputation.',
        intervals='95% percentile intervals from 1000 donor-cluster bootstrap resamples within each group, seed 20261019. Resample PatientIDs with replacement, retain all models of sampled donors with multiplicities, and refit every nuisance model in every resample. No partner reselection; intervals are conditional on the exploratory partner choice.',
        already_corrected='The 25Q3 released gene effects already incorporate copy-number and screen-quality corrections. These are additional robustness checks, not a reconstruction of the original pipeline.',
        interpretation='Descriptive partial correlations, not tests of differences between groups. Model-level estimates; related models retained. Selection of partners in LOY can inflate their LOY correlations. No claim of causal regulation or new LOY-selective dependencies.',
        validation={'original_six_raw_correlations_and_counts_reproduced':True,'same_models_all_analyses':True,
                    'all_efficacy_prefixes_validated':True,'max_OLS_normal_equation_error':max(checks),
                    'no_missing_covariates':True},
        sources={str(p):sha256(p) for p in sources})
    (DATA/'manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    print(stats.pivot(index=['partner','group'],columns='analysis',values='pearson_r').round(4).to_string())
    print(json.dumps({k:manifest[k] for k in ['groups','excluded_females','common_essential_controls_used','validation']},indent=2))


if __name__ == '__main__':
    main()
