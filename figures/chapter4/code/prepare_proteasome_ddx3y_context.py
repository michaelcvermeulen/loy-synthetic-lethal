"""Predefined proteasome set: DDX3Y moderation and frozen HT1080 evidence."""
from pathlib import Path
import json
import re
import numpy as np
import pandas as pd
from scipy import stats
from prepare_ddx3x_partner_adjustment import nuisance_design
from prepare_ddx3y_expression_dependency import test_term, independent_check
from prepare_xy_ranked_partners_by_y_status import sha256

THESIS = Path(__file__).resolve().parents[3]
OUT = THESIS / 'figures/chapter4'
DATA = OUT / 'source_data/proteasome_ddx3y_context'
PREV = OUT / 'source_data/ddx3y_expression_dependency'
DEP = Path('A:/DEPMAP/DepMap')
RES = THESIS.parent / 'RNAseq/analysis/HT1080/SEPT_11_2026_REANALYZE/resources'
SCREEN = THESIS.parent / 'crispr_screen/ht1080_loy_screen_2'
SET = 'GOCC_PROTEASOME_COMPLEX'
SUBSETS = {
    'GOCC_PROTEASOME_REGULATORY_PARTICLE_LID_SUBCOMPLEX': 'Lid annotation',
    'GOCC_PROTEASOME_REGULATORY_PARTICLE_BASE_SUBCOMPLEX': 'Base annotation',
    'GOCC_PROTEASOME_CORE_COMPLEX': 'Core annotation'}
GROUPS = ['XY_low_DDX3Y', 'XY_high_DDX3Y', 'LOY', 'Female']
SAMPLES = ['D0_CTRL', 'D0_LOY', 'POST_CTRL_REP1', 'POST_CTRL_REP2',
           'POST_LOY_REP1', 'POST_LOY_REP2']
ALIASES = {'ECPAS': 'KIAA0368', 'SEM1': 'SHFM1'}
FOLLOWUP = ['PSMD8', 'PSMD11', 'PSMD2', 'PSMC1', 'PSMC2', 'PSMC4']


def natural(gene):
    return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', gene)]


def main():
    DATA.mkdir(parents=True, exist_ok=True)
    membership_path = RES / 'msigdb_all_membership.csv.gz'
    members = pd.concat([
        chunk.loc[chunk.pathway.isin([SET] + list(SUBSETS))]
        for chunk in pd.read_csv(membership_path, chunksize=250000)], ignore_index=True)
    meta_sets = pd.read_csv(RES / 'msigdb_all_set_metadata.csv')
    meta_sets = meta_sets.loc[meta_sets.gs_name.isin([SET] + list(SUBSETS))]
    assert meta_sets.db_version.eq('2026.1.Hs').all()
    genes = sorted(members.loc[members.pathway.eq(SET), 'symbol'].unique(), key=natural)
    assert len(genes) == 66 and 'PSMD8' in genes
    protocol = dict(
        analysis_date='2026-10-20',
        gene_set=f'{SET}; GO:0000502; MSigDB 2026.1.Hs, 66 genes.',
        scope='Proteasome complex and associated proteins; not a claim that every '
              'annotated member is a constitutive 26S structural subunit.',
        selection='Freeze the external set before examining new dependency results. '
                  'The proteasome focus was motivated post hoc by PSMD8.',
        cohort='Reuse expression-covered 25Q3 effects, 26Q1 Chapter 3 annotations, '
               '305 XY / 216 LOY / 448 female models. Do not update membership.',
        split='Reuse XY DDX3Y RNA median 5.547863131926883 log2(TPM+1); '
              '153 lower/equal, 152 higher.',
        eligibility='At least 100 finite DDX3X/partner observations in each of the '
                    'four display strata and nonconstant residual profiles; retain '
                    'all missing or ineligible set members in the coverage table.',
        primary='Within XY, DDX3X effect ~ lineage + screen efficacy + library fractions '
                '+ centred DDX3Y RNA + standardized partner effect + their product. '
                'Partner standardized once by its sample SD in complete-case XY. '
                'Expression centred at the fixed full-XY median.',
        estimate='Product coefficient: change in DDX3X-versus-partner slope per '
                 'one log2(TPM+1) DDX3Y unit, with partner in XY SD units. '
                 'Negative means a less positive/more negative slope at higher RNA.',
        uncertainty='PatientID-clustered CR1 sandwich; two-sided t(G-1) P values '
                    'and nominal 95% intervals.',
        correction='BH across all eligible predefined-set interaction tests; '
                   'BY sensitivity for the same family. This separate exploratory '
                   'family does not replace the previous 31-partner family.',
        sensitivities=['Unadjusted interaction with both main effects',
                       'Primary plus quadratic DDX3Y RNA',
                       'Primary plus DDX3X RNA',
                       'Primary plus mean common-essential dependency excluding '
                       'DDX3X, all proteasome-set genes and all original 31 partners'],
        display='Adjusted Pearson residual correlations fit separately within each '
                'stratum. Fixed GO subset grouping, natural gene order; missing '
                'values distinct from zero. All 66 set members displayed.',
        ht1080='Reuse saved primary drugZ normZ and directional FDR without refitting. '
               'CPM uses all 77441 guides per library; guide contrast is mean '
               'endpoint log2(CPM+1) LOY minus CTRL. Retain D0 counts, CPM, '
               'zero/detection flags, and previously flagged D0 outliers. '
               'No guide exclusions, no new gene-level testing, no D0 correction '
               'of the primary drugZ results.',
        guide_display='PSMD8 and the other five proteasome genes in the original '
                      '31-partner set: PSMD11, PSMD2, PSMC1, PSMC2, PSMC4.',
        aliases=ALIASES,
        limitations='Post hoc pathway follow-up; no independent validation, causal '
                    'claim, drug prediction or claim of double-knockout suppression.')
    path = DATA / 'protocol.json'
    if path.exists():
        assert json.loads(path.read_text()) == protocol
    else:
        path.write_text(json.dumps(protocol, indent=2), encoding='utf-8')
    members.to_csv(DATA / 'gene_set_memberships.csv', index=False)
    meta_sets.to_csv(DATA / 'gene_set_metadata.csv', index=False)
    annotations = {}
    for name, label in SUBSETS.items():
        for gene in members.loc[members.pathway.eq(name), 'symbol']:
            annotations.setdefault(gene, label)
    grouped = sorted(genes, key=lambda g: (
        list(SUBSETS.values()).index(annotations[g]) if g in annotations else 3, natural(g)))
    coverage = pd.DataFrame(dict(gene=grouped))
    coverage['annotation_group'] = coverage.gene.map(annotations).fillna('Other associated')
    coverage['display_order'] = np.arange(len(coverage))
    cohort_path = PREV / 'analysis_model_values.csv'
    frame = pd.read_csv(cohort_path, index_col='ModelID')
    controls_path = DEP / 'AchillesCommonEssentialControls.csv'
    controls = pd.read_csv(controls_path).Gene.str.replace(r' \(\d+\)$', '', regex=True).tolist()
    old_partners = json.loads((PREV / 'manifest.json').read_text())['partner_order']
    wanted = set(genes + controls + ['DDX3X'])
    effect_path = DEP / 'CRISPRGeneEffect.csv'
    effects = pd.read_csv(effect_path, usecols=lambda c: c.startswith('Unnamed:')
                          or c.split(' (')[0] in wanted, index_col=0)
    effects.columns = effects.columns.str.replace(r' \(\d+\)$', '', regex=True)
    assert effects.columns.is_unique and effects.index.is_unique
    effects = effects.reindex(frame.index)
    np.testing.assert_allclose(frame.DDX3X, effects.DDX3X)
    exclusion = set(genes + old_partners + ['DDX3X'])
    clean_controls = [g for g in controls if g in effects and g not in exclusion
                      and effects[g].notna().all()]
    assert len(clean_controls) >= 400
    pd.DataFrame(dict(gene=controls,
                     used=[g in clean_controls for g in controls],
                     excluded_panel=[g in exclusion for g in controls])).to_csv(
                         DATA / 'common_essential_control_audit.csv', index=False)
    frame['essential_mean_excluding_panel'] = effects[clean_controls].mean(axis=1)
    corr_rows, interaction_rows, validations = [], [], []
    previous = pd.read_csv(PREV / 'expression_partner_interactions.csv')
    overlap_checks = []
    for gene in grouped:
        ix = coverage.gene.eq(gene)
        coverage.loc[ix, 'depmap_present'] = gene in effects
        counts = {}
        correlations = []
        for group in GROUPS:
            idx = frame.index[frame.display_group.eq(group)]
            if gene in effects:
                idx = idx[effects.loc[idx, gene].notna()]
            else:
                idx = idx[:0]
            counts[group] = len(idx)
            coverage.loc[ix, 'n_' + group] = len(idx)
            sub = frame.loc[idx]
            if len(sub) >= 100:
                vals = np.column_stack([sub.DDX3X, effects.loc[idx, gene]])
                z = nuisance_design(sub, 'lineage_efficacy').to_numpy(float)
                residual = vals - z @ np.linalg.lstsq(z, vals, rcond=None)[0]
                r = float(np.corrcoef(residual, rowvar=False)[0, 1])
            else:
                r = np.nan
            correlations.append(r)
            corr_rows.append(dict(gene=gene, stratum=group, n_models=len(idx), pearson_r=r))
        eligible = min(counts.values()) >= 100 and np.isfinite(correlations).all()
        coverage.loc[ix, 'eligible'] = eligible
        if not eligible:
            continue
        idx = frame.index[frame.group.eq('XY') & effects[gene].notna()]
        sub = frame.loc[idx]
        e = sub.DDX3Y_expression.to_numpy(float) - 5.547863131926883
        partner = effects.loc[idx, gene].to_numpy(float)
        sd = partner.std(ddof=1)
        p = (partner - partner.mean()) / sd
        for analysis in ['primary', 'raw', 'quadratic_expression', 'add_DDX3X_expression',
                         'add_essential_mean']:
            method = 'raw' if analysis == 'raw' else 'lineage_efficacy'
            z = np.column_stack([nuisance_design(sub, method).to_numpy(float), e, p])
            if analysis == 'quadratic_expression':
                z = np.column_stack([z, e**2])
            elif analysis == 'add_DDX3X_expression':
                z = np.column_stack([z, sub.DDX3X_expression])
            elif analysis == 'add_essential_mean':
                z = np.column_stack([z, sub.essential_mean_excluding_panel])
            result, _, _ = test_term(sub.DDX3X, e*p, z, sub.PatientID)
            row = {k: v.item() if isinstance(v, np.ndarray) else v for k, v in result.items()}
            row.update(gene=gene, analysis=analysis, partner_sd=sd,
                       effect_raw_partner_units=row['effect']/sd)
            interaction_rows.append(row)
            if analysis == 'primary':
                b, se = independent_check(sub.DDX3X.to_numpy(float), e*p, z,
                                          sub.PatientID.to_numpy())
                validations.append(dict(gene=gene, effect_error=abs(b-row['effect']),
                                        se_error=abs(se-row['se_cluster'])))
            match = previous.loc[previous.gene.eq(gene) & previous.analysis.eq(analysis)]
            if len(match):
                old = match.iloc[0]
                overlap_checks.append(dict(gene=gene, analysis=analysis,
                    effect_error=abs(row['effect_raw_partner_units']-old.effect),
                    p_error=abs(row['p_value']-old.p_value)))
    corr = pd.DataFrame(corr_rows)
    results = pd.DataFrame(interaction_rows)
    for method in ['bh', 'by']:
        results['fdr_' + method + '_set'] = results.groupby('analysis').p_value.transform(
            lambda p: stats.false_discovery_control(p.to_numpy(), method=method))
    results.to_csv(DATA / 'interaction_results.csv', index=False)
    corr.to_csv(DATA / 'stratum_correlations.csv', index=False)
    validation = pd.DataFrame(validations)
    checks = pd.DataFrame(overlap_checks)
    assert validation.effect_error.max() < 1e-10 and validation.se_error.max() < 1e-10
    assert checks.effect_error.max() < 1e-10 and checks.p_error.max() < 1e-10
    validation.to_csv(DATA / 'independent_fit_checks.csv', index=False)
    checks.to_csv(DATA / 'previous_results_reproduction.csv', index=False)
    print('DepMap:', len(results.query("analysis == 'primary'")), 'genes tested', flush=True)
    # Frozen screen outputs and raw guide representation.
    source_manifest = json.loads((OUT / 'source_data/ht1080_candidate_discovery/source_manifest.json').read_text())
    counts_path = Path(source_manifest['counts']['path'])
    drugz_path = Path(source_manifest['drugz']['path'])
    library_path = Path(source_manifest['library']['path'])
    for name, source in [('counts', counts_path), ('drugz', drugz_path), ('library', library_path)]:
        assert sha256(source) == source_manifest[name]['sha256']
    counts = pd.read_csv(counts_path, sep='\t').set_index('sgRNA')
    assert len(counts) == 77441 and counts.index.is_unique
    cpm = counts[SAMPLES].div(counts[SAMPLES].sum(), axis=1) * 1e6
    log = np.log2(cpm + 1)
    previous_log = pd.read_csv(OUT / 'source_data/ht1080_figure1_qc/guide_log2_cpm.csv').set_index('sgRNA')
    np.testing.assert_allclose(log, previous_log.loc[log.index, SAMPLES], atol=1e-12)
    dz = pd.read_csv(drugz_path, sep='\t').set_index('GENE')
    library = pd.read_csv(library_path, sep='\t')
    library_pairs = set(zip(library.gene, library.Seq))
    outlier_path = OUT / 'source_data/ht1080_figure1_qc/d0_outlier_guides.csv'
    outlier_ids = set(pd.read_csv(outlier_path).sgRNA)
    guide_rows, abundance_rows, screen_rows = [], [], []
    for gene in grouped:
        symbol = ALIASES.get(gene, gene)
        ids = sorted(counts.index[counts.Gene.eq(symbol)])
        record = dict(gene=gene, screen_symbol=symbol, n_guides=len(ids),
                      drugz_present=symbol in dz.index)
        if symbol in dz.index:
            z = dz.loc[symbol]
            record.update(normZ=z.normZ, fdr_synth=z.fdr_synth, fdr_supp=z.fdr_supp,
                          directional_fdr=z.fdr_synth if z.normZ < 0 else z.fdr_supp,
                          numObs=z.numObs)
        for number, guide in enumerate(ids, 1):
            assert (symbol, guide.rsplit('_', 1)[-1]) in library_pairs
            row = dict(gene=gene, screen_symbol=symbol, sgrna=guide, guide_number=number,
                       sequence=guide.rsplit('_', 1)[-1],
                       all_samples_zero=bool(counts.loc[guide, SAMPLES].sum() == 0),
                       d0_outlier_previously_flagged=guide in outlier_ids)
            for sample in SAMPLES:
                row[sample + '_reads'] = int(counts.loc[guide, sample])
                row[sample + '_cpm'] = float(cpm.loc[guide, sample])
                abundance_rows.append(dict(
                    gene=gene, sgrna=guide, guide_number=number, sample=sample,
                    genotype='CTRL' if 'CTRL' in sample else 'LOY',
                    stage='D0' if sample.startswith('D0') else 'Endpoint',
                    raw_count=int(counts.loc[guide, sample]), cpm=float(cpm.loc[guide, sample]),
                    log2_cpm_plus_1=float(log.loc[guide, sample]),
                    d0_outlier_previously_flagged=guide in outlier_ids))
            row['endpoint_loy_minus_ctrl'] = float(
                log.loc[guide, SAMPLES[4:]].mean() - log.loc[guide, SAMPLES[2:4]].mean())
            row['d0_loy_minus_ctrl'] = float(log.loc[guide, 'D0_LOY']-log.loc[guide, 'D0_CTRL'])
            row['ctrl_endpoint_minus_d0'] = float(
                log.loc[guide, SAMPLES[2:4]].mean()-log.loc[guide, 'D0_CTRL'])
            row['loy_endpoint_minus_d0'] = float(
                log.loc[guide, SAMPLES[4:]].mean()-log.loc[guide, 'D0_LOY'])
            row['endpoint_detected_samples'] = int(counts.loc[guide, SAMPLES[2:]].gt(0).sum())
            guide_rows.append(row)
        screen_rows.append(record)
    guide_table = pd.DataFrame(guide_rows)
    screen_table = pd.DataFrame(screen_rows)
    summaries = guide_table.groupby('gene').agg(
        zero_guides=('all_samples_zero', 'sum'),
        flagged_d0_guides=('d0_outlier_previously_flagged', 'sum'),
        median_guide_endpoint_contrast=('endpoint_loy_minus_ctrl', 'median'),
        median_guide_ctrl_dropout=('ctrl_endpoint_minus_d0', 'median'),
        median_guide_loy_dropout=('loy_endpoint_minus_d0', 'median'),
        min_d0_ctrl_cpm=('D0_CTRL_cpm', 'min'), min_d0_loy_cpm=('D0_LOY_cpm', 'min'))
    screen_table = screen_table.join(summaries, on='gene')
    screen_table.to_csv(DATA / 'ht1080_gene_evidence.csv', index=False)
    guide_table.to_csv(DATA / 'ht1080_guide_evidence.csv', index=False)
    pd.DataFrame(abundance_rows).to_csv(DATA / 'ht1080_guide_abundances.csv', index=False)
    coverage['screen_symbol'] = coverage.gene.map(lambda g: ALIASES.get(g, g))
    coverage['ht1080_scored'] = coverage.screen_symbol.isin(dz.index)
    coverage.to_csv(DATA / 'gene_set_coverage.csv', index=False)
    joined = coverage.merge(screen_table.drop(columns=['screen_symbol']), on='gene', validate='one_to_one')
    matrix = corr.pivot(index='gene', columns='stratum', values='pearson_r')
    joined = joined.join(matrix, on='gene').merge(
        results.loc[results.analysis.eq('primary')], on='gene', how='left', validate='one_to_one')
    joined.to_csv(DATA / 'combined_evidence.csv', index=False)
    shared = joined.loc[joined.eligible.astype(bool)]
    summary = dict(
        annotated_genes=len(genes), tested_genes=len(shared),
        missing_depmap=coverage.loc[~coverage.depmap_present.astype(bool), 'gene'].tolist(),
        ht1080_scored=int(coverage.ht1080_scored.sum()),
        primary_interaction_hits=results.loc[results.analysis.eq('primary') &
                                             results.fdr_bh_set.lt(.05), 'gene'].tolist(),
        sensitivity_hits={a: rows.loc[rows.fdr_bh_set.lt(.05), 'gene'].tolist()
                          for a, rows in results.groupby('analysis')},
        lower_positive_higher_negative=int(((shared.XY_low_DDX3Y > 0) &
                                           (shared.XY_high_DDX3Y < 0)).sum()),
        high_negative_loy_female_positive=int(((shared.XY_high_DDX3Y < 0) &
                                              (shared.LOY > 0) & (shared.Female > 0)).sum()),
        ht1080_depleted_025=screen_table.loc[screen_table.fdr_synth.lt(.25), 'gene'].tolist(),
        ht1080_enriched_025=screen_table.loc[screen_table.fdr_supp.lt(.25), 'gene'].tolist(),
        clean_essential_controls=len(clean_controls),
        validation_max_effect_error=float(validation.effect_error.max()),
        validation_max_se_error=float(validation.se_error.max()))
    manifest = dict(**protocol, summary=summary, source_hashes={
        str(p): sha256(p) for p in [membership_path, RES / 'msigdb_all_set_metadata.csv',
                                    cohort_path, effect_path, controls_path, counts_path,
                                    drugz_path, library_path, outlier_path, path, Path(__file__)]},
        gene_set_reference='https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/GOCC_PROTEASOME_COMPLEX.html',
        alias_references=['https://www.ncbi.nlm.nih.gov/gene/7979',
                          'https://www.ncbi.nlm.nih.gov/gene/23392'])
    (DATA / 'manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    print(json.dumps(summary, indent=2))
    fields = ['gene', 'XY_low_DDX3Y', 'XY_high_DDX3Y', 'LOY', 'Female',
              'effect', 'p_value', 'fdr_bh_set', 'normZ', 'directional_fdr']
    print(joined.sort_values('p_value')[fields].head(12).to_string(index=False))
    print('HT1080 follow-up genes')
    print(screen_table.loc[screen_table.gene.isin(FOLLOWUP)].to_string(index=False))


if __name__ == '__main__':
    main()
