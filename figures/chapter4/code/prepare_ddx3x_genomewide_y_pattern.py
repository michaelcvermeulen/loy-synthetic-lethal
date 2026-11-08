"""Exploratory genome-wide DDX3X partners shared by LOY and female models.

Freeze effect-size criteria before evaluating the genome-wide results. These
criteria describe a pattern and are not significance or differential tests.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from prepare_ddx3x_partner_adjustment import nuisance_design
from prepare_xy_ranked_partners_by_y_status import sha256

THESIS = Path(__file__).resolve().parents[3]
OUT = THESIS/'figures/chapter4'
PREVIOUS = OUT/'source_data/ddx3x_partner_adjustment'
DATA = OUT/'source_data/ddx3x_genomewide_y_pattern'
SOURCE = Path('A:/DEPMAP/DepMap/CRISPRGeneEffect.csv')
GROUPS = ['XY','LOY','Female']
METHODS = ['raw','lineage','lineage_efficacy']
REFERENCE_PARTNERS = ['IGBP1','GBF1','IK']
MIN_MODELS = 100
MIN_POSITIVE = .25
MAX_ABS_XY = .15
MIN_DELTA = .20


def correlations_for_group(ge,meta):
    """Batch genes by identical complete-case masks, fitting pairwise designs."""
    anchor = ge.DDX3X.to_numpy(float)
    assert np.isfinite(anchor).all()
    y = ge.to_numpy(float)
    finite = np.isfinite(y)
    packed = np.packbits(finite,axis=0).T
    _,patterns = np.unique(packed,axis=0,return_inverse=True)
    arrays = {a:np.full(y.shape[1],np.nan) for a in METHODS}
    counts = finite.sum(axis=0)
    donor_counts = np.zeros(y.shape[1],int)
    rank = np.zeros(y.shape[1],int)
    max_error = 0.
    for pattern in np.unique(patterns):
        columns = np.flatnonzero(patterns==pattern)
        mask = finite[:,columns[0]]
        assert (finite[:,columns] == mask[:,None]).all()
        frame = meta.iloc[np.flatnonzero(mask)]
        donor_counts[columns] = frame.PatientID.nunique()
        if len(frame) < MIN_MODELS:
            continue
        pair = np.column_stack([anchor[mask],y[np.ix_(mask,columns)]])
        for method in METHODS:
            x = nuisance_design(frame,method).to_numpy(float)
            beta,_,q,_ = np.linalg.lstsq(x,pair,rcond=None)
            residual = pair-x@beta
            if method=='lineage_efficacy':
                rank[columns] = q
            assert len(frame)-q > 3
            max_error = max(max_error,float(np.max(np.abs(x.T@residual))))
            num = residual[:,0]@residual[:,1:]
            den = np.sqrt(np.sum(residual[:,0]**2)*np.sum(residual[:,1:]**2,axis=0))
            r = np.divide(num,den,out=np.full(len(columns),np.nan),where=den>1e-10)
            arrays[method][columns] = np.clip(r,-1,1)
    assert max_error < 1e-7
    return arrays,counts,donor_counts,rank,int(len(np.unique(patterns))),max_error


def select_pattern(frame,min_positive=MIN_POSITIVE,max_xy=MAX_ABS_XY,min_delta=MIN_DELTA):
    return (frame.eligible & frame.loy_r.ge(min_positive) & frame.female_r.ge(min_positive) &
            frame.xy_r.abs().le(max_xy) & frame.loy_minus_xy.ge(min_delta) &
            frame.female_minus_xy.ge(min_delta))


def main():
    DATA.mkdir(parents=True,exist_ok=True)
    previous = json.loads((PREVIOUS/'manifest.json').read_text())
    digest = sha256(SOURCE)
    assert digest == previous['sources'][str(SOURCE)]
    meta = pd.read_csv(PREVIOUS/'model_values.csv').set_index('ModelID')
    ge = pd.read_csv(SOURCE,index_col=0)
    original_names = ge.columns.tolist()
    ge.columns = ge.columns.str.replace(r' \(\d+\)$','',regex=True)
    assert ge.columns.is_unique and ge.index.is_unique
    result = pd.DataFrame({'gene':ge.columns,'source_feature':original_names})
    patterns,errors = {},{}
    for group in GROUPS:
        frame = meta.loc[meta.group.eq(group)]
        values = ge.loc[frame.index]
        arrays,n,nd,q,npats,err = correlations_for_group(values,frame)
        prefix = group.lower()
        result[f'{prefix}_r'] = arrays['lineage_efficacy']
        result[f'{prefix}_raw_r'] = arrays['raw']
        result[f'{prefix}_lineage_r'] = arrays['lineage']
        result[f'{prefix}_n_models'] = n
        result[f'{prefix}_n_patients'] = nd
        result[f'{prefix}_design_rank'] = q
        patterns[group],errors[group] = npats,err
        print(group,'complete-case patterns',npats,'estimable genes',int(np.isfinite(arrays['lineage_efficacy']).sum()),flush=True)
    result['eligible'] = (result.gene.ne('DDX3X') &
        result[[f'{g.lower()}_n_models' for g in GROUPS]].ge(MIN_MODELS).all(axis=1) &
        result[['xy_r','loy_r','female_r']].notna().all(axis=1))
    result['loy_minus_xy'] = result.loy_r-result.xy_r
    result['female_minus_xy'] = result.female_r-result.xy_r
    result['min_y_absent_r'] = result[['loy_r','female_r']].min(axis=1)
    result['min_delta'] = result[['loy_minus_xy','female_minus_xy']].min(axis=1)
    result['selected'] = select_pattern(result)
    result['reference_partner'] = result.gene.isin(REFERENCE_PARTNERS)
    result.to_csv(DATA/'all_gene_correlations.csv',index=False)
    selected = result.loc[result.selected].sort_values(
        ['min_y_absent_r','min_delta','gene'],ascending=[False,False,True]).copy()
    selected.insert(0,'display_rank',np.arange(1,len(selected)+1))
    selected.to_csv(DATA/'selected_partners.csv',index=False)
    sensitivity = []
    for positive in [.20,.25,.30]:
        for abs_xy in [.10,.15,.20]:
            for delta in [.15,.20,.25]:
                keep = select_pattern(result,positive,abs_xy,delta)
                sensitivity.append(dict(min_positive_r=positive,max_abs_xy_r=abs_xy,min_delta=delta,
                                        genes=int(keep.sum()),other_than_original_three=int((keep&~result.reference_partner).sum()),
                                        primary=(positive==MIN_POSITIVE and abs_xy==MAX_ABS_XY and delta==MIN_DELTA)))
    pd.DataFrame(sensitivity).to_csv(DATA/'threshold_sensitivity_counts.csv',index=False)
    old = pd.read_csv(PREVIOUS/'correlations.csv')
    check_errors = []
    for row in old.loc[old.analysis.isin(METHODS)].itertuples():
        suffix = {'raw':'raw_r','lineage':'lineage_r','lineage_efficacy':'r'}[row.analysis]
        new = result.loc[result.gene.eq(row.partner)].iloc[0]
        assert new[f'{row.group.lower()}_n_models']==row.n_models
        check_errors.append(abs(new[f'{row.group.lower()}_{suffix}']-row.pearson_r))
    assert max(check_errors) < 1e-10
    # Independently refit each selected gene, including its pair-specific anchor fit.
    selected_errors = []
    for gene in selected.gene:
        for group in GROUPS:
            frame = meta.loc[meta.group.eq(group)]
            pair = ge.loc[frame.index,['DDX3X',gene]].dropna()
            frame = frame.loc[pair.index]
            x = nuisance_design(frame,'lineage_efficacy').to_numpy(float)
            residual = pair.to_numpy()-x@np.linalg.lstsq(x,pair.to_numpy(),rcond=None)[0]
            independent = np.corrcoef(residual.T)[0,1]
            saved = selected.loc[selected.gene.eq(gene),f'{group.lower()}_r'].iloc[0]
            selected_errors.append(abs(independent-saved))
    assert not selected_errors or max(selected_errors)<1e-10
    manifest = dict(
        seed='DDX3X',dependency_release=previous['dependency_release'],annotation_release=previous['annotation_release'],
        cohort='Exact same XY, LOY and eligible female models as the preceding adjustment figure.',groups=previous['groups'],
        total_gene_profiles=len(result),eligible_nonself_genes=int(result.eligible.sum()),
        min_paired_models_each_group=MIN_MODELS,
        adjustment=previous['primary_adjustment'],
        missingness='Genes grouped by identical complete-case masks for computation; both DDX3X and partner nuisance fits use the same pair-complete models. No imputation.',
        selection={'loy_and_female_min_r':MIN_POSITIVE,'xy_max_absolute_r':MAX_ABS_XY,'min_increase_vs_xy_in_both':MIN_DELTA,
                   'type':'Exploratory effect-size pattern, not a statistical-significance or differential-correlation test.'},
        selected_genes=len(selected),other_than_original_three=int((selected.selected&~selected.reference_partner).sum()),
        original_three_retained=selected.loc[selected.reference_partner,'gene'].tolist(),
        order='Descending minimum of LOY and female adjusted r; tie break minimum increase over XY, then gene symbol.',
        display_labels={'XY':'Y-retaining males','X0':'Primary Y-negative/LOY males','XX':'Eligible recorded-female models'},
        karyotype_note='XY, X0 and XX are requested display shorthands for the analysis groups, not verified chromosome counts.',
        inference='No P values, FDR calls or claim of biological absence in XY. Selection on these same data inflates contrasts. Number depends on the stated descriptive cutoffs; nearby thresholds are saved.',
        validation={'previous_coefficients_reproduced':len(check_errors),'max_previous_r_error':max(check_errors),
                    'independent_selected_pair_checks':len(selected_errors),'max_selected_r_error':max(selected_errors,default=0),
                    'complete_case_patterns':patterns,'max_normal_equation_errors':errors},
        sources={str(SOURCE):digest,**{str(p):sha256(p) for p in [PREVIOUS/'manifest.json',PREVIOUS/'model_values.csv',PREVIOUS/'correlations.csv']}})
    (DATA/'manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    print(json.dumps({k:manifest[k] for k in ['total_gene_profiles','eligible_nonself_genes','selected_genes','other_than_original_three','original_three_retained','validation']},indent=2))
    print(selected[['gene','xy_r','loy_r','female_r','loy_minus_xy','female_minus_xy']].round(3).to_string(index=False))


if __name__=='__main__':
    main()
