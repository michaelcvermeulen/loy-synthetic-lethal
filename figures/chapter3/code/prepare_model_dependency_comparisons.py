"""Freeze descriptive RNA and source-matched DNA dependency comparisons."""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np
import pandas as pd

GENES = ['DDX3X','EIF1AX','RPS4X','ZFX']
LINES = {'ACH-000054':'HT1080','ACH-000681':'A549','ACH-000979':'DU145','ACH-000092':'H2452'}


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream,'sha256').hexdigest()


def main(args):
    root=args.data_dir
    sources={
        'gene_effect':root/'raw/depmap_public_26q1/crispr/CRISPRGeneEffect.csv',
        'rna_calls':root/'derived/rna_y_vulnerabilities_v1/final/rna_y_calls.csv',
        'dna_cohort':root/'derived/loy_dependencies_v1/primary_cohort.csv',
        'dna_audit':root/'derived/loy_dependencies_v1/cohort_audit.csv',
        'dna_statistics':root/'derived/loy_dependencies_v1/loy_dependency_results.csv',
        'dna_design':root/'derived/loy_dependencies_v1/design_primary.csv',
    }
    initial={k:(p.stat().st_size,p.stat().st_mtime_ns) for k,p in sources.items()}
    header=pd.read_csv(sources['gene_effect'],nrows=0).columns
    cols=[c for c in header if c.split(' (')[0] in GENES]
    effect=pd.read_csv(sources['gene_effect'],usecols=[header[0]]+cols).rename(columns={header[0]:'ModelID',**{c:c.split(' (')[0] for c in cols}})
    assert not effect.ModelID.duplicated().any()
    assert effect.loc[effect.ModelID.isin(LINES),GENES].notna().all().all()
    rna=pd.read_csv(sources['rna_calls'])
    dna=pd.read_csv(sources['dna_cohort'])
    audit=pd.read_csv(sources['dna_audit'],low_memory=False)
    results=pd.read_csv(sources['dna_statistics'],low_memory=False)
    records={}
    for kind in ['rna','dna']:
        if kind=='rna':
            cohort=rna.loc[rna.included_donor_crispr_integrated.eq(True) & rna.y_state.isin(['EXPRESSED','DEPLETED'])].copy()
            assert cohort.donor_key.is_unique
            assert cohort.Sex.eq('Male').all()
            cohort['LOY_status']=cohort.y_state.map({'EXPRESSED':'XY','DEPLETED':'LOY'})
            names=rna[['ModelID','CellLineName']].copy()
            group_labels={'XY':'Y-expressing','LOY':'Y-depleted'}
        else:
            cohort=dna.copy()
            assert cohort.donor_key.is_unique and cohort.Sex.eq('Male').all()
            cohort['LOY_status']=cohort.LOY.map({0.:'XY',1.:'LOY'})
            assert cohort.LOY_status.notna().all()
            names=audit[['ModelID','CellLineName']].copy()
            group_labels={'XY':'Y-retaining','LOY':'LOY'}
        subset=cohort[['ModelID','LOY_status']].copy()
        subset['in_reference']=True
        extra=pd.DataFrame({'ModelID':[x for x in LINES if x not in set(subset.ModelID)],'LOY_status':'outside_reference','in_reference':False})
        subset=pd.concat([subset,extra],ignore_index=True).merge(names,on='ModelID',validate='one_to_one')
        subset['cell_line']=subset.ModelID.map(LINES).fillna(subset.CellLineName)
        subset=subset.merge(effect,on='ModelID',validate='one_to_one')
        long=subset.melt(id_vars=['ModelID','cell_line','LOY_status','in_reference'],value_vars=GENES,var_name='gene',value_name='gene_effect')
        long['highlight']=long.ModelID.isin(LINES)
        summaries=[]
        for gene in GENES:
            sub=long.loc[long.gene.eq(gene) & long.in_reference]
            bygroup={g:sub.loc[sub.LOY_status.eq(g),'gene_effect'].dropna() for g in ['XY','LOY']}
            s={'gene':gene,'n_XY':len(bygroup['XY']),'n_LOY':len(bygroup['LOY']),
               'mean_XY':bygroup['XY'].mean(),'mean_LOY':bygroup['LOY'].mean()}
            s['mean_difference']=s['mean_LOY']-s['mean_XY']
            if kind=='rna':
                s['plot_annotation']=f"Mean difference = {s['mean_difference']:+.2f}"
                s['method']='Descriptive mean difference, Y-depleted minus Y-expressing; no inferential test'
            else:
                ref=results.loc[results.symbol.eq(gene)].iloc[0]
                assert np.isclose(s['mean_LOY'],ref.primary_mean_LOY,atol=1e-10)
                assert np.isclose(s['mean_XY'],ref.primary_mean_Y_present,atol=1e-10)
                assert s['n_XY']==ref.primary_n_Y_present and s['n_LOY']==ref.primary_n_LOY
                for c in ['primary_delta','primary_se','primary_ci_low','primary_ci_high','primary_p_value','primary_fdr']:
                    s[c]=ref[c]
                p=ref.primary_fdr
                pf=f'{p:.2g}'
                s['plot_annotation']=f"Adjusted difference = {ref.primary_delta:+.2f}\nFDR-adjusted P = {pf}"
                s['method']='Existing primary dependency regression; source-matched cohort and design'
            summaries.append(s)
        dest=args.output_dir/kind
        dest.mkdir(parents=True,exist_ok=False)
        long.to_csv(dest/'dependency_comparison.csv',index=False)
        pd.DataFrame(summaries).to_csv(dest/'dependency_statistics.csv',index=False)
        cohort.to_csv(dest/'reference_cohort.csv',index=False)
        config={
            'release_label':'DepMap Public 26Q1','group_labels':group_labels,
            'comparison_label':'RNA-based Y status; all eligible donor representatives' if kind=='rna' else 'DNA-supported primary cohort; grey models lie outside this cohort',
            'interpretation':'RNA depletion is not proof of physical LOY or a diploid XY karyotype' if kind=='rna' else 'Same 151-model primary cohort as the saved adjusted dependency analysis',
            'output_stem':'ch03_main_selected_models_paralog_dependencies' if kind=='rna' else 'ch03_supp_selected_models_dna_paralog_dependencies',
            'statistics_kind':kind,
        }
        (dest/'figure_config.json').write_text(json.dumps(config,indent=2)+'\n',encoding='utf-8')
        records[kind]={'counts':cohort.LOY_status.value_counts().to_dict(),'statistics':summaries}
    evidence=rna.loc[rna.ModelID.isin(LINES),['ModelID','CellLineName','y_state','RNA_expression_score','validation_WGS_call','validation_WGS_copy_estimate','DNA_first_y_state','included_donor_crispr_integrated']]
    evidence=evidence.merge(audit[['ModelID','review_LOY_recommended_call','review_LOY_status','review_LOY_confidence','review_LOY_reason','primary_exclusion_reason','cohort_primary']],on='ModelID',validate='one_to_one')
    evidence.to_csv(args.output_dir/'selected_model_y_status_evidence.csv',index=False)
    provenance={}
    for label,path in sources.items():
        digest=sha(path)
        assert initial[label]==(path.stat().st_size,path.stat().st_mtime_ns),f'Transfer still modifying {path}'
        provenance[label]={'path':str(path),'bytes':path.stat().st_size,'sha256':digest}
    (args.output_dir/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(records,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--data-dir',type=Path,required=True)
    p.add_argument('--output-dir',type=Path,required=True)
    main(p.parse_args())
