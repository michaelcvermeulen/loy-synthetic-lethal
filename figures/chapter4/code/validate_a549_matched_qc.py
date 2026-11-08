"""Independent checks for the displayed A549 QC data and integration links."""
from pathlib import Path
import hashlib
import json
import re
import sys
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT.parent/'.figure_python_deps'))
import numpy as np
import pandas as pd
DATA=ROOT/'figures/chapter4/source_data/a549_matched_qc'
INITIAL=ROOT/'figures/chapter4/source_data/a549_initial_qc'

def main():
    manifest=json.loads((DATA/'analysis_manifest.json').read_text())
    summary=json.loads((DATA/'summary.json').read_text())
    for key,item in manifest['inputs'].items():
        assert hashlib.sha256(Path(item['path']).read_bytes()).hexdigest()==item['sha256'],key
    for name in ['baseline_gene_qc.csv','sample_qc.csv','sample_metadata.csv']:
        assert (DATA/name).read_bytes()==(INITIAL/name).read_bytes(),name
    retained=pd.read_csv(manifest['inputs']['retained_counts']['path'],sep='\t').set_index('sgRNA')
    samples=summary['samples']
    assert len(samples)==7 and not any('_P3_' in s for s in samples)
    log=pd.read_csv(DATA/'guide_log2_cpm.csv',index_col=0)
    assert log.columns.tolist()==samples
    expected=np.log2(retained[samples].div(retained[samples].sum(),axis=1)*1e6+1)
    np.testing.assert_allclose(log,expected,atol=1e-12)
    corr=pd.read_csv(DATA/'spearman_guide_counts.csv',index_col=0)
    pd.testing.assert_frame_equal(corr,pd.read_csv(INITIAL/'spearman_all_samples.csv',index_col=0))
    genes=pd.read_csv(DATA/'baseline_gene_qc.csv')
    for g in ['CTRL','LOY']:
        d=genes.loc[genes.genotype.eq(g)&genes.category.isin(['Essential','Non-essential'])]
        positive=d.category.eq('Essential').to_numpy()
        ranks=pd.Series(-d.lfc.to_numpy()).rank(method='average').to_numpy()
        n,m=positive.sum(),(~positive).sum()
        auc=(ranks[positive].sum()-n*(n+1)/2)/(n*m)
        assert abs(auc-summary['roc'][g]['auc'])<1e-12
    guides=pd.read_csv(manifest['inputs']['mageck_guides']['path'],sep='\t').set_index('sgrna')
    plotted=pd.read_csv(DATA/'endpoint_reference_density_guides.csv').set_index('sgrna')
    np.testing.assert_array_equal(plotted.LFC,guides.loc[plotted.index,'LFC'])
    total=retained[[s for s in samples if '_P4_' in s]].sum(axis=1)
    np.testing.assert_array_equal(plotted.included,total.loc[plotted.index].gt(0))
    # Saved means have five significant digits; three nearly-zero LFCs round to tied means.
    sign_mismatch=(guides.LFC>0)!=(guides.treat_mean>guides.control_mean)
    rounding_ties=guides.treat_mean.eq(guides.control_mean)&guides.LFC.abs().lt(2e-5)
    assert (~sign_mismatch | rounding_ties).all()
    curves=pd.read_csv(DATA/'endpoint_reference_density_curves.csv')
    for cat in ['Non-targeting','Essential','Non-essential']:
        assert abs(np.trapezoid(curves[cat],curves.LFC)-1)<.0001
        d=plotted.loc[plotted.category.eq(cat)&plotted.included,'LFC']
        assert len(d)==summary['density']['groups'][cat]['included_guides']
        assert curves.LFC.min()<d.min() and curves.LFC.max()>d.max()
    result={
        'all_original_and_derived_inputs_unchanged':True,
        'seven_d0_p4_libraries_only':True,
        'guide_cpm_and_correlations_verified':True,
        'baseline_gene_scores_preserved_and_roc_rank_verified':True,
        'endpoint_lfc_values_direction_and_read_eligibility_verified':True,
        'rounded_mean_ties_with_abs_lfc_below_2e_5':int(sign_mismatch.sum()),
        'density_integrals_and_full_observed_ranges_verified':True,
        'figure_export_checks':json.loads((DATA/'figure_validation.json').read_text()),
    }
    for rel in ['Chapter4_figure_captions.md','draftv1/chapter4/chapter4_results.md','research/chapter4_a549_processing_20261019/README.md']:
        path=ROOT/rel
        text=path.read_text(encoding='utf-8')
        for dest in re.findall(r'\]\(([^)]+)\)',text):
            if 'a549' not in dest.lower() or dest.startswith('http'): continue
            dest=dest.split('#')[0]
            assert (Path(dest) if re.match(r'^[A-Z]:/',dest) else path.parent/dest).exists(),(rel,dest)
    method=(ROOT/'draftv1/chapter4/chapter4_methods.md').read_text(encoding='utf-8')
    a549=method.split('### 4.2.3 CRISPR screening in A549')[1].split('### 4.2.4')[0]
    assert 'no variance scaling' not in a549 and '307 constant' in a549
    assert 'PCR jackpotting' in a549
    result['active_a549_links_and_methods_verified']=True
    (DATA/'independent_validation.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
