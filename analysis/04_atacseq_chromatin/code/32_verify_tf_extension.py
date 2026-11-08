"""Cross-check sample effects, matching, multiplicity, and finished TF artifacts."""
from pathlib import Path
import json, sys, re
import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control

ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/'motifs/TF_extension'
preliminary='--before-footprints' in sys.argv
checks=[]
samples=pd.read_csv(ROOT/'config/samples.tsv',sep='\t')
groups={g:samples.loc[samples.condition.eq(g),'sample'].tolist() for g in ['CTRL','LOY']}
for kind,n in [('motifs',879),('archetypes',422)]:
    z=pd.read_csv(OUT/f'chromVAR/{kind}_zscores.csv').set_index('id')
    assert z.columns.tolist()==samples['sample'].tolist()
    assert len(z)==n and z.index.is_unique and np.isfinite(z.to_numpy()).all()
    tests=pd.read_csv(OUT/f'chromVAR/{kind}_independent_limma.csv').set_index('id').loc[z.index]
    effect=z[groups['LOY']].mean(axis=1)-z[groups['CTRL']].mean(axis=1)
    assert np.allclose(effect,tests.LOY_minus_CTRL_zscore,rtol=1e-9,atol=1e-9)
    assert np.allclose(false_discovery_control(tests['P.Value']),tests['adj.P.Val'],rtol=1e-10,atol=0)
    checks.append(f'PASS: {kind}: six-sample order, finite scores, effect direction, and independent-model BH correction ({n} tests)')

h=pd.read_csv(OUT/'HOMER/known_motif_all_results.csv')
cv=pd.read_csv(OUT/'chromVAR/motifs_independent_limma.csv')
for arm,tests in h.groupby('arm'):
    assert len(tests)==1758
    assert np.allclose(false_discovery_control(tests.PValue),tests.BH_FDR_both_directions,rtol=1e-10,atol=0)
    for direction in ['gained','lost']:
        folder=OUT/f'HOMER/{arm}'
        assert (folder/direction/'COMPLETE.json').exists()
        t=pd.read_csv(folder/f'{direction}_target.csv');b=pd.read_csv(folder/f'{direction}_background.csv')
        assert t.peak_id.is_unique and b.peak_id.is_unique and set(t.peak_id).isdisjoint(b.peak_id)
        assert t.FDR.lt(.05).all() and t.logFC.abs().gt(1).all() and b.FDR.ge(.2).all()
        assert (t.logFC.gt(0) if direction=='gained' else t.logFC.lt(0)).all()
        assert t.N_fraction.eq(0).all() and b.N_fraction.eq(0).all()
        strata=['gc_bin','accessibility_bin']
        tn=t.groupby(strata).size();bn=b.groupby(strata).size().reindex(tn.index)
        ratio=bn/tn
        assert ratio.notna().all() and ratio.nunique()==1
        assert bn.sum()==len(b)
        m=tests[tests.direction.eq(direction)]
        assert set(m.motif_id)==set(cv.id) and m.motif_id.is_unique
        checks.append(f'PASS: {arm} {direction}: {len(t)} targets and {len(b)} disjoint matched controls; all 879 motif identities agree')
    checks.append(f'PASS: {arm}: HOMER BH across both directions (1,758 tests)')
ref=json.loads((ROOT/'config/TF_extension/reference_verified.json').read_text())
assert ref['verified_peak_sequences']==84510 and ref['BAM_chromosome_lengths_agree'] and ref['line_endings']=='LF'
checks.append('PASS: recorded independent reference audit covers every tested peak and BAM chromosome lengths')
joint=h[h.arm.eq('primary_TMM')].merge(cv[['id','LOY_minus_CTRL_zscore','adj.P.Val']],left_on='motif_id',right_on='id',validate='many_to_one')
both=joint[(joint.BH_FDR_both_directions<.05)&(joint['adj.P.Val']<.05)].copy()
both['direction_agrees']=np.where(both.direction.eq('gained'),both.LOY_minus_CTRL_zscore>0,both.LOY_minus_CTRL_zscore<0)
both.to_csv(OUT/'concordance/HOMER_chromVAR_both_FDR005.csv',index=False)
checks.append(f'RESULT: {len(both)} motif models pass both primary HOMER and chromVAR FDR 0.05; {int(both.direction_agrees.sum())} have agreeing directions')
tf_rna=pd.read_csv(OUT/'concordance/chromVAR_TF_RNA_all.csv',dtype={'rna_gene_id':str})
source_rna=pd.read_csv(ROOT/'config/RNAseq/ht1080_batches_LOY_v_CTRL_DESeq2_batch_treatment.csv',dtype={'gene_id':str}).set_index('gene_id')
assert source_rna.index.is_unique
mapped=tf_rna[tf_rna.rna_gene_id.notna()]
original=source_rna.loc[mapped.rna_gene_id]
assert np.allclose(original.logFC,mapped.rna_log2FC,rtol=1e-12,atol=1e-14,equal_nan=True)
assert np.allclose(original.FDR,mapped.rna_FDR,rtol=1e-12,atol=0,equal_nan=True)
assert original.chr.astype(str).tolist()==mapped.rna_chromosome.astype(str).tolist()
finite=np.isfinite(mapped.rna_log2FC)&np.isfinite(mapped.rna_FDR)
assert mapped.RNA_statistic_status.eq('finite').to_numpy().tolist()==finite.to_numpy().tolist()
checks.append('PASS: TF RNA values and chromosomes match the unchanged selected source; nonfinite statistics remain distinguishable from absent mappings')

if not preliminary:
    import pyBigWig, pysam
    fp=OUT/'TOBIAS'
    assert (fp/'ANALYSIS_COMPLETE.json').exists() and (fp/'29_PlotAggregate.complete.json').exists()
    pools=pd.read_csv(fp/'replicate_pool_accounting.csv')
    assert set(pools['sample'])==set(samples['sample']) and pools['sample'].is_unique
    assert np.allclose(pools.actual_reads,pools.actual_fragments*2)
    totals=[];correction_qc=[]
    for condition in ['CTRL','LOY']:
        pool=json.loads((fp/f'{condition}_matched_pool.complete.json').read_text())
        assert pool['reads']==pools.loc[pools.condition.eq(condition),'actual_reads'].sum()
        totals.append(pool['reads'])
        for stage in ['ATACorrect','ScoreBigwig']:assert (fp/f'23_{stage}_{condition}.complete.json').exists()
        log=(ROOT/f'logs/23_ATACorrect_{condition}.log').read_text()
        metrics={key:float(re.search(r'\t'+key+r'\t([0-9.eE+-]+)',log).group(1))
                 for key in ['TOTAL_READS','PEAK_READS','NONPEAK_READS','FRiP']}
        assert metrics['TOTAL_READS']==pool['reads']
        assert metrics['PEAK_READS']+metrics['NONPEAK_READS']==metrics['TOTAL_READS']
        variances={f'{when}_{strand}':float(value) for when,strand,value in
                   re.findall(r'BIAS\t(pre|post)-bias variance (forward|reverse):\t([0-9.eE+-]+)',log)}
        assert len(variances)==4
        correction_qc.append(dict(condition=condition,input_reads=int(metrics['TOTAL_READS']),
            cut_site_FRiP=metrics['FRiP'],**variances,
            post_to_pre_bias_variance_forward=variances['post_forward']/variances['pre_forward'],
            post_to_pre_bias_variance_reverse=variances['post_reverse']/variances['pre_reverse']))
    pd.DataFrame(correction_qc).to_csv(fp/'ATACorrect_summary.csv',index=False)
    assert abs(totals[0]-totals[1])/min(totals)<.002
    track_qc=[]
    with pysam.FastaFile(ref['linux_cache']) as genome:
        for condition in ['CTRL','LOY']:
            for kind in ['uncorrected','corrected','footprints']:
                path=fp/f'{condition}/{condition}_{kind}.bw'
                with pyBigWig.open(str(path)) as bw:
                    header=bw.header();chroms=bw.chroms()
                    assert header['nBasesCovered']>0
                    assert all(np.isfinite(header[k]) for k in ['minVal','maxVal','sumData','sumSquared'])
                    for chrom,length in chroms.items():assert genome.get_reference_length(chrom)==length
                    track_qc.append(dict(condition=condition,track=kind,**header))
    pd.DataFrame(track_qc).to_csv(fp/'signal_track_QC.csv',index=False)
    tests=pd.read_csv(fp/'TF_scores_with_chromVAR.csv')
    assert tests.motif_id.is_unique and set(tests.motif_id)==set(cv.id)
    valid=np.isfinite(tests.LOY_CTRL_pvalue)&np.isfinite(tests.LOY_CTRL_change)
    assert np.allclose(false_discovery_control(tests.loc[valid,'LOY_CTRL_pvalue']),tests.loc[valid,'pooled_site_BH_FDR'],rtol=1e-10,atol=0)
    for name in ['TOBIAS_differential_volcano','CTCF_aggregate']:
        for ext in ['png','pdf']:assert (fp/f'{name}.{ext}').stat().st_size>1000
    assert (OUT/'concordance/TF_ATAC_footprint_RNA.csv').exists()
    checks.append(f'PASS: condition pools conserve reads ({totals[0]:,} CTRL; {totals[1]:,} LOY), with <0.2% depth difference')
    checks.append('PASS: all footprint stages and QC complete; motif identities and pooled-site BH correction agree')
    checks.append('PASS: all six uncorrected/corrected/footprint tracks have covered bases, finite summaries, and matching reference lengths')
else:
    checks.append('PENDING: footprint stages and aggregate QC; this is a preliminary audit')
checks.append('LIMITATION: technical-batch assignments remain unconfirmed; passage sensitivity is not labeled batch adjustment')
name='verification_preliminary.txt' if preliminary else 'verification.txt'
(OUT/name).write_text('\n'.join(checks)+'\n')
print('\n'.join(checks))
