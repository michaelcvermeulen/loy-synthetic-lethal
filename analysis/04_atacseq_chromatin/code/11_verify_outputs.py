"""Check data conservation, coordinate invariants, and completed analysis outputs."""
import hashlib
import base64
import json
from pathlib import Path
import re
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


if __name__ == '__main__':
    samples = pd.read_csv(ROOT/'config/samples.tsv',sep='\t')['sample'].tolist()
    checks = []
    for sample in samples:
        qc = json.loads((ROOT/f'qc/{sample}_preparation.json').read_text())
        get = lambda key: qc.get(key,0)
        excluded_records = sum(get(k) for k in ['unmatched_primary_records','unmapped_records',
            'noncanonical_records','nonprimary_records','improper_records'])
        assert qc['source_records'] == 2*qc['complete_primary_pairs']+excluded_records, sample
        rejected_pairs = sum(get(k) for k in ['duplicate_pairs','qc_failed_pairs','low_mapq_pairs',
            'invalid_geometry_pairs','inconsistent_template_length_pairs','blacklist_pairs'])
        assert qc['complete_primary_pairs'] == qc['accepted_fragments']+rejected_pairs, sample
        assert sum(qc['chromosome_fragments'].values()) == qc['accepted_fragments'], sample
        log = (ROOT/f'logs/{sample}_macs3.log').read_text()
        parsed = re.search(r'total fragments in treatment:\s*(\d+)',log)
        assert parsed and int(parsed.group(1)) == qc['accepted_fragments'], sample
        checks.append(f'{sample}: record/pair conservation and MACS3 ingestion passed')
    peaks = pd.read_csv(ROOT/'peaks/consensus_membership.tsv',sep='\t')
    counts = pd.read_csv(ROOT/'counts/peak_counts.tsv',sep='\t')
    assert peaks.peak_id.is_unique and peaks.peak_id.equals(counts.peak_id)
    assert counts.columns.tolist() == ['peak_id']+samples
    values = counts[samples].to_numpy()
    assert np.all(values >= 0) and np.all(values == np.floor(values))
    assert (peaks.end-peaks.start).eq(501).all()
    assert ((peaks.support_CTRL >= 2) | (peaks.support_LOY >= 2)).all()
    for chrom, group in peaks.groupby('chrom'):
        ordered = group.sort_values('start')
        assert np.all(ordered.start.to_numpy()[1:] >= ordered.end.to_numpy()[:-1]), chrom
    qc = pd.read_csv(ROOT/'qc/peak_qc.tsv',sep='\t')
    assert qc.frip_own.between(0,1).all() and qc.frip_consensus.between(0,1).all()
    assert (qc.extra_consensus_assignments >= 0).all()
    checks.append('Consensus coordinates, replicate support, integer counts, and unique-fragment FRiP passed')
    primary_ids = None
    for arm in ['primary_TMM','background_TMM','passage_TMM']:
        result = pd.read_csv(ROOT/f'differential/{arm}_all.csv')
        assert result.peak_id.is_unique and result.FDR.between(0,1).all()
        assert result.PValue.between(0,1).all() and np.isfinite(result.logFC).all()
        assert result['TREAT_FDR_lfc0.5'].between(0,1).all()
        assert not result.chrom.eq('chrY').any() and not result.is_PAR.any()
        if primary_ids is None: primary_ids = result.peak_id.tolist()
        assert result.peak_id.tolist() == primary_ids
    checks.append('Primary/sensitivity models share tested regions and exclude chrY/PAR; statistics are finite')
    depth = pd.read_csv(ROOT/'qc/matched_depth_peak_counts.tsv',sep='\t')
    assert sorted(depth['sample']) == sorted(samples) and depth.matched_fragments.nunique() == 1
    target_depth = int(depth.matched_fragments.iloc[0])
    original_depths = {s: json.loads((ROOT/f'qc/{s}_preparation.json').read_text())[
        'accepted_fragments'] for s in samples}
    assert target_depth == min(original_depths.values())
    for sample in samples:
        log = (ROOT/f'logs/{sample}_matched_depth_macs3.log').read_text()
        parsed = re.search(r'total fragments in treatment:\s*(\d+)', log)
        assert parsed and int(parsed.group(1)) == target_depth, sample
        if original_depths[sample] == target_depth:
            full = ROOT/f'peaks/{sample}/{sample}_peaks.narrowPeak'
            matched = ROOT/f'peaks/matched_depth/{sample}/{sample}_peaks.narrowPeak'
            assert full.read_bytes() == matched.read_bytes(), sample
    checks.append('Matched-depth comparison includes all six samples at exactly equal fragment depth')
    checks.append('MACS3 matched-depth ingestion verified; unchanged-depth BAMPE/BEDPE peaks are byte-identical')
    required = ['annotation/annotated_differential_peaks.csv','annotation/promoter_peak_gene_links.csv',
      'annotation/promoter_gene_summary.csv',
      'differential/csaw_regions.csv','differential/primary_TMM_csaw_support.csv',
      'motifs/matching_accounting.csv','motifs/METHODS.txt','report/ATACseq_report.html']
    for name in required: assert (ROOT/name).is_file(), name
    genes = pd.read_csv(ROOT/'annotation/promoter_gene_summary.csv')
    assert genes.gene_id.is_unique
    motif_file = ROOT/'motifs/motif_enrichment.csv'
    if motif_file.exists():
        motifs = pd.read_csv(motif_file)
        assert motifs.PValue.between(0,1).all() and motifs.FDR.between(0,1).all()
        assert not motifs.duplicated(['motif_id','direction']).any()
    checks.append('Required downstream tables and standalone report exist')
    rna_dir = ROOT/'integration/RNAseq'
    if (rna_dir/'summary.json').exists():
        manifest = pd.read_csv(rna_dir/'source_manifest.csv')
        for row in manifest.itertuples():
            assert hashlib.sha256((ROOT/row.snapshot).read_bytes()).hexdigest() == row.sha256
        source_path = ROOT/manifest.loc[manifest.role == 'result','snapshot'].iloc[0]
        source = pd.read_csv(source_path, dtype={'entrezid':str})
        source = source[source.entrezid.notna() & ~source.entrezid.duplicated(keep=False)].set_index('entrezid')
        genes = pd.read_csv(rna_dir/'ATAC_RNA_gene_universe.csv', dtype={'entrez_id':str})
        assert genes.entrez_id.is_unique and genes.mapping_status.eq('eligible').all()
        assert genes.atac_chrom.eq(genes.rna_chrom).all() and not genes.rna_PAR.any()
        original = source.loc[genes.entrez_id]
        for old, new in [('logFC','rna_logFC'), ('FDR','rna_FDR'), ('PValue','rna_PValue'), ('logCPM','rna_baseMean')]:
            assert np.allclose(original[old].to_numpy(), genes[new].to_numpy(), rtol=1e-12, atol=1e-14)
        audit = pd.read_csv(rna_dir/'gene_mapping_audit.csv')
        assert len(audit) == len(pd.read_csv(ROOT/'annotation/promoter_gene_summary.csv'))
        joint = pd.read_csv(rna_dir/'ATAC_RNA_joint_significant_genes.csv')
        assert joint.rna_FDR.lt(.05).all() and joint.atac_any_significant_promoter_peak.all()
        assert len(pd.read_csv(rna_dir/'ATAC_RNA_promoter_peak_links.csv')) == len(pd.read_csv(ROOT/'annotation/promoter_peak_gene_links.csv'))
        checks.append('RNA source hashes, unchanged RNA statistics, unique common gene universe, and promoter-link accounting passed')
    tf_dir = ROOT/'motifs/TF_extension'
    if (tf_dir/'TOBIAS/summary.json').exists():
        assert (tf_dir/'verification.txt').exists()
        html=(ROOT/'report/ATACseq_report.html').read_text(encoding='utf-8')
        assert 'TOBIAS footprinting is pending' not in html
        assert 'Original motifmatchr/JASPAR2022 enrichment screen' in html
        tf_figures=['HOMER/primary_TMM/DA_volcano_MA.png','HOMER/known_motif_enrichment.png',
            'chromVAR/motifs_independent_heatmap.png','chromVAR/archetypes_independent_heatmap.png',
            'concordance/chromVAR_TF_RNA.png','concordance/HOMER_chromVAR.png',
            'TOBIAS/TOBIAS_differential_volcano.png','concordance/chromVAR_TOBIAS.png','TOBIAS/CTCF_aggregate.png']
        for name in tf_figures:
            image=tf_dir/name
            assert image.is_file() and base64.b64encode(image.read_bytes()).decode() in html,name
        checks.append('All nine current TF extension figures are embedded in the standalone report; footprint completion status agrees')
    records = []
    for directory in ['code','config']:
        for path in sorted((ROOT/directory).rglob('*')):
            if path.is_file() and '__pycache__' not in path.parts:
                digest=hashlib.sha256()
                with path.open('rb') as handle:
                    for block in iter(lambda:handle.read(1024*1024),b''):digest.update(block)
                records.append({'path':str(path.relative_to(ROOT)),
                    'sha256':digest.hexdigest()})
    pd.DataFrame(records).to_csv(ROOT/'logs/code_config_checksums.csv',index=False)
    (ROOT/'report/verification.txt').write_text('\n'.join(checks)+'\n')
    print('\n'.join(checks))
