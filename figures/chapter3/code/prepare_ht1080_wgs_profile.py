"""Freeze matched HT1080 WGS gene/segment copy ratios and ploidy for a figure.

Copy estimates follow the saved Y-dosage analysis: relative ratio * ploidy.
The original relative ratios remain in the output. No integer rounding is used.
"""
from pathlib import Path
import argparse
import hashlib
import json
import shutil
import numpy as np
import pandas as pd

MODEL = 'ACH-000054'
THESIS = Path(__file__).resolve().parents[3]


def digest(path):
    sha = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            sha.update(block)
    return sha.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--data-dir', type=Path, required=True)
    parser.add_argument('--baseline-dir', type=Path, default=THESIS/'figures/chapter3/source_data/cell_line_profiles')
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    raw = args.data_dir/'raw/depmap_public_26q1'
    inputs = {
        'gene_cn':raw/'copy_number/OmicsCNGeneWGS.csv',
        'segments':raw/'copy_number/OmicsCNSegmentsWGS.csv',
        'signatures':raw/'signatures/OmicsGlobalSignatures.csv',
        'profiles':raw/'metadata/OmicsProfiles.csv',
        'readme':raw/'documentation/README.txt',
        'saved_wgs_reference':args.data_dir/'derived/rna_y_vulnerabilities_v1/final/wgs_reference_profiles.csv',
    }
    manifest = pd.read_csv(args.data_dir/'docs/data_file_manifest.csv')
    records = {}
    for name, path in inputs.items():
        stat = path.stat()
        records[name] = {'path':str(path),'bytes':stat.st_size,'mtime_ns':stat.st_mtime_ns,'sha256':digest(path)}
        rel = path.relative_to(args.data_dir).as_posix()
        match = manifest.loc[manifest.new_relative_path.eq(rel)]
        if len(match):
            assert len(match)==1 and match.iloc[0]['bytes']==stat.st_size
            assert match.iloc[0].sha256==records[name]['sha256'], name
            records[name]['matches_hpc_manifest'] = True
    coords = pd.read_csv(args.baseline_dir/'reference_material/y_gene_coordinates_grch38.csv').sort_values('start')
    genes = coords.gene.tolist()
    metadata_cols = ['SequencingID','ModelConditionID','ModelID','IsDefaultEntryForMC','IsDefaultEntryForModel']
    cn = pd.read_csv(inputs['gene_cn'],usecols=lambda col: col in metadata_cols or col.split(' (')[0] in genes)
    cn = cn.loc[cn.ModelID.eq(MODEL) & cn.IsDefaultEntryForModel.eq('Yes')]
    assert len(cn)==1 and cn.IsDefaultEntryForMC.eq('Yes').all()
    row = cn.iloc[0]
    sequencing, condition = row.SequencingID, row.ModelConditionID
    profiles = pd.read_csv(inputs['profiles'])
    profile = profiles.loc[profiles.SequencingID.eq(sequencing)]
    assert len(profile)==1 and profile.DataType.eq('wgs').all()
    assert profile.ModelID.eq(MODEL).all() and profile.ModelConditionID.eq(condition).all()
    expr = pd.read_csv(args.baseline_dir/'ht1080_y_expression.csv')
    rna = profiles.loc[profiles.SequencingID.isin(expr.SequencingID.unique())]
    assert len(rna)==1 and rna.DataType.eq('rna').all() and rna.ModelConditionID.eq(condition).all()
    signatures = pd.read_csv(inputs['signatures'])
    signature = signatures.loc[signatures.SequencingID.eq(sequencing)]
    assert len(signature)==1 and signature.ModelID.eq(MODEL).all()
    ploidy = float(signature.iloc[0].Ploidy)
    assert np.isfinite(ploidy) and ploidy>0
    saved = pd.read_csv(inputs['saved_wgs_reference'])
    saved = saved.loc[saved.ModelID.eq(MODEL)&saved.SequencingID.eq(sequencing)&saved.chromosome.eq('Y')]
    assert len(saved)==1 and saved.iloc[0].ploidy_basis=='same_sequencing'
    assert np.isclose(saved.iloc[0].ploidy,ploidy)
    assert np.isclose(saved.iloc[0].copy_estimate,saved.iloc[0].ratio*ploidy)
    gene_rows = []
    for gene in genes:
        cols = [c for c in cn.columns if c.split(' (')[0]==gene]
        assert len(cols)==1
        value = float(row[cols[0]])
        assert np.isfinite(value) and value>=0
        gene_rows.append({'gene':gene,'source_column':cols[0],'ModelID':MODEL,
                          'SequencingID':sequencing,'ModelConditionID':condition,
                          'relative_copy_ratio':value,'matched_ploidy':ploidy,
                          'estimated_copy_number':value*ploidy})
    gene_table = pd.DataFrame(gene_rows)
    segments = pd.read_csv(inputs['segments'],low_memory=False)
    segments = segments.loc[segments.SequencingID.eq(sequencing)].copy()
    assert segments.ModelID.eq(MODEL).all() and segments.ModelConditionID.eq(condition).all()
    chroms = ['chr'+str(i) for i in range(1,23)]+['chrX','chrY']
    segments = segments.loc[segments.Chrom.isin(chroms)].drop(columns=['Unnamed: 0'],errors='ignore')
    assert set(segments.Chrom)==set(chroms)
    assert (segments.End>segments.Start).all() and np.isfinite(segments.RelCopyRatio).all()
    assert segments.RelCopyRatio.ge(0).all()
    segments['matched_ploidy'] = ploidy
    segments['estimated_copy_number'] = segments.RelCopyRatio*ploidy
    bands = pd.read_csv(args.baseline_dir/'reference_material/cytoBand_hg38.txt.gz',sep='\t',header=None,
                        names=['chromosome','start','end','band','stain'])
    lengths = bands.groupby('chromosome')['end'].max().reindex(chroms).rename('length_bp').reset_index()
    assert lengths.length_bp.notna().all()
    for chrom in chroms:
        s = segments[segments.Chrom.eq(chrom)].sort_values('Start')
        assert (s.Start.to_numpy()[1:]>=s.End.to_numpy()[:-1]).all(), chrom
        assert s.Start.min()>=0 and s.End.max()<=lengths.set_index('chromosome').loc[chrom,'length_bp']
    for name,path in inputs.items():
        assert (path.stat().st_size,path.stat().st_mtime_ns)==(records[name]['bytes'],records[name]['mtime_ns']), name
    out = args.output_dir
    out.mkdir(parents=True,exist_ok=False)
    for name in ['ht1080_background.json','ht1080_y_expression.csv','ht1080_dependency_context.csv','figure_config.json']:
        shutil.copy2(args.baseline_dir/name,out/name)
    shutil.copytree(args.baseline_dir/'reference_material',out/'reference_material')
    gene_table.to_csv(out/'ht1080_y_gene_wgs_cn.csv',index=False)
    segments.to_csv(out/'ht1080_genome_wgs_segments.csv',index=False)
    lengths.to_csv(out/'chromosome_lengths_grch38.csv',index=False)
    pd.concat([profile,rna],ignore_index=True).to_csv(out/'matched_omics_profiles.csv',index=False)
    signature.to_csv(out/'matched_wgs_signatures.csv',index=False)
    saved.to_csv(out/'saved_y_reference_check.csv',index=False)
    meta = {'model':MODEL,'WGS_SequencingID':sequencing,'RNA_SequencingID':rna.iloc[0].SequencingID,
            'ModelConditionID':condition,'same_culture_condition':True,'matched_ploidy':ploidy,
            'copy_number_scale':'Estimated copies = released linear relative copy ratio multiplied by ploidy from the same WGS sequencing',
            'interpretation':'Continuous estimates, not integer or allele-specific calls; RNA and WGS were collected at different times',
            'gene_order':genes,'n_genes':len(gene_table),'n_segments':len(segments),'n_chromosomes':24,
            'genome_reference':'GRCh38','inputs':records}
    (out/'wgs_cn_provenance.json').write_text(json.dumps(meta,indent=2)+'\n',encoding='utf-8')
    print(gene_table[['gene','relative_copy_ratio','estimated_copy_number']].to_string(index=False))
    print(segments.groupby('Chrom').estimated_copy_number.agg(['count','min','max']).to_string())
    print('Ploidy:',ploidy,'; segments:',len(segments))


if __name__=='__main__':
    main()
