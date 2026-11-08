"""Freeze public 26Q1 data for the compact A549/DU145 baseline profiles."""
from pathlib import Path
import json
import shutil
import numpy as np
import pandas as pd
from prepare_cell_line_profiles import digest, symbol

THESIS = Path(__file__).resolve().parents[3]
ROOT = Path('H:/PRISM/data/CCLE/26Q1')
RAW = ROOT/'raw/depmap_public_26q1'
BASE = THESIS/'figures/chapter3/source_data/cell_line_profiles'
OUT = THESIS/'figures/chapter3/source_data/a549_du145_profiles'
MODELS = {'A549':'ACH-000681', 'DU145':'ACH-000979'}

def main():
    OUT.mkdir(exist_ok=False)
    inputs = {'gene_cn':RAW/'copy_number/OmicsCNGeneWGS.csv',
              'segments':RAW/'copy_number/OmicsCNSegmentsWGS.csv',
              'profiles':RAW/'metadata/OmicsProfiles.csv',
              'signatures':RAW/'signatures/OmicsGlobalSignatures.csv',
              'effects':RAW/'crispr/CRISPRGeneEffect.csv',
              'expression':BASE/'selected_models_y_expression.csv',
              'metadata':BASE/'selected_model_metadata.csv'}
    manifest = pd.read_csv(ROOT/'docs/data_file_manifest.csv')
    provenance = {'release':'DepMap Public 26Q1','inputs':{},'models':{}}
    for key,path in inputs.items():
        rec = {'path':str(path),'bytes':path.stat().st_size,'sha256':digest(path)}
        if path.is_relative_to(ROOT):
            expected = manifest.loc[manifest.new_relative_path.eq(path.relative_to(ROOT).as_posix())]
            assert len(expected)==1 and expected.iloc[0].sha256==rec['sha256'], key
            rec['matches_hpc_manifest']=True
        provenance['inputs'][key]=rec
    shutil.copytree(BASE/'reference_material',OUT/'reference_material',ignore=shutil.ignore_patterns('*.jpg'))
    shutil.copy2(THESIS/'figures/chapter3/source_data/ht1080_profile_wgs_cn/chromosome_lengths_grch38.csv',OUT)
    coords = pd.read_csv(OUT/'reference_material/y_gene_coordinates_grch38.csv').sort_values('start')
    genes = coords.gene.tolist()
    ids = ['SequencingID','ModelConditionID','ModelID','IsDefaultEntryForModel','IsDefaultEntryForMC']
    cn = pd.read_csv(inputs['gene_cn'],usecols=lambda c:c in ids or symbol(c) in genes)
    cn = cn.loc[cn.ModelID.isin(MODELS.values()) & cn.IsDefaultEntryForModel.eq('Yes')]
    profiles = pd.read_csv(inputs['profiles'])
    sig = pd.read_csv(inputs['signatures'])
    segments = pd.read_csv(inputs['segments'],low_memory=False)
    expr = pd.read_csv(inputs['expression'])
    expr = expr.loc[expr.ModelID.isin(MODELS.values())]
    expr.to_csv(OUT/'y_expression.csv',index=False)
    pd.read_csv(inputs['metadata']).query('ModelID in @MODELS.values()').to_csv(OUT/'model_metadata.csv',index=False)
    effects = pd.read_csv(inputs['effects'],index_col=0)
    effects.columns = [symbol(c) for c in effects.columns]
    assert effects.index.is_unique and effects.columns.is_unique
    gene_tables, segment_tables, dep_tables, profile_tables = [], [], [], []
    chroms = ['chr'+str(i) for i in range(1,23)]+['chrX','chrY']
    for line,model in MODELS.items():
        rows = cn.loc[cn.ModelID.eq(model)]
        assert len(rows)==1
        row = rows.iloc[0]
        signature = sig.loc[sig.SequencingID.eq(row.SequencingID)]
        assert len(signature)==1 and signature.ModelID.eq(model).all()
        ploidy = float(signature.iloc[0].Ploidy)
        assert np.isfinite(ploidy) and ploidy>0
        rna_id = expr.loc[expr.ModelID.eq(model),'SequencingID'].unique()
        assert len(rna_id)==1
        matched = profiles.loc[profiles.SequencingID.isin([row.SequencingID,rna_id[0]])]
        assert len(matched)==2 and matched.ModelID.eq(model).all()
        assert matched.ModelConditionID.eq(row.ModelConditionID).all()
        assert set(matched.DataType)=={'rna','wgs'}
        profile_tables.append(matched)
        for gene in genes:
            col = [c for c in rows.columns if symbol(c)==gene]
            assert len(col)==1
            value=float(row[col[0]])
            assert np.isfinite(value) and value>=0
            gene_tables.append({'cell_line':line,'ModelID':model,'SequencingID':row.SequencingID,
                                'gene':gene,'relative_copy_ratio':value,'matched_ploidy':ploidy,
                                'estimated_copy_number':value*ploidy})
        seg = segments.loc[segments.SequencingID.eq(row.SequencingID) & segments.Chrom.isin(chroms)].copy()
        assert set(seg.Chrom)==set(chroms) and seg.ModelID.eq(model).all()
        assert seg.ModelConditionID.eq(row.ModelConditionID).all()
        assert seg.RelCopyRatio.ge(0).all() and (seg.End>seg.Start).all()
        seg['cell_line']=line
        seg['matched_ploidy']=ploidy
        seg['estimated_copy_number']=seg.RelCopyRatio*ploidy
        segment_tables.append(seg)
        ref = effects.drop(index=model)
        deps = pd.DataFrame({'gene_effect':effects.loc[model], 'reference_median':ref.median(),
             'reference_q25':ref.quantile(.25),'reference_q75':ref.quantile(.75),
             'n_reference':ref.count(),'reference_dependency_fraction':ref.lt(-.5).sum()/ref.count()})
        deps.index.name='gene'
        deps['difference_from_reference_median']=deps.gene_effect-deps.reference_median
        deps=deps.sort_values('difference_from_reference_median')
        candidates=deps.loc[deps.gene_effect.lt(-.5)&deps.reference_dependency_fraction.le(.20)]
        deps.to_csv(OUT/f'{line.lower()}_dependency_ranking.csv')
        top=candidates.head(6).reset_index()
        assert len(top)==6
        top['cell_line']=line
        dep_tables.append(top)
        provenance['models'][line]={'ModelID':model,'WGS_SequencingID':row.SequencingID,
          'RNA_SequencingID':rna_id[0],'ModelConditionID':row.ModelConditionID,
          'matched_ploidy':ploidy,'same_model_condition':True,'n_segments':len(seg),
          'n_dependency_candidates':len(candidates),'displayed_dependencies':top.gene.tolist()}
    pd.DataFrame(gene_tables).to_csv(OUT/'y_gene_wgs_cn.csv',index=False)
    pd.concat(segment_tables).to_csv(OUT/'genome_wgs_segments.csv',index=False)
    pd.concat(dep_tables).to_csv(OUT/'selected_dependencies.csv',index=False)
    pd.concat(profile_tables).to_csv(OUT/'matched_omics_profiles.csv',index=False)
    provenance['dependency_selection']='Per model: gene effect < -0.5; <=20% of other screened models with effect < -0.5; six largest negative differences from reference median. Both sexes/all lineages; no common-essential exclusion; descriptive only.'
    provenance['copy_number']='Released linear relative copy ratio multiplied by ploidy from the same WGS sequencing profile. Continuous estimates, not integer calls. RNA and WGS share model condition but not necessarily collection date.'
    (OUT/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(provenance['models'],indent=2))

if __name__=='__main__':
    main()
