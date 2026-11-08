"""Extract five default public WGS profiles without pooling assays or stocks."""
from pathlib import Path
import json
import shutil
import numpy as np
import pandas as pd
from prepare_cell_line_profiles import digest, symbol

THESIS=Path(__file__).resolve().parents[3]
ROOT=Path('H:/PRISM/data/CCLE/26Q1')
RAW=ROOT/'raw/depmap_public_26q1'
OUT=THESIS/'figures/chapter3/source_data/selected_models_y_cn'
MODELS={'HT1080':'ACH-000054','A549':'ACH-000681','DU145':'ACH-000979',
        'H2452':'ACH-000092','HepG2':'ACH-000739'}

def main():
    OUT.mkdir(exist_ok=False)
    paths={'profiles':RAW/'metadata/OmicsProfiles.csv','signatures':RAW/'signatures/OmicsGlobalSignatures.csv',
      'gene_cn':RAW/'copy_number/OmicsCNGeneWGS.csv','segments':RAW/'copy_number/OmicsCNSegmentsWGS.csv',
      'metadata':RAW/'metadata/Model.csv',
      'reviewed_calls':ROOT/'derived/sex_chromosome_review_v1/sex_chromosome_calls_reviewed.csv'}
    manifest=pd.read_csv(ROOT/'docs/data_file_manifest.csv')
    provenance={'release':'DepMap Public 26Q1','inputs':{},'models':{}}
    for key,path in paths.items():
        rec={'path':str(path),'bytes':path.stat().st_size,'sha256':digest(path)}
        if path.is_relative_to(RAW):
            expected=manifest.loc[manifest.new_relative_path.eq(path.relative_to(ROOT).as_posix())]
            assert len(expected)==1 and expected.iloc[0].sha256==rec['sha256']
            rec['matches_hpc_manifest']=True
        provenance['inputs'][key]=rec
    reference=THESIS/'figures/chapter3/source_data/cell_line_profiles/reference_material'
    for file in ['y_gene_coordinates_grch38.csv','chrY_cytobands_grch38.csv']:
        shutil.copy2(reference/file,OUT/file)
    coords=pd.read_csv(OUT/'y_gene_coordinates_grch38.csv').sort_values('start')
    genes=coords.gene.tolist()
    profiles=pd.read_csv(paths['profiles'])
    sig=pd.read_csv(paths['signatures'])
    segments=pd.read_csv(paths['segments'],low_memory=False)
    ids=['SequencingID','ModelID','ModelConditionID','IsDefaultEntryForModel']
    cn=pd.read_csv(paths['gene_cn'],usecols=lambda c:c in ids or symbol(c) in genes)
    gene_rows=[];segment_rows=[];profile_rows=[];signature_rows=[]
    for line,model in MODELS.items():
        prof=profiles.loc[profiles.ModelID.eq(model)&profiles.DataType.eq('wgs')&profiles.IsDefaultEntryForModel.eq('Yes')]
        assert len(prof)==1
        p=prof.iloc[0];s=sig.loc[sig.SequencingID.eq(p.SequencingID)]
        assert len(s)==1 and s.ModelID.eq(model).all()
        ploidy=float(s.iloc[0].Ploidy);assert np.isnan(ploidy) or ploidy>0
        c=cn.loc[cn.SequencingID.eq(p.SequencingID)]
        assert len(c)==1 and c.ModelID.eq(model).all() and c.ModelConditionID.eq(p.ModelConditionID).all()
        for gene in genes:
            col=[x for x in c if symbol(x)==gene];assert len(col)==1
            value=float(c.iloc[0][col[0]]);assert np.isfinite(value) and value>=0
            gene_rows.append({'cell_line':line,'ModelID':model,'SequencingID':p.SequencingID,
                              'gene':gene,'relative_copy_ratio':value,'matched_ploidy':ploidy,
                              'estimated_copy_number':value*ploidy})
        seg=segments.loc[segments.SequencingID.eq(p.SequencingID)&segments.Chrom.eq('chrY')].copy()
        assert len(seg)>0 and seg.ModelID.eq(model).all() and seg.RelCopyRatio.ge(0).all()
        seg=seg.sort_values('Start')
        assert (seg.Start.to_numpy()[1:]>=seg.End.to_numpy()[:-1]).all()
        seg['cell_line']=line;seg['matched_ploidy']=ploidy
        seg['estimated_copy_number']=seg.RelCopyRatio*ploidy
        segment_rows.append(seg);profile_rows.append(prof);signature_rows.append(s)
        provenance['models'][line]={'ModelID':model,'SequencingID':p.SequencingID,
           'ModelConditionID':p.ModelConditionID,'sequencing_date':p.SequencingDate,
           'ploidy':ploidy if np.isfinite(ploidy) else None,'segment_count':len(seg),'covered_start':int(seg.Start.min()),
           'covered_end':int(seg.End.max())}
    pd.DataFrame(gene_rows).to_csv(OUT/'y_gene_wgs_cn.csv',index=False)
    pd.concat(segment_rows).to_csv(OUT/'chrY_wgs_segments.csv',index=False)
    pd.concat(profile_rows).to_csv(OUT/'wgs_profiles.csv',index=False)
    pd.concat(signature_rows).to_csv(OUT/'wgs_signatures.csv',index=False)
    for key in ['metadata','reviewed_calls']:
        df=pd.read_csv(paths[key],low_memory=False)
        df.loc[df.ModelID.isin(MODELS.values())].to_csv(OUT/f'{key}.csv',index=False)
    provenance['units']='Plot uses released linear relative copy ratio for all five lines; matched HepG2 ploidy is unavailable. Estimated copies retained in source tables only where same-sequencing ploidy exists; no integer rounding or ploidy imputation.'
    provenance['scope']='Public baseline models, not the experimental stocks or edited populations. Gene-level and segment-level outputs shown separately.'
    (OUT/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(provenance['models'],indent=2))
    print(pd.DataFrame(gene_rows).pivot(index='cell_line',columns='gene',values='estimated_copy_number').round(3).to_string())

if __name__=='__main__':main()
