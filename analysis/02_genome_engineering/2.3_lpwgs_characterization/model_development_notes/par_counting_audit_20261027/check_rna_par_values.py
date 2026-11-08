"""Verify atlas PAR values against recorded featureCounts files; no input edits."""
from pathlib import Path
import json
import numpy as np
import pandas as pd

OUT=Path(__file__).resolve().parent
ROOT=OUT.parents[2]
sources=ROOT/'figures/chapter3/source_data'
checks=[]
for cell in ['ht1080','a549','du145']:
    data=pd.read_csv(sources/f'{cell}_y_expression_confirmation/relative_expression_per_sample.csv')
    par=data[data.gene_name.isin(['SLC25A6','CD99','VAMP7'])]
    for sample,rows in par.groupby('sample'):
        source=Path(rows.source_file.iloc[0])
        if not source.is_absolute():
            source=ROOT.parent/'thesis/Ben_Thesis'/source
        counts=pd.read_csv(source,sep='\t',comment='#')
        assert counts.Geneid.is_unique
        values=counts.set_index('Geneid').iloc[:,-1]
        for row in rows.itertuples():
            assert row.raw_count==values.loc[row.gene_id]
            assert np.isclose(row.cpm,1e6*row.raw_count/values.sum())
            assert np.isclose(row.relative_expression_percent,100*row.cpm/row.control_mean_cpm)
            checks.append({'cell':cell,'sample':sample,'gene':row.gene_name,
                           'raw_count':row.raw_count,'CPM':row.cpm,
                           'percent_of_sgCTRL':row.relative_expression_percent,
                           'count_chromosome':row.count_chromosome,'source_file':str(source)})
pd.DataFrame(checks).to_csv(OUT/'rna_PAR_verified_observations.csv',index=False)
result={'n_observations':len(checks),'n_samples':26,'genes':['SLC25A6','CD99','VAMP7'],
        'raw_counts_match':True,'CPM_and_ratios_match':True,
        'all_count_annotations_on_chrX':all(x['count_chromosome']=='chrX' for x in checks),
        'scope':'Checks table extraction and arithmetic; does not establish correctness of original alignment or PAR re-quantification.'}
(OUT/'rna_PAR_checks.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
print(json.dumps(result,indent=2))
