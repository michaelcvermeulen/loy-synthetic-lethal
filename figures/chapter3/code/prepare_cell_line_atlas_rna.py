"""Rebuild descriptive RNA tables from recorded raw-count paths and gene mappings.
Run after initial source registration; preserves metadata/sample selection.
"""
from pathlib import Path
import pandas as pd
import numpy as np
ROOT=Path(__file__).resolve().parents[3]
for cell in ['a549','du145']:
 folder=ROOT/f'figures/chapter3/source_data/{cell}_y_expression_confirmation'
 meta=pd.read_csv(folder/'sample_metadata.csv')
 mapping=pd.read_csv(folder/'gene_mapping_plotted.csv')
 result=[]
 for row in meta.itertuples():
  counts=pd.read_csv(row.source_file,sep='\t',comment='#')
  assert counts.Geneid.is_unique
  value=counts.columns[-1]
  total=counts[value].sum()
  assert total==row.library_size_all_genes and (counts[value]>=0).all()
  selected=counts.set_index('Geneid').loc[mapping.gene_id].reset_index()
  part=mapping.copy()
  part['count_coordinates_chr']=selected.Chr.to_numpy()
  assert all(set(ch.split(';'))=={expected} for ch,expected in zip(part.count_coordinates_chr,part.count_chromosome))
  part['raw_count']=selected[value].to_numpy()
  part['cpm']=part.raw_count/total*1e6
  part['log2_cpm_plus_1']=np.log2(part.cpm+1)
  for key in meta.columns:part[key]=getattr(row,key)
  result.append(part)
 combined=pd.concat(result,ignore_index=True)
 previous=pd.read_csv(folder/'Y_PAR_per_sample.csv')
 keys=['sample','gene_name']
 old=previous.set_index(keys).sort_index();new=combined.set_index(keys).sort_index()
 assert old.index.equals(new.index)
 assert np.allclose(old[['raw_count','cpm','log2_cpm_plus_1']],new[['raw_count','cpm','log2_cpm_plus_1']])
 combined.to_csv(folder/'Y_PAR_per_sample.csv',index=False)
 print(cell,'raw-count reconstruction verified',len(combined),'rows')
