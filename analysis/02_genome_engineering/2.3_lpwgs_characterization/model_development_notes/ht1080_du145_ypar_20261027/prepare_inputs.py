from pathlib import Path
import argparse,json,shutil
import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
PILOT=ROOT.parent/'WGS/reanalysis/Y_PAR_masked_20261027'
ap=argparse.ArgumentParser();ap.add_argument('cell',choices=['HT1080','DU145']);args=ap.parse_args()
cell=args.cell
samples={'HT1080':['HT1080-sgCTRL-D0','HT1080-LOY-D0'],'DU145':['LNCaP-sgCTRL','LNCaP-LOY-2X']}[cell]
for s in samples:assert (PILOT/'qc'/f'{s}.REALIGNMENT_COMPLETE.txt').exists()
out=HERE/cell/'plot_inputs';out.mkdir(parents=True,exist_ok=True)
d=pd.read_csv(HERE/cell/'masked/bins.csv')
assert set(d['sample'])==set(samples)
d['chromosome']='chr'+d.chromosome.astype(str)
d['label']=d['sample'].map(dict(zip(samples,['sgCTRL','LOY'])))
d['midpoint']=(d.start+d.end)/2e6
d['log2']=np.log2(d.copynumber.clip(lower=.001))
d['log2_seg']=np.log2(d.segmented.clip(lower=.001))
d.to_csv(out/'qdnaseq_bins.csv',index=False)
cov=PILOT/'atlas_pair_coverage'/cell/'masked'
for name in ['binned_coverage.csv','autosomal_normalization.csv','regional_PAR_coverage.csv']:
    shutil.copy2(cov/name,out/name)
ann=ROOT/'research/model_development/wgs_filter_rerun_20261027/bin_annotation.csv'
shutil.copy2(ann,out/'bin_annotation.csv')
(out/'sample_map.json').write_text(json.dumps({'cell':cell,'control':samples[0],'LOY':samples[1]},indent=2))
print('Plot inputs ready:',cell)
