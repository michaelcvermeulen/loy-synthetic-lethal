from pathlib import Path
import json
import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
PILOT=ROOT.parent/'WGS/reanalysis/Y_PAR_masked_20261027'
FIG=ROOT/'figures/chapter3/drafts/a549_dense_bins_20261027'
rows=[];allbins=[]
for size,minimum in [(100,50),(50,25)]:
    snap=FIG/f'{size}kb/source_data/regional_y_atlases_normalised/a549'
    z=pd.read_csv(snap/'paired_Y_bins_and_exclusions.csv')
    coverage=pd.read_csv(snap/'binned_coverage.csv')
    q=coverage.query("track=='Y' and sample=='A549-LOY'")[['end','read_count','aligned_bases']]
    z=z.merge(q,on='end',validate='one_to_one')
    z['size_kb']=size;allbins.append(z)
    for region in ['MSY Yp','euchromatic Yq']:
        v=z[z.included & z.region.eq(region)]
        rows.append(dict(size_kb=size,region=region,n_bins=len(v),genomic_span_Mb=(v.end-v.start+1).sum()/1e6,
            control_reads_min=v.control_primary_reads.min(),control_reads_median=v.control_primary_reads.median(),
            LOY_reads_median=v.read_count.median(),LOY_zero_read_bins=sum(v.read_count.eq(0)),
            LOY_zero_read_percent=100*v.read_count.eq(0).mean(),
            corrected_near_zero_bins=sum(v.copynumber_LOY.le(.001)),
            median_corrected_percent=v.percent_of_sgCTRL.median(),
            ratio_Q25=v.percent_of_sgCTRL.quantile(.25),ratio_Q75=v.percent_of_sgCTRL.quantile(.75),
            ratio_max=v.percent_of_sgCTRL.max()))
    if size==50:
        for region in ['MSY Yp','euchromatic Yq']:
            v=z[z.included & z.region.eq(region) & z.control_primary_reads.ge(50)]
            rows.append(dict(size_kb='50 (>=50 control reads)',region=region,n_bins=len(v),
                genomic_span_Mb=(v.end-v.start+1).sum()/1e6,
                control_reads_min=v.control_primary_reads.min(),control_reads_median=v.control_primary_reads.median(),
                LOY_reads_median=v.read_count.median(),LOY_zero_read_bins=sum(v.read_count.eq(0)),
                LOY_zero_read_percent=100*v.read_count.eq(0).mean(),
                corrected_near_zero_bins=sum(v.copynumber_LOY.le(.001)),
                median_corrected_percent=v.percent_of_sgCTRL.median(),
                ratio_Q25=v.percent_of_sgCTRL.quantile(.25),ratio_Q75=v.percent_of_sgCTRL.quantile(.75),
                ratio_max=v.percent_of_sgCTRL.max()))
pd.DataFrame(rows).to_csv(HERE/'resolution_comparison.csv',index=False)
pd.concat(allbins).to_csv(HERE/'all_Y_bins_with_resolution.csv',index=False)

# Exact count conservation check: sums of 50-kb read-start counts reproduce
# the saved 100-kb counts in every chromosome/bin/sample.
c100=pd.read_csv(ROOT/'research/model_development/a549_ypar_realign_20261027/masked_raw_counts.csv')
c50=pd.read_csv(HERE/'raw_counts_50kb.csv')
coordinates100=pd.read_csv(ROOT/'research/model_development/a549_ypar_realign_20261027/masked_bins.csv').query("sample=='A549-CTRL2'")[['bin_id','chromosome','start','end']]
coordinates50=pd.read_csv(HERE/'bins_50kb.csv').query("sample=='A549-CTRL2'")[['bin_id','chromosome','start','end']]
samples=['A549-CTRL2','A549-LOY']
a=coordinates100.merge(c100,on='bin_id',validate='one_to_one').set_index(['chromosome','start'])[samples].sort_index()
b=coordinates50.merge(c50,on='bin_id',validate='one_to_one')
b['start']=(b.start-1)//100000*100000+1
b=b.groupby(['chromosome','start'])[samples].sum().sort_index()
assert a.index.equals(b.index)
assert (a.to_numpy()==b.to_numpy()).all()
v100=pd.read_csv(PILOT/'a549_panel_tracks/binned_coverage.csv').query("track=='Y'")
v50=pd.read_csv(PILOT/'a549_50kb_comparison/binned_coverage.csv').query("track=='Y'")
v50['start0']=v50.start0//100000*100000
x=v100.set_index(['sample','start0']).aligned_bases.sort_index()
y=v50.groupby(['sample','start0']).aligned_bases.sum().sort_index()
assert x.index.equals(y.index) and (x.to_numpy()==y.to_numpy()).all()
(HERE/'resolution_validation.json').write_text(json.dumps({'all_100kb_raw_counts_reproduced_exactly_by_50kb_sums':True,
    'genomic_bins_checked':len(a),'samples':samples,'Y_primary_aligned_bases_conserved':True,
    'annotation_validation':json.loads((PILOT/'a549_50kb_comparison/annotation_validation.json').read_text())},indent=2))
print(pd.DataFrame(rows).to_string(index=False))
print('Exact raw-count conservation verified for',len(a),'100-kb bins and both samples.')
