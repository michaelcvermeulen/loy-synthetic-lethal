"""Compare reference bin filters without selecting on the LOY outcome."""
from pathlib import Path
import pandas as pd

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
a=pd.read_csv(ROOT/'research/model_development/wgs_filter_rerun_20261027/bin_annotation.csv')
c=pd.read_csv(ROOT/'research/model_development/a549_ypar_realign_20261027/masked_raw_counts.csv')
d=a.merge(c,on='bin_id',validate='one_to_one')
y=d.loc[d.chromosome.eq('Y')].copy()
y['region']='excluded'
y.loc[(y.start>2458320)&(y.end<=10565750),'region']='MSY Yp'
y.loc[(y.start>10883085)&(y.end<=27449937),'region']='euchromatic Yq'
y=y[y.region.ne('excluded')]
rows=[]
for m in [55,40,30,25,10,0]:
    for bl in [20,50,80,100.01]:
        keep=(y.bases>=60)&(y.mappability>=m)&(y.blacklist<bl)&(y['A549-CTRL2']>=20)
        for region in ['MSY Yp','euchromatic Yq']:
            z=y[keep & y.region.eq(region)]
            rows.append(dict(mappability_min=m,difficulty_max_exclusive=bl,control_reads_min=20,
                region=region,n_bins=len(z),control_min=z['A549-CTRL2'].min(),
                control_median=z['A549-CTRL2'].median()))
pd.DataFrame(rows).to_csv(HERE/'filter_eligibility_comparison.csv',index=False)
y.to_csv(HERE/'all_candidate_Y_bins.csv',index=False)
print(pd.DataFrame(rows).query('mappability_min in [55,30]').to_string(index=False))
