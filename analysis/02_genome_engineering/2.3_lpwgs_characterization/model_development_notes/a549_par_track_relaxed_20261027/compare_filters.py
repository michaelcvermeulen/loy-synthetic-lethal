from pathlib import Path
import pandas as pd
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
PILOT=ROOT.parent/'WGS/reanalysis/Y_PAR_masked_20261027'
ann=pd.read_csv(ROOT/'research/model_development/wgs_filter_rerun_20261027/bin_annotation.csv')
coverage=pd.read_csv(PILOT/'a549_panel_tracks/binned_coverage.csv')
ctrl=coverage.query("sample=='A549-CTRL2' and track=='Y'")[['end','read_count','depth']]
specs=[('strict55_20',ROOT/'research/model_development/a549_ypar_realign_20261027/masked_bins.csv')]+[(name,HERE/name/'bins.csv') for name in ['map30_difficulty50','map30_difficulty80','map30_no_difficulty']]
rows=[];audit=[]
for name,path in specs:
    d=pd.read_csv(path).query("chromosome == 'Y'")
    z=d[d['sample'].eq('A549-CTRL2')].merge(d[d['sample'].eq('A549-LOY')],on=['bin_id','chromosome','start','end'],suffixes=('_control','_LOY'),validate='one_to_one')
    z=z.merge(ctrl,on='end',validate='one_to_one')
    z['region']='excluded'
    z.loc[(z.start>2458320)&(z.end<=10565750),'region']='MSY Yp'
    z.loc[(z.start>10883085)&(z.end<=27449937),'region']='euchromatic Yq'
    z['included']=z.use_control & z.use_LOY & z.region.ne('excluded') & z.read_count.ge(50) & z.copynumber_control.gt(.001)
    z['ratio']=100*z.copynumber_LOY.clip(lower=0)/z.copynumber_control
    z['filter']=name
    audit.append(z)
    for region in ['MSY Yp','euchromatic Yq']:
        zz=z[z.included & z.region.eq(region)]
        assert np.isfinite(zz.ratio).all()
        rows.append({'filter':name,'region':region,'n_bins':len(zz),'median_percent':zz.ratio.median(),
            'control_reads_min':zz.read_count.min(),'control_reads_median':zz.read_count.median(),
            'ratio_Q25':zz.ratio.quantile(.25),'ratio_Q75':zz.ratio.quantile(.75),
            'ratio_max':zz.ratio.max(),'n_above_100_percent':sum(zz.ratio>100)})
pd.DataFrame(rows).to_csv(HERE/'filter_result_comparison.csv',index=False)
pd.concat(audit).to_csv(HERE/'filter_bin_audit.csv',index=False)
print(pd.DataFrame(rows).to_string(index=False))
