from pathlib import Path
import json
import numpy as np
import pandas as pd
import requests

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'integration/interactive_atac_rna'
source=(ROOT/'report/interactive_assets/data.js').read_text(encoding='utf-8')
data=json.loads(source[len('window.ATAC_RNA_DATA='):-1])
checks=[]
for key,panel in data['panels'].items():
    rows=panel['rows'];ids=[r[0] for r in rows]
    assert len(ids)==len(set(ids)) and all(i in data['expression'] for i in ids)
    overlap=sum(r[3]<.1 and r[4]<.05 for r in rows)
    dataset,feature=key.split(':')
    if dataset=='Combined':
        old=pd.read_csv(ROOT/'integration/feature_overlaps/overlap_counts.csv')
        expected=old[(old.scheme=='RNA_fixed_005')&(old.feature==feature)&np.isclose(old.ATAC_FDR,.1)].iloc[0]
        assert overlap==expected.overlap and len(rows)==expected.universe
    else:
        old=pd.read_csv(ROOT/'integration/batch_RNA_overlaps/overlap_summary.csv')
        expected=old[(old.batch==dataset)&(old.feature==feature)&(old.universe_type=='batch-specific')].iloc[0]
        assert overlap==expected.Both and len(rows)==expected.eligible_genes
    checks.append({'panel':key,'genes':len(rows),'overlap':overlap})
response=requests.get('http://127.0.0.1:8765/api/locus',params={'gene':'25850','feature':'promoter','view':'best'},timeout=30)
response.raise_for_status();locus=response.json()
old=np.load(ROOT/'integration/gene_atlas/coverage/25850_peak1.npz')
assert locus['bin']==25 and locus['start']<=old['starts'][0] and locus['end']>=old['ends'][-1]
edge=np.r_[locus['starts'],locus['ends'][-1]];v=np.array(locus['coverage'])
maximum_error=0
for i in range(6):
    integral=np.r_[0,np.cumsum(v[i]*np.diff(edge))]
    computed=(np.interp(old['ends'],edge,integral)-np.interp(old['starts'],edge,integral))/(old['ends']-old['starts'])
    maximum_error=max(maximum_error,float(abs(computed-old['values'][i]).max()))
    assert np.allclose(computed,old['values'][i],rtol=2e-6,atol=2e-6)
raw=pd.read_csv(ROOT/'counts/peak_counts.tsv',sep='\t').set_index('peak_id')
norm=pd.read_csv(ROOT/'differential/normalization_factors.csv')
eff=norm[norm.arm=='primary_TMM'].set_index('sample').loc[locus['samples'],'effective_size']
for peak in locus['linked_peaks']:
    assert np.array_equal(peak['raw'],raw.loc[peak['id'],locus['samples']])
    assert np.allclose(peak['cpm'],raw.loc[peak['id'],locus['samples']].values/eff.values*1e6,atol=1e-6)
assert requests.get('http://127.0.0.1:8765/api/locus',params={'gene':'unknown'},timeout=10).status_code==400
assert requests.get('http://127.0.0.1:8765/config/samples.tsv',timeout=10).status_code==404
rna=pd.read_csv(OUT/'RNA_normalization_verification.csv')
assert rna.max_relative_baseMean_difference.max()<1e-8
result={'status':'passed','panels':checks,'RNA_expression_genes':len(data['expression']),
        'on_demand_coverage_vs_original_BigWig_max_absolute_difference':maximum_error,
        'checks':['All nine overlap panels reproduce original eligible universes and counts',
                  'RNA normalized counts reproduce the combined DESeq2 baseMean for all 16,541 input genes',
                  'All six on-demand ATAC tracks match original BigWig-derived atlas coverage within numerical precision',
                  'Peak raw counts and TMM CPM match the original count matrix and normalization',
                  'Invalid gene requests and non-report file requests are rejected']}
(OUT/'data_verification.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result,indent=2))
