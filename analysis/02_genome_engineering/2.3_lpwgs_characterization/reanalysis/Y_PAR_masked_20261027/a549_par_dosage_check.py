"""Check local A549 X/Y dosage without assuming the published stock karyotype."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import csv,json,statistics,subprocess,hashlib

ROOT=Path(__file__).resolve().parent
THESIS=ROOT.parents[2]/'SEPT_11_2026_THESIS'
OUT=ROOT/'a549_par_dosage_20261028';OUT.mkdir(exist_ok=True)
SAMTOOLS=str(ROOT/'runtime/env/bin/samtools')
ANNOTATION=THESIS/'research/model_development/wgs_filter_rerun_20261027/bin_annotation.csv'
FIT=THESIS/'research/model_development/a549_par_track_relaxed_20261027/map30_no_difficulty/bins.csv'
def read(path):
    with path.open(newline='',encoding='utf-8-sig') as f:return list(csv.DictReader(f))
def save(path,data):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
annotation=read(ANNOTATION)
selected=[]
for r in annotation:
    if r['use'].upper()!='TRUE':continue
    c=r['chromosome'];s=int(float(r['start']))-1;e=int(float(r['end']))
    if c=='X' and s>=2394410 and e<=153925834:region='X_nonPAR'
    elif c=='Y' and s>=2458320 and e<=10565750:region='MSY_Yp'
    elif c=='Y' and s>=10883085 and e<=27449937:region='euchromatic_Yq'
    else:continue
    selected.append({'chromosome':'chr'+c,'start0':s,'end':e,'region':region,'bin_id':r['bin_id']})
bed=OUT/'fixed_strict_nonPAR_XY_bins.bed'
bed.write_text(''.join(f"{r['chromosome']}\t{r['start0']}\t{r['end']}\t{r['region']}:{r['bin_id']}\n" for r in selected))
normal={r['sample']:float(r['median_autosomal_depth']) for r in read(ROOT/'a549_pair_comparison/autosomal_normalization.csv') if r['alignment']=='Y_PAR_masked'}
fit={(r['sample'],r['bin_id']):float(r['copynumber']) for r in read(FIT) if r['copynumber'] and r['use'].upper()=='TRUE'}
par=[r for r in read(ROOT/'a549_pair_comparison/combined_PAR_coverage.csv') if r['alignment']=='Y_PAR_masked']

def process(sample):
    bam=ROOT/'alignment'/f'{sample}.YPARmasked.bam'
    dest=OUT/f'{sample}.nonPAR_XY.bedcov.tsv'
    with dest.open('w') as f:
        subprocess.run([SAMTOOLS,'bedcov','-j','-Q','37','-G','0x800',str(bed),str(bam)],stdout=f,check=True)
    data=[]
    for r,line in zip(selected,dest.read_text().splitlines(),strict=True):
        c,s,e,name,total=line.split('\t');assert (c,int(s),int(e))==(r['chromosome'],r['start0'],r['end'])
        depth=int(total)/(int(e)-int(s))
        data.append(dict(sample=sample,**r,primary_depth=depth,depth_relative_to_autosomal_median=depth/normal[sample],corrected_CN=fit[(sample,r['bin_id'])]))
    summary=[]
    for region in ['X_nonPAR','MSY_Yp','euchromatic_Yq']:
        z=[r for r in data if r['region']==region]
        summary.append({'sample':sample,'region':region,'n_fixed_bins':len(z),
            'median_primary_depth':statistics.median(r['primary_depth'] for r in z),
            'autosomal_median_primary_depth':normal[sample],
            'median_relative_primary_depth':statistics.median(r['depth_relative_to_autosomal_median'] for r in z),
            'median_corrected_CN_same_bins':statistics.median(r['corrected_CN'] for r in z)})
    print(sample,summary,flush=True)
    return data,summary

with ThreadPoolExecutor(max_workers=2) as pool:result=list(pool.map(process,['A549-CTRL2','A549-LOY']))
save(OUT/'per_bin_dosage.csv',[r for data,summary in result for r in data])
summary=[r for data,s in result for r in s]
save(OUT/'nonPAR_XY_dosage.csv',summary)
save(OUT/'PAR_dosage.csv',par)
ctrl={r['region']:r for r in summary if r['sample']=='A549-CTRL2'}
pred=[]
for y_region in ['MSY_Yp','euchromatic_Yq']:
    x=ctrl['X_nonPAR'];y=ctrl[y_region]
    pred.append({'control_Y_region':y_region,
        'raw_X_to_Y_depth_ratio':x['median_relative_primary_depth']/y['median_relative_primary_depth'],
        'predicted_PAR_retention_X_over_XplusY_raw':x['median_relative_primary_depth']/(x['median_relative_primary_depth']+y['median_relative_primary_depth']),
        'corrected_X_to_Y_CN_ratio':x['median_corrected_CN_same_bins']/y['median_corrected_CN_same_bins'],
        'predicted_PAR_retention_X_over_XplusY_corrected':x['median_corrected_CN_same_bins']/(x['median_corrected_CN_same_bins']+y['median_corrected_CN_same_bins'])})
save(OUT/'dosage_model_predictions.csv',pred)
(OUT/'provenance.json').write_text(json.dumps({'samples':['A549-CTRL2','A549-LOY'],
    'fixed_bin_rule':'Original strict reference mask; exclude PAR and region-boundary overlaps; no sample-dependent selection for this check',
    'read_rule':'MAPQ >=37; unmapped, secondary, supplementary, QC-failed and marked duplicate records excluded',
    'normalizer':'Same 21446 strict autosomal bins used in PAR plots',
    'raw_depth_is_not_GC_corrected':True,'corrected_CN_source':str(FIT),
    'model_assumptions':'Representative nonPAR X/Y dosage reflects the respective PAR contributions; X dosage remains unchanged; Y contribution is removed. Exact karyotype and rearrangements are not identified by this check.',
    'annotation_sha256':hashlib.sha256(ANNOTATION.read_bytes()).hexdigest(),
    'fit_sha256':hashlib.sha256(FIT.read_bytes()).hexdigest()},indent=2))
print('Predictions',pred,flush=True)
