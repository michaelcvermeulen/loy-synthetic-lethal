"""Produce matched 100/50-kb preview renderers from the last approved preview layout."""
from pathlib import Path
import argparse,hashlib,json,shutil
import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
PILOT=ROOT.parent/'WGS/reanalysis/Y_PAR_masked_20261027'
ap=argparse.ArgumentParser();ap.add_argument('size',type=int,choices=[100,50]);args=ap.parse_args()
size=args.size;minimum=50 if size==100 else 25
variant=f'{size}kb'
W=HERE/variant;W.mkdir(exist_ok=True)
OUT=ROOT/'figures/chapter3/drafts/a549_dense_bins_20261027'/variant
CODE=OUT/'code';CODE.mkdir(parents=True,exist_ok=True)
source=ROOT/'figures/chapter3/drafts/a549_xpar_relaxed_20261027/code/make_a549_xpar_preview.py'
original=source.read_text(encoding='utf-8')
shutil.copy2(source,CODE/'previous_preview_renderer.py')
binsfile=ROOT/'research/model_development/a549_par_track_relaxed_20261027/map30_no_difficulty/bins.csv' if size==100 else HERE/'bins_50kb.csv'
covfile=PILOT/('a549_panel_tracks' if size==100 else 'a549_50kb_comparison')/'binned_coverage.csv'
d=pd.read_csv(binsfile)
d['chromosome']='chr'+d.chromosome.astype(str)
d['label']=d['sample'].map({'A549-CTRL2':'sgCTRL','A549-LOY':'LOY'})
d['midpoint']=(d.start+d.end)/2e6
d['log2']=np.log2(d.copynumber.clip(lower=.001))
d['log2_seg']=np.log2(d.segmented.clip(lower=.001))
d.to_csv(W/'qdnaseq_bins.csv',index=False)
shutil.copy2(covfile,W/'binned_coverage.csv')
shutil.copy2(PILOT/'a549_pair_comparison/autosomal_normalization.csv',W/'autosomal_normalization.csv')
ctrl=d.query("label=='sgCTRL' and chromosome=='chrY'").merge(pd.read_csv(covfile).query("sample=='A549-CTRL2' and track=='Y'")[['end','read_count']],on='end',validate='one_to_one')
ctrl['region']='excluded'
ctrl.loc[(ctrl.start>2458320)&(ctrl.end<=10565750),'region']='MSY Yp'
ctrl.loc[(ctrl.start>10883085)&(ctrl.end<=27449937),'region']='euchromatic Yq'
keep=ctrl.use&ctrl.copynumber.gt(.001)&ctrl.read_count.ge(minimum)&ctrl.region.ne('excluded')
n=ctrl[keep].groupby('region').size().to_dict()
assert set(n)=={'MSY Yp','euchromatic Yq'}
parbins=len(pd.read_csv(covfile).query("sample=='A549-CTRL2' and track!='Y'"))
text=original
def replace(a,b):
    global text
    assert text.count(a)==1,(a,text.count(a))
    text=text.replace(a,b)
replace("W = Path('A:/Dropbox/LOY_paralogs/projects/SEPT_11_2026_THESIS/research/model_development/a549_par_track_relaxed_20261027')",f'W = Path({W.as_posix()!r})')
replace('assert len(paired) == 625',f'assert len(paired) == {625 if size==100 else 1250}')
replace('paired.control_primary_reads.ge(50)',f'paired.control_primary_reads.ge({minimum})')
replace('assert len(plot_data) == 81',f'assert len(plot_data) == {sum(n.values())}')
replace('(10 bins)',f'({n["MSY Yp"]} bins)')
replace('(71 bins)',f'({n["euchromatic Yq"]} bins)')
replace("label='Y bin'",f"label='{size}-kb Y bin'")
replace('assert px.included.all() and len(px)==28',f'assert px.included.all() and len(px)=={parbins}')
replace("'status': 'Revised preview: X-PAR in panel a, moderately relaxed Y-bin filter; not adopted into draft 2'",
        f"'status': '{size}-kb broad-coverage preview; X-PAR in panel a; not adopted into draft 2'")
replace("'bin_rule': 'bases >=60%; mappability >=30%; difficulty overlap <50%; finite CN and segmentation in both samples; control corrected CN >0.001; >=50 primary MAPQ>=37 control reads per full Y bin; fully inside MSY Yp/euchromatic Yq; exclude centromere and Yq12'",
        f"'bin_rule': 'bases >=60%; mappability >=30%; no blanket difficulty-overlap exclusion; finite CN and segmentation; control CN >0.001; >={minimum} primary MAPQ>=37 control reads per {size}-kb Y bin; fully inside MSY Yp/euchromatic Yq; centromere/Yq12 excluded'")
replace("'mappability_operational_definition': 'The same annotation scores with relaxed map>=30/difficulty<50 thresholds; QDNAseq refitted. X-PAR assessed from observed primary-read coverage, not the original duplicated-reference mappability score.'",
        f"'mappability_operational_definition': 'Reference mappability >=30% at {size} kb. Difficulty overlaps retained as annotations rather than blanket exclusions. Read support required in control; no LOY abundance filter.',\n    'bin_size_kb': {size},\n    'minimum_control_reads': {minimum}")
(CODE/'make_a549_dense_preview.py').write_text(text,encoding='utf-8')
compile(text,str(CODE/'make_a549_dense_preview.py'),'exec')
begin='# b: The identical coordinate transform is applied to cytobands and every locus.'
end='assert np.allclose(ax_a.get_xlim(), ax_b.get_xlim())'
assert original.split(begin,1)[1].split(end,1)[0]==text.split(begin,1)[1].split(end,1)[0]
(CODE/'derivation.json').write_text(json.dumps({'source':str(source),'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
    'bin_size_kb':size,'minimum_control_reads':minimum,'Y_region_bins':n,'PAR_bins':parbins,
    'RNA_and_ideogram_rendering_unchanged':True},indent=2),encoding='utf-8')
print(CODE/'make_a549_dense_preview.py');print(n)
