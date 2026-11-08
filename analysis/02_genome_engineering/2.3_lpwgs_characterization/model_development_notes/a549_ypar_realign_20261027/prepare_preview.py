"""Build an isolated A549 preview renderer from the adopted atlas layout."""
from pathlib import Path
import hashlib,json,shutil
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[3]
HERE=Path(__file__).resolve().parent
OUT=ROOT/'figures/chapter3/drafts/a549_ypar_masked_20261027'
OUT.mkdir(parents=True,exist_ok=True)
CODE=OUT/'code';CODE.mkdir(exist_ok=True)
source=ROOT/'figures/chapter3/code/make_regional_y_atlas.py'
original=source.read_text(encoding='utf-8')
snapshot=CODE/'adopted_renderer_snapshot.py'
if snapshot.exists():
    assert snapshot.read_text(encoding='utf-8')==original
else:
    shutil.copy2(source,snapshot)

if (HERE/'masked_bins.csv').exists():
    d=pd.read_csv(HERE/'masked_bins.csv')
    d['chromosome']='chr'+d.chromosome.astype(str)
    d['label']=d['sample'].map({'A549-CTRL2':'sgCTRL','A549-LOY':'LOY'})
    d['midpoint']=(d.start+d.end)/2e6
    d['log2']=np.log2(d.copynumber.clip(lower=.001))
    d['log2_seg']=np.log2(d.segmented.clip(lower=.001))
    assert d.loc[~d.use,['copynumber','segmented','log2','log2_seg']].isna().all().all()
    d.to_csv(HERE/'qdnaseq_bins.csv',index=False)
else:
    print('Preparing renderer only; awaiting realigned pair analysis for plot data.')

text=original
def replace(old,new):
    global text
    assert text.count(old)==1,(old,text.count(old))
    text=text.replace(old,new)
replace("ROOT = Path(__file__).resolve().parents[3]",f"ROOT = Path({str(ROOT.as_posix())!r})")
replace("W = FIG / f'source_data/{SLUG}_lpwgs_confirmation'",f"assert CELL == 'A549'\nW = Path({str(HERE.as_posix())!r})")
replace("FLOOR = .001  # Saved numerical log transformation floor, not a detection limit.",
        "FLOOR = .001  # Control-denominator eligibility threshold; no floor added to displayed LOY DNA.")
replace("atol=1e-12)","atol=1e-12, equal_nan=True)")
replace("100 * np.exp2(paired.loc[ok, 'log2_LOY'] - paired.loc[ok, 'log2_sgCTRL']))",
        "100 * paired.loc[ok, 'copynumber_LOY'].clip(lower=0) / paired.loc[ok, 'copynumber_sgCTRL'])\nassert paired.loc[ok,'copynumber_LOY'].ge(-1e-12).all()  # only numerical subnormal negatives may become zero")
replace("assert len(plot_data) == {'HT1080': 237, 'A549': 248, 'DU145': 222}[CELL]",
        "assert len(plot_data) == 58\nassert set(plot_data.region) == {'MSY Yp','euchromatic Yq'}")
replace("for a, b in [(0, q0), (p2, length)]:", "for a, b in [(p1, q0)]:")
replace("ax_a.axvspan(gx(a), gx(b), color=SHADE, zorder=0)",
        "ax_a.axvspan(gx(a), gx(b), facecolor='#F1F1F1', edgecolor='#D0D0D0', hatch='///', lw=.4, zorder=0)")
replace("    median_label('PAR1', p1/2, 78, f'{summary.loc[\"PAR1\", \"median_percent_of_sgCTRL\"]:.0f}%', TEAL)",
        "    ax_a.text(p1/2, 78, 'Y PAR\\nmasked', ha='center', va='center', rotation=90, fontsize=7.0, color='#666666')")
replace("f'{summary.loc[\"MSY Yp\", \"median_percent_of_sgCTRL\"]:.1f}%', COL['LOY'])",
        "f'{summary.loc[\"MSY Yp\", \"median_percent_of_sgCTRL\"]:.1f}% (7 bins)', COL['LOY'])")
replace("median_label('euchromatic Yq', 19.0, 18, f'{summary.loc[\"euchromatic Yq\", \"median_percent_of_sgCTRL\"]:.1f}%', COL['LOY'])",
        "median_label('euchromatic Yq', 19.0, 18, f'{summary.loc[\"euchromatic Yq\", \"median_percent_of_sgCTRL\"]:.1f}% (51 bins)', COL['LOY'])")
replace("label='100-kb bin'", "label='Quality-filtered 100-kb bin'")
replace("ylabel='Y-depleted DNA\\n(% of sgCTRL)'", "ylabel='Male-specific Y DNA\\n(% of sgCTRL)'")
replace("# b: The identical coordinate transform is applied to cytobands and every locus.", '''# PAR is displayed separately: these are combined X/Y measurements, not Y retention.
par_source=ROOT.parent/'WGS/reanalysis/Y_PAR_masked_20261027/a549_pair_comparison/PAR_ratios.csv'
par=pd.read_csv(par_source).query("alignment == 'Y_PAR_masked'").set_index('region')
shutil.copy2(par_source,SNAP/'PAR_ratios.csv')
values=par.loc[['PAR1','PAR2'],'combined_PAR_percent_of_control'].to_numpy()
ax_par=fig.add_axes([.607,.790,.203,.100],facecolor='white')
ax_par.axhline(100,color=COL['sgCTRL'],lw=.7,ls='--')
ax_par.scatter([0,1],values,color=TEAL,s=16,zorder=3)
for i,value in enumerate(values):
    ax_par.text(i,value-10,f'{value:.0f}%',ha='center',va='top',fontsize=6.6,color=TEAL)
ax_par.set(xlim=(-.5,1.5),ylim=(0,max(120,float(values.max())+25)),yticks=[0,50,100],xticks=[0,1],xticklabels=['PAR1','PAR2'])
ax_par.tick_params(labelsize=6,length=2,pad=2)
ax_par.set_ylabel('% of sgCTRL',fontsize=6,labelpad=2)
ax_par.spines[['top','right']].set_visible(False)
ax_par.set_title('Combined X/Y PAR DNA',fontsize=7.0,color=TEAL,pad=5)
fig.text(.032,.992,'A549 | Y-PAR-masked alignment preview',fontsize=7,color='#666666',va='top')

# b: The identical coordinate transform is applied to cytobands and every locus.''')
# b and c remain exactly the adopted rendering code; only a and metadata change.
replace("'status': 'Author-approved three-panel layout, adopted 2026-10-27'",
        "'status': 'Preview only: matched A549 Y-PAR-masked realignment, not adopted into draft 2'")
replace("'DNA_formula': '100 * 2**(saved log2_LOY - saved log2_sgCTRL)'",
        "'DNA_formula': '100 * max(corrected_CN_LOY,0) / corrected_CN_sgCTRL; no positive numerator floor; negative subnormal roundoff treated as zero',\n    'subnormal_negative_values_rounded_to_zero': int((plot_data.copynumber_LOY.lt(0)).sum())")
replace("'separate_mappability_mask_available': False", "'separate_mappability_mask_available': True")
replace("'mappability_operational_definition': 'source use flag plus anatomical exclusion and valid-control-coverage rules; no new mappability score inferred'",
        "'mappability_operational_definition': 'Frozen T2T mask: bases >=60%, mappability >=55%, difficulty overlap <20%; preserved through QDNAseq. All Y-PAR bins excluded. PAR inset uses combined X/Y coverage normalized to median usable-autosomal depth.'")
(CODE/'make_a549_ypar_preview.py').write_text(text,encoding='utf-8')
compile(text,str(CODE/'make_a549_ypar_preview.py'),'exec')
start='# b: The identical coordinate transform is applied to cytobands and every locus.'
end='assert np.allclose(ax_a.get_xlim(), ax_b.get_xlim())'
assert original.split(start,1)[1].split(end,1)[0] == text.split(start,1)[1].split(end,1)[0]
(CODE/'renderer_derivation.json').write_text(json.dumps({'template':str(source),
    'template_sha256':hashlib.sha256(original.encode('utf-8')).hexdigest(),
    'new_renderer_sha256':hashlib.sha256(text.encode('utf-8')).hexdigest(),
    'note':'Isolated preview renderer derived from adopted a-c layout. Panels b/c unchanged. Panel a uses matched masked alignments and corrected quality mask, with a separate combined-PAR depth inset.'},indent=2),encoding='utf-8')
print(CODE/'make_a549_ypar_preview.py')
