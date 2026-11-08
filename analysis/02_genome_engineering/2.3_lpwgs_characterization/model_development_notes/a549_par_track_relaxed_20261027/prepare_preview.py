"""A549 preview with X-PAR tracks in a and a moderately relaxed Y-bin mask."""
from pathlib import Path
import hashlib,json,shutil
import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
OUT=ROOT/'figures/chapter3/drafts/a549_xpar_relaxed_20261027'
CODE=OUT/'code';CODE.mkdir(parents=True,exist_ok=True)
PILOT=ROOT.parent/'WGS/reanalysis/Y_PAR_masked_20261027'
source=ROOT/'figures/chapter3/drafts/a549_ypar_masked_20261027/code/make_a549_ypar_preview.py'
original=source.read_text(encoding='utf-8')
shutil.copy2(source,CODE/'previous_preview_renderer.py')
d=pd.read_csv(HERE/'map30_difficulty50/bins.csv')
d['chromosome']='chr'+d.chromosome.astype(str)
d['label']=d['sample'].map({'A549-CTRL2':'sgCTRL','A549-LOY':'LOY'})
d['midpoint']=(d.start+d.end)/2e6
d['log2']=np.log2(d.copynumber.clip(lower=.001))
d['log2_seg']=np.log2(d.segmented.clip(lower=.001))
d.to_csv(HERE/'qdnaseq_bins.csv',index=False)
shutil.copy2(PILOT/'a549_panel_tracks/binned_coverage.csv',HERE/'binned_coverage.csv')
shutil.copy2(PILOT/'a549_pair_comparison/autosomal_normalization.csv',HERE/'autosomal_normalization.csv')
text=original
def replace(a,b):
    global text
    assert text.count(a)==1,(a,text.count(a))
    text=text.replace(a,b)
replace("W = Path('A:/Dropbox/LOY_paralogs/projects/SEPT_11_2026_THESIS/research/model_development/a549_ypar_realign_20261027')",f'W = Path({HERE.as_posix()!r})')
replace("R / 'sample_metadata.csv']","R / 'sample_metadata.csv', W / 'binned_coverage.csv', W / 'autosomal_normalization.csv']")
replace("paired['included'] = (paired.valid_pair & paired.control_above_floor &",
'''coverage=pd.read_csv(SNAP/'binned_coverage.csv')
control_counts=coverage.query("track == 'Y' and sample == 'A549-CTRL2'").set_index('end').read_count
paired['control_primary_reads']=paired.end.map(control_counts)
assert paired.control_primary_reads.notna().all()
paired['included'] = (paired.valid_pair & paired.control_above_floor & paired.control_primary_reads.ge(50) &''')
replace('assert len(plot_data) == 58','assert len(plot_data) == 81')
replace("hatch='///', lw=.4, zorder=0)","hatch=None, lw=0, zorder=0)")
replace("facecolor='#F1F1F1', edgecolor='#D0D0D0'", "facecolor=SHADE, edgecolor=SHADE")
replace("for a, b in [(p1, q0)]:", "for a, b in [(0, q0), (p2, length)]:")
replace("ylabel='Male-specific Y DNA\\n(% of sgCTRL)'", "ylabel='Y-depleted DNA\\n(% of sgCTRL)'")
replace("    ax_a.text(p1/2, 78, 'Y PAR\\nmasked', ha='center', va='center', rotation=90, fontsize=7.0, color='#666666')",'    pass  # X-PAR measurements are drawn below.')
replace('(7 bins)', '(10 bins)')
replace('(51 bins)', '(71 bins)')
replace("(p1/2, 'PAR1', TEAL, 'center')", "(p1/2, 'PAR1 (X)', TEAL, 'center')")
replace("(gx(length)+.13, 'PAR2', TEAL, 'right')", "(gx(length)+.13, 'PAR2 (X)', TEAL, 'right')")
replace("label='Quality-filtered 100-kb bin'", "label='Y bin'")
replace("label='100% of sgCTRL')]", "label='100% of sgCTRL'),\n           Line2D([], [], marker='o', ls='none', color=TEAL, markersize=3, label='X-PAR bin')]")
replace("bbox_to_anchor=(0, -.025), ncol=3,", "bbox_to_anchor=(0, -.025), ncol=4,")
start=text.index('# PAR is displayed separately:')
stop=text.index('# b: The identical coordinate transform',start)
text=text[:start]+'''# X-PAR bins are displayed in the PAR spans of panel a. These are explicitly
# schematic projections of fractional X-PAR positions, not coordinate liftovers.
norm=pd.read_csv(SNAP/'autosomal_normalization.csv').query("alignment == 'Y_PAR_masked'").set_index('sample').median_autosomal_depth
px=coverage.query("track != 'Y'")
keys=['chromosome','start0','end','track']
px=px.query("sample == 'A549-CTRL2'").merge(px.query("sample == 'A549-LOY'"),on=keys,
    suffixes=('_control','_LOY'),validate='one_to_one')
px['width']=px.end-px.start0
px['included']=(px.read_count_control>=20)&(px.read_count_control/px.width*100000>=50)&(px.depth_control>0)
px['percent_of_sgCTRL']=100*(px.depth_LOY/norm['A549-LOY'])/(px.depth_control/norm['A549-CTRL2'])
assert px.included.all() and len(px)==28
assert px.percent_of_sgCTRL.between(0,DNA_TOP).all()
par_summary=[]
for region,x0,x1,y0,y1 in [('PAR1_X',0,2394410,0,p1_bp),('PAR2_X',153925834,154259566,p2_bp,length_bp)]:
    keep=px.track.eq(region)
    xx=px.loc[keep]
    px.loc[keep,'projected_Y_span_midpoint']=(y0+((xx.start0+xx.end)/2-x0)/(x1-x0)*(y1-y0))/1e6
    xx=px.loc[keep]
    value=float(xx.percent_of_sgCTRL.median())
    ax_a.scatter(gx(xx.projected_Y_span_midpoint),xx.percent_of_sgCTRL,s=8,color=TEAL,alpha=.65,linewidth=0,zorder=3)
    ax_a.plot([gx(y0/1e6),gx(y1/1e6)],[value]*2,color=TEAL,lw=2,zorder=4)
    par_summary.append({'region':region,'n_bins':len(xx),'median_percent_of_sgCTRL':value,
                        'minimum_control_reads':int(xx.read_count_control.min())})
    if region=='PAR1_X':
        median_label('PAR1_X',p1/2,112,f'{value:.0f}%',TEAL)
    else:
        ax_a.annotate(f'{value:.0f}%',xy=(gx((y0+y1)/2e6),value),xytext=(gx(length)+.08,112),
            ha='right',va='center',fontsize=7.8,fontweight='semibold',color=TEAL,
            bbox=dict(facecolor='white',edgecolor='none',alpha=.9,pad=1),
            arrowprops=dict(arrowstyle='-',lw=.6,color=TEAL),zorder=7)
px.to_csv(SNAP/'X_PAR_bin_ratios.csv',index=False)
pd.DataFrame(par_summary).to_csv(SNAP/'X_PAR_regional_medians.csv',index=False)
fig.text(left,.617,'PAR: X-mapped combined X/Y signal; X positions projected within the labelled PAR spans.',
         fontsize=6.3,color=TEAL,va='bottom')

'''+text[stop:]
replace("'status': 'Preview only: matched A549 Y-PAR-masked realignment, not adopted into draft 2'",
        "'status': 'Revised preview: X-PAR in panel a, moderately relaxed Y-bin filter; not adopted into draft 2'")
replace("'bin_rule': 'use TRUE and finite raw/segmented values in both samples; control corrected CN > 0.001; fully inside a nonexcluded region; exclude every bin overlapping centromere or Yq12'",
        "'bin_rule': 'bases >=60%; mappability >=30%; difficulty overlap <50%; finite CN and segmentation in both samples; control corrected CN >0.001; >=50 primary MAPQ>=37 control reads per full Y bin; fully inside MSY Yp/euchromatic Yq; exclude centromere and Yq12'")
replace("'mappability_operational_definition': 'Frozen T2T mask: bases >=60%, mappability >=55%, difficulty overlap <20%; preserved through QDNAseq. All Y-PAR bins excluded. PAR inset uses combined X/Y coverage normalized to median usable-autosomal depth.'",
        "'mappability_operational_definition': 'The same annotation scores with relaxed map>=30/difficulty<50 thresholds; QDNAseq refitted. X-PAR assessed from observed primary-read coverage, not the original duplicated-reference mappability score.',\n    'X_PAR_bins': len(px),\n    'X_PAR_regions': par_summary,\n    'X_PAR_coordinate_note': 'X PAR relative positions projected to the corresponding labelled PAR spans; not a liftover; Y coordinates apply to the MSY track and ideogram.'")
(CODE/'make_a549_xpar_preview.py').write_text(text,encoding='utf-8')
compile(text,str(CODE/'make_a549_xpar_preview.py'),'exec')
begin='# b: The identical coordinate transform is applied to cytobands and every locus.'
end='assert np.allclose(ax_a.get_xlim(), ax_b.get_xlim())'
assert original.split(begin,1)[1].split(end,1)[0]==text.split(begin,1)[1].split(end,1)[0]
(CODE/'derivation.json').write_text(json.dumps({'source':str(source),'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
    'RNA_and_ideogram_rendering_unchanged':True,'selected_filter':'map30_difficulty50',
    'filter_choice':'Moderate relaxation retains difficulty-overlap exclusion and >=50 control reads per full 100-kb Y bin; broader alternatives tabulated separately.'},indent=2),encoding='utf-8')
print(CODE/'make_a549_xpar_preview.py')
