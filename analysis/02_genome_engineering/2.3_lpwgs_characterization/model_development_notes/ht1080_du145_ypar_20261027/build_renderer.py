"""Adapt the cell-specific atlas template; preserve its entire ideogram/RNA code."""
from pathlib import Path
import hashlib,json

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
source=HERE/'template_snapshot.py'
original=source.read_text(encoding='utf-8')
text=original
def replace(a,b):
    global text
    assert text.count(a)==1,(a,text.count(a))
    text=text.replace(a,b)
replace('ROOT = Path(__file__).resolve().parents[3]',f'ROOT = Path({ROOT.as_posix()!r})')
replace("DNA_TOP = 260 if CELL == 'DU145' else 156",'DNA_TOP = 156  # Shared starting scale; expand below if any observed point requires it.')
replace("W = FIG / f'source_data/{SLUG}_lpwgs_confirmation'",f"assert CELL in ['HT1080','DU145']\nW = Path({HERE.as_posix()!r}) / CELL / 'plot_inputs'")
replace("R / 'sample_metadata.csv']","R / 'sample_metadata.csv', W / 'binned_coverage.csv', W / 'autosomal_normalization.csv', W / 'bin_annotation.csv', W / 'sample_map.json']")
replace('atol=1e-12)','atol=1e-12, equal_nan=True)')
replace("paired['included'] = (paired.valid_pair & paired.control_above_floor &",'''sample_map=json.loads((SNAP/'sample_map.json').read_text())
coverage=pd.read_csv(SNAP/'binned_coverage.csv')
control_counts=coverage.loc[coverage.track.eq('Y') & coverage['sample'].eq(sample_map['control'])].set_index('end').read_count
paired['control_primary_reads']=paired.end.map(control_counts)
assert paired.control_primary_reads.notna().all()
annotation=pd.read_csv(SNAP/'bin_annotation.csv').set_index('bin_id')
for field in ['bases','mappability','blacklist']:
    paired['reference_'+field]=paired.bin_id.map(annotation[field])
paired['original_strict_mask']=paired.bin_id.map(annotation['use'])
paired['included'] = (paired.valid_pair & paired.control_above_floor & paired.control_primary_reads.ge(50) &''')
replace("100 * np.exp2(paired.loc[ok, 'log2_LOY'] - paired.loc[ok, 'log2_sgCTRL']))",
        "100 * paired.loc[ok,'copynumber_LOY'].clip(lower=0) / paired.loc[ok,'copynumber_sgCTRL'])\nassert paired.loc[ok,'copynumber_LOY'].ge(-1e-12).all()")
replace("assert len(plot_data) == {'HT1080': 237, 'A549': 248, 'DU145': 222}[CELL]",
        "assert len(plot_data)>0 and not plot_data.region.str.startswith('PAR').any()")
replace("assert plot_data.percent_of_sgCTRL.between(0, DNA_TOP).all()",'assert np.isfinite(plot_data.percent_of_sgCTRL).all()')
replace("fig = plt.figure(figsize=(183 / 25.4, 169 / 25.4))",'''norm=pd.read_csv(SNAP/'autosomal_normalization.csv').set_index('sample').median_autosomal_depth
px=coverage.loc[coverage.track.ne('Y')]
keys=['chromosome','start0','end','track']
px=px.loc[px['sample'].eq(sample_map['control'])].merge(px.loc[px['sample'].eq(sample_map['LOY'])],on=keys,
    suffixes=('_control','_LOY'),validate='one_to_one')
px['width']=px.end-px.start0
px['included']=(px.read_count_control>=20)&(px.read_count_control/px.width*100000>=50)&px.depth_control.gt(0)
px['percent_of_sgCTRL']=100*(px.depth_LOY/norm[sample_map['LOY']])/(px.depth_control/norm[sample_map['control']])
assert px.loc[px.included,'percent_of_sgCTRL'].notna().all()
par_summary=px.loc[px.included].groupby('track').agg(n_bins=('start0','size'),
    median_percent_of_sgCTRL=('percent_of_sgCTRL','median'),minimum_control_reads=('read_count_control','min'))
assert set(par_summary.index)=={'PAR1_X','PAR2_X'}
par_summary.to_csv(SNAP/'X_PAR_regional_medians.csv')
maximum=max(plot_data.percent_of_sgCTRL.max(),px.loc[px.included,'percent_of_sgCTRL'].max())
if maximum>DNA_TOP-6: DNA_TOP=float(np.ceil((maximum+12)/50)*50)
assert plot_data.percent_of_sgCTRL.between(0,DNA_TOP).all()
assert px.loc[px.included,'percent_of_sgCTRL'].between(0,DNA_TOP).all()
fig = plt.figure(figsize=(183 / 25.4, 169 / 25.4))''')
start=text.index('# a: One corrected bin-level DNA ratio;')
end=text.index('# b: The identical coordinate transform',start)
text=text[:start]+(HERE/'panel_a.py.txt').read_text(encoding='utf-8')+'\n'+text[end:]
replace("'status': 'Author-approved three-panel layout, adopted 2026-10-27'",
        "'status': 'New-alignment preview on the cell-specific template; broad 100-kb filter; not adopted into draft 2'")
replace("'DNA_formula': '100 * 2**(saved log2_LOY - saved log2_sgCTRL)'",
        "'DNA_formula': '100 * max(corrected_CN_LOY,0)/corrected_CN_control; no positive numerator floor',\n    'X_PAR_formula': '100 * (LOY primary-read depth / LOY autosomal median) / (control primary-read depth / control autosomal median)'")
replace("'bin_rule': 'use TRUE and finite raw/segmented values in both samples; control corrected CN > 0.001; fully inside a nonexcluded region; exclude every bin overlapping centromere or Yq12'",
        "'bin_rule': 'Reference non-N>=60%, mappability>=30%; no blanket difficulty-overlap exclusion; >=50 primary MAPQ>=37 control reads per full Y bin; finite CN/segmentation, control CN>0.001; whole bins within MSY regions; centromere/Yq12 excluded'")
replace("'separate_mappability_mask_available': False","'separate_mappability_mask_available': True")
replace("'mappability_operational_definition': 'source use flag plus anatomical exclusion and valid-control-coverage rules; no new mappability score inferred'",
        "'mappability_operational_definition': 'Original reference annotation with a >=30% threshold; explicit sample control-read support; difficulty overlap retained as metadata',\n    'PAR_coordinate_note': 'X-relative positions projected within labelled PAR spans, not a coordinate liftover',\n    'X_PAR_bins_displayed': int(px.included.sum()),\n    'X_PAR_regions': json.loads(par_summary.reset_index().to_json(orient='records')),\n    'DU145_Yp_split_bp': 4800000 if CELL=='DU145' else None")
begin='# b: The identical coordinate transform is applied to cytobands and every locus.'
finish='assert np.allclose(ax_a.get_xlim(), ax_b.get_xlim())'
assert original.split(begin,1)[1].split(finish,1)[0]==text.split(begin,1)[1].split(finish,1)[0]
replace("ax_c.text(1, -.36, 'Teal shading: PAR genes (not Y-specific)', transform=ax_c.transAxes,\n          ha='right', va='center', color=TEAL, fontsize=6.6)", '')
compile(text,str(HERE/'make_atlas_preview.py'),'exec')
(HERE/'make_atlas_preview.py').write_text(text,encoding='utf-8')
(HERE/'template_snapshot.py').write_text(original,encoding='utf-8')
(HERE/'renderer_derivation.json').write_text(json.dumps({'source':str(source),
    'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'ideogram_RNA_data_and_plotting_unchanged':True,
    'removed_annotations':['Blue explanatory note beneath panel a','Blue explanatory note beneath panel c'],
    'cells':['HT1080','DU145'],'DU145_Yp_split_preserved':4800000},indent=2),encoding='utf-8')
print(HERE/'make_atlas_preview.py')
