"""Record comparisons and a draft-only caption after the new atlas is rendered."""
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
OUT=ROOT/'figures/chapter3/drafts/a549_ypar_masked_20261027'
SNAP=OUT/'source_data/regional_y_atlases_normalised/a549'
PILOT=ROOT.parent/'WGS/reanalysis/Y_PAR_masked_20261027'
source_pairs=SNAP/'paired_Y_bins_and_exclusions.csv'
new=pd.read_csv(source_pairs)
mask_old=pd.read_csv(ROOT/'research/model_development/wgs_filter_rerun_20261027/A549/corrected_bins.csv')
old=mask_old.loc[mask_old.chromosome.eq('Y')].pivot(index='bin_id',columns='sample',values='copynumber')
new_summary=pd.read_csv(SNAP/'regional_median_paired_ratios.csv').set_index('region')
historical=pd.read_csv(ROOT/'research/model_development/wgs_filter_rerun_20261027/regional_comparison.csv')
hist=historical.query("cell == 'A549' and analysis == 'current'").set_index('region')
rows=[]
for region in ['MSY Yp','euchromatic Yq']:
    ids=new.loc[new.included & new.region.eq(region),'bin_id']
    assert old.loc[ids,'A549-LOY'].ge(-1e-12).all()
    rr=100*old.loc[ids,'A549-LOY'].clip(lower=0)/old.loc[ids,'A549-CTRL2']
    assert rr.notna().all()
    rows.append({'region':region,'historical_figure_bins':int(hist.loc[region,'n_bins']),
        'historical_figure_percent':float(hist.loc[region,'median_percent']),
        'filtered_bins':len(ids),'old_alignment_filtered_percent':float(rr.median()),
        'masked_alignment_filtered_percent':float(new_summary.loc[region,'median_percent_of_sgCTRL'])})
comparison=pd.DataFrame(rows);comparison.to_csv(HERE/'regional_comparison.csv',index=False)
par=pd.read_csv(PILOT/'a549_pair_comparison/PAR_ratios.csv')
par.to_csv(SNAP/'PAR_ratios.csv',index=False)
parnew=par.query("alignment == 'Y_PAR_masked'").set_index('region')

before=json.loads((HERE/'selected_figure_before.json').read_text(encoding='utf-8-sig'))
assert all(hashlib.sha256(Path(r['path']).read_bytes()).hexdigest().lower()==r['sha256'].lower() for r in before)
assert (PILOT/'qc/A549-LOY.REALIGNMENT_COMPLETE.txt').exists()
for sample in ['A549-CTRL2','A549-LOY']:
    a=json.loads((PILOT/f'qc/{sample}.original_read_qc.json').read_text())
    b=json.loads((PILOT/f'qc/{sample}.masked_read_qc.json').read_text())
    assert a['read_fingerprint']==b['read_fingerprint']
    f=json.loads((PILOT/f'qc/{sample}.fastq_validation.json').read_text())
    assert f['matches_original_primary_reads']
valid=json.loads((SNAP/'validation.json').read_text(encoding='utf-8'))
assert valid['RNA_values_and_means_unchanged'] and valid['Y_bins_displayed']==58
assert new.loc[new.region.str.startswith('PAR'),'included'].eq(False).all()
assert new.loc[new.included,'percent_of_sgCTRL'].notna().all()
assert valid['no_points_clipped']
comparison_lines=['| Region | Current figure | Corrected filter, old alignment | Corrected filter, masked alignment |',
                  '|---|---:|---:|---:|']
for r in rows:
    comparison_lines.append(f"| {r['region']} | {r['historical_figure_percent']:.2f}% ({r['historical_figure_bins']} bins) | {r['old_alignment_filtered_percent']:.2f}% ({r['filtered_bins']} bins) | {r['masked_alignment_filtered_percent']:.2f}% ({r['filtered_bins']} bins) |")
caption='''**Fig 3.12 draft preview. Regional Y-chromosome depletion and expression in A549 after Y-PAR-masked realignment.**

**a,** lpWGS from the sgCTRL and Y-depleted populations was realigned to the same CHM13v2.0/HG002-Y reference with the Y copies of PAR1 and PAR2 masked. The main track shows the ratio of GC- and mappability-corrected, normalised QDNAseq copy-number estimates in Y-depleted versus sgCTRL cells, as a percentage, in quality-filtered 100-kb bins. Red lines show regional medians (seven male-specific Yp bins; 51 euchromatic Yq bins). The dashed line marks 100% of sgCTRL. The frozen bin mask requires at least 60% non-N bases, at least 55% mappability and less than 20% difficulty-region overlap. Bins overlapping region boundaries, the centromere or Yq12, and bins without a valid control denominator, are excluded. Y-PAR bins are not quantified on this track: their reference copies were masked by design. The Yq12 span is compressed. The inset reports combined X/Y PAR coverage, now mapped to X, separately from male-specific Y retention. Each inset value is total aligned PAR bases per X-PAR base, normalised to the median depth across 21,446 usable autosomal bins, then expressed relative to sgCTRL (MAPQ >=37, primary mapped reads passing QC and duplicate filters). This measurement cannot distinguish X- from Y-derived PAR DNA. One DNA sample was sequenced per population. Regional estimates are descriptive and do not estimate the percentage of Y-negative cells; zero/floor-limited bins do not establish absolute DNA absence.

**b,** Chromosome Y ideogram showing RNA-seq genes above and ddPCR target genes below, including AMELY, DDX3Y and KDM5D. Panels a and b share HG002 T2T-Y coordinates and the same Yq12 compression. Marked ddPCR loci represent target-gene spans.

**c,** RNA-seq expression normalised to mean sgCTRL expression. Points represent four samples per condition; horizontal marks show group means. The dashed line marks 100% of sgCTRL. Blue denotes sgCTRL and red denotes Y-depleted cells. Teal shading and labels identify PAR genes, whose expression is not Y-specific. Samples derive from one established population per condition. The ideogram coordinates and all RNA measurements are unchanged from the current figure.
'''
(OUT/'preview_caption.md').write_text(caption,encoding='utf-8')
report=f'''# A549 Fig 3.12 Y-PAR-masked alignment preview

Status: preview only. Selected figure exports and draft 2 prose/captions were not replaced.

Both A549-CTRL2 and A549-LOY were recovered from their complete original BAMs and realigned using the same verified Y-PAR-masked reference, BWA-MEM2 2.2.1 settings and samtools 1.23.1 sorting. Every primary read name, mate, sequence and quality was preserved. These are recovered fastp-processed sequences, not the original sequencer FASTQs. The reusable native-WSL index was not rebuilt.

## Regional DNA comparison

{chr(10).join(comparison_lines)}

The old-filtered and new-filtered columns use the identical 58 bins and the linear ratio of corrected copy-number estimates. The current figure used its historical set of bins and a 0.001 numerator floor. The new linear plot does not introduce that numerator floor; a control estimate >0.001 is still required. Subnormal negative roundoff values (approximately -1.65e-309) are displayed as zero, with raw corrected values retained in the source tables. Thus the first-to-last comparison includes filtering/normalisation changes as well as realignment; the middle-to-last comparison isolates the alignment update with the corrected pipeline held constant.

## PAR interpretation

After masking, the combined PAR1 signal is {parnew.loc['PAR1','combined_PAR_percent_of_control']:.2f}% of control and PAR2 is {parnew.loc['PAR2','combined_PAR_percent_of_control']:.2f}%. These values combine homologous X/Y DNA, now mapped to X; they are not Y-PAR retention percentages. The old figure's Y-only PAR1/PAR2 percentages therefore are not a like-for-like comparison with this inset. All masked Y-PAR intervals have zero direct mapped coverage by design and are excluded rather than plotted as biological loss.

## Analysis and validation

- QDNAseq 1.46.0 counts use the historical settings (paired reads, MAPQ >=37; defaults otherwise, including supplementary alignments); identical settings were retained for comparison. The separate PAR depth diagnostic uses primary mapped reads, excludes secondary/supplementary/QC-fail/marked-duplicate records and excludes deletions/skips from autosomal bed coverage.
- The original annotation remains valid for retained bins because only Y-PAR bases changed, and every bin overlapping those intervals already fails the mask. The original low mappability of duplicated PAR sequence is not repurposed as a quality score for the masked-reference PAR inset.
- GC/mappability correction, normalisation, outlier smoothing, square-root CBS segmentation and segmented normalisation were rerun for the matched pair with seed 20261027. No copy-number calling was required for this descriptive figure.
- Read fingerprints match for both samples; quality masks are preserved at each analysis checkpoint; 58 Y bins and all 88 RNA measurements are verified; selected exports are unchanged.
- See `regional_comparison.csv`, `raw_count_alignment_comparison.csv`, `mask_checkpoints.csv`, `masked_bins.csv`, and the preview's source-data `validation.json`. PAR source tables are in `{PILOT.as_posix()}/a549_pair_comparison`.

## Outputs

Preview PDF/PNG/SVG and a standalone draft caption: `{OUT.as_posix()}`.

Reproduction order, from the projects directory: `recover_align_A549_LOY.sh` (WSL, one-time; refuses overwrite), `count_new_bams.R` for each new BAM (one-time), `fit_new_pair.R`, `compare_a549_pair.py` (WSL), `prepare_preview.py`, generated `make_a549_ypar_preview.py A549 --output-dir <preview folder>`, then `finish_preview_report.py`. The scripts reuse the completed control alignment and never overwrite original BAMs or selected figures.
'''
(HERE/'README.md').write_text(report,encoding='utf-8')
(OUT/'README.md').write_text('Preview only.\n\nAnalysis, comparisons and reproduction notes: '+str(HERE/'README.md')+'\n\nCaption: preview_caption.md\n',encoding='utf-8')
(HERE/'final_validation.json').write_text(json.dumps({'read_content_preserved_both_samples':True,
    'selected_figure_unchanged':True,'RNA_values_unchanged':True,'Y_bins':58,
    'Y_PAR_not_interpreted_as_loss':True,'regions':rows,
    'combined_PAR_ratios':parnew.reset_index().to_dict('records')},indent=2),encoding='utf-8')
print(comparison.to_string(index=False));print(parnew.to_string())
