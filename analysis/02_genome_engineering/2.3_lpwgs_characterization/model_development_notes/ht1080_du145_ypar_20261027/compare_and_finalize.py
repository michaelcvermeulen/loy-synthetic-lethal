"""Write the per-line comparison, caption and validation after rendering."""
from pathlib import Path
import argparse,hashlib,json,shutil
import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
PILOT=ROOT.parent/'WGS/reanalysis/Y_PAR_masked_20261027'
ap=argparse.ArgumentParser();ap.add_argument('cell',choices=['HT1080','DU145']);args=ap.parse_args()
cell=args.cell;slug=cell.lower()
samples={'HT1080':['HT1080-sgCTRL-D0','HT1080-LOY-D0'],'DU145':['LNCaP-sgCTRL','LNCaP-LOY-2X']}[cell]
OUT=ROOT/'figures/chapter3/drafts/ht1080_du145_ypar_20261027'/slug
SNAP=OUT/'source_data/regional_y_atlases_normalised'/slug
new=pd.read_csv(SNAP/'paired_Y_bins_and_exclusions.csv')
summary=pd.read_csv(SNAP/'regional_median_paired_ratios.csv').set_index('region')
par=pd.read_csv(SNAP/'X_PAR_regional_medians.csv').set_index('track')
old=pd.read_csv(HERE/cell/'original/bins.csv').query("chromosome=='Y'").pivot(index='bin_id',columns='sample',values='copynumber')
rows=[]
for region in summary.index:
    ids=new.loc[new.included & new.region.eq(region),'bin_id']
    ctrl=old.loc[ids,samples[0]];loy=old.loc[ids,samples[1]]
    assert ctrl.gt(.001).all() and loy.ge(-1e-12).all()
    values=100*loy.clip(lower=0)/ctrl
    rows.append(dict(region=region,n_bins=len(ids),old_alignment_same_bins_percent=float(values.median()),
        new_alignment_percent=float(summary.loc[region,'median_percent_of_sgCTRL'])))
pd.DataFrame(rows).to_csv(HERE/cell/'regional_comparison.csv',index=False)

rawold=pd.read_csv(HERE/cell/'original/raw_counts.csv').set_index('bin_id')
rawnew=pd.read_csv(HERE/cell/'masked/raw_counts.csv').set_index('bin_id')
ann=pd.read_csv(ROOT/'research/model_development/wgs_filter_rerun_20261027/bin_annotation.csv')
autobins=ann.loc[ann['use'] & ann.chromosome.isin([str(i) for i in range(1,23)]),'bin_id']
assert len(autobins)==21446
counts=[]
for sample in samples:
    a=rawold.loc[autobins,sample];b=rawnew.loc[autobins,sample]
    counts.append(dict(sample=sample,autosomal_bins=len(autobins),raw_count_correlation=float(a.corr(b)),
        identical_bins=int((a==b).sum()),median_absolute_change=float((a-b).abs().median()),maximum_absolute_change=int((a-b).abs().max())))
pd.DataFrame(counts).to_csv(HERE/cell/'autosomal_count_comparison.csv',index=False)

combined=[]
for mode in ['original','masked']:
    folder=PILOT/'atlas_pair_coverage'/cell/mode
    z=pd.read_csv(folder/'regional_PAR_coverage.csv')
    norm=pd.read_csv(folder/'autosomal_normalization.csv').set_index('sample').median_autosomal_depth
    for region in ['PAR1','PAR2']:
        result={}
        reads={}
        for sample in samples:
            v=z[z['sample'].eq(sample)&z.region.eq(region)]
            x=v[v.chromosome.eq('chrX')].iloc[0]
            depth=float((v.mean_depth*(v.end-v.start0)).sum()/(x.end-x.start0))
            result[sample]=depth/norm[sample]
            reads[sample]=int(v.read_count.sum())
        combined.append(dict(alignment=mode,region=region,control_relative_PAR_depth=result[samples[0]],
            LOY_relative_PAR_depth=result[samples[1]],combined_PAR_percent_of_control=100*result[samples[1]]/result[samples[0]],
            control_high_quality_read_count=reads[samples[0]],LOY_high_quality_read_count=reads[samples[1]]))
pd.DataFrame(combined).to_csv(HERE/cell/'combined_PAR_comparison.csv',index=False)

readchecks=[]
for sample in samples:
    q=PILOT/'qc'
    assert (q/f'{sample}.REALIGNMENT_COMPLETE.txt').exists()
    oldq=json.loads((q/f'{sample}.original_read_qc.json').read_text())
    newq=json.loads((q/f'{sample}.masked_read_qc.json').read_text())
    fq=json.loads((q/f'{sample}.fastq_validation.json').read_text())
    assert oldq['read_fingerprint']==newq['read_fingerprint'] and fq['matches_original_primary_reads']
    readchecks.append(dict(sample=sample,pairs=fq['pairs'],primary_read_ends=newq['counts']['primary'],
        primary_mapping_percent=100*(1-newq['counts']['primary_unmapped']/newq['counts']['primary']),
        read_content_preserved=True,original_sha256=(q/f'{sample}.original.sha256').read_text().split()[0]))
pd.DataFrame(readchecks).to_csv(HERE/cell/'read_validation.csv',index=False)
validation=json.loads((SNAP/'validation.json').read_text(encoding='utf-8'))
assert validation['RNA_values_and_means_unchanged'] and validation['no_points_clipped']
assert new.loc[new.included,'control_primary_reads'].ge(50).all()
assert new.loc[new.included,'reference_bases'].ge(60).all()
assert new.loc[new.included,'reference_mappability'].ge(30).all()
support=[]
for region,v in new[new.included].groupby('region',sort=False):
    support.append(dict(region=region,bins=len(v),minimum_control_primary_reads=int(v.control_primary_reads.min()),
        median_control_primary_reads=float(v.control_primary_reads.median()),
        bins_also_passing_original_strict_mask=int(v.original_strict_mask.sum()),
        bins_with_difficulty_overlap_ge20=int(v.reference_blacklist.ge(20).sum()),
        minimum_reference_mappability=float(v.reference_mappability.min()),
        maximum_percent_of_control=float(v.percent_of_sgCTRL.max())))
pd.DataFrame(support).to_csv(HERE/cell/'display_bin_support.csv',index=False)
expected_rna=132 if cell=='HT1080' else 66
assert validation['n_RNA_observations']==expected_rna
before=json.loads((HERE/'selected_before.json').read_text(encoding='utf-8-sig'))
assert all(hashlib.sha256(Path(r['path']).read_bytes()).hexdigest().lower()==r['sha256'].lower() for r in before)
for name in ['make_atlas_preview.py','panel_a.py.txt','template_snapshot.py','renderer_derivation.json']:
    dest=OUT/'code';dest.mkdir(exist_ok=True);shutil.copy2(HERE/name,dest/name)
fig_num='3.10' if cell=='HT1080' else '3.14'
sampling=('Six RNA samples per condition are displayed: three from August 2025 (circles) and three from May 2026 (triangles).' if cell=='HT1080' else
          'Three RNA samples per condition are displayed. NA marks NLGN4Y, whose control mean expression is zero.')
split=('DU145 retains the pre-existing 4.8-Mb distal/proximal Yp split; this is a fixed comparison boundary, not a newly inferred breakpoint.' if cell=='DU145' else '')
caption_split='Distal and proximal Yp are separated at 4.8 Mb.' if cell=='DU145' else ''
caption=f'''**Fig {fig_num} preview. Regional Y depletion and X-mapped PAR coverage in {cell}.**

**a,** Matched lpWGS after realignment to CHM13v2.0/HG002-Y with the Y copies of PAR masked. Grey points show corrected 100-kb male-specific Y bins as percentages of sgCTRL; red lines show regional medians. Bins require non-N sequence >=60%, mappability >=30% and >=50 primary MAPQ >=37 control reads. Centromere and Yq12 bins are excluded; Yq12 is compressed. {caption_split}

Teal points and lines show X-mapped PAR bins and median ratios, using autosome-normalised primary-read depth. This is combined X/Y signal, not Y-specific PAR retention. X-relative positions are projected into the labelled PAR spans. Dashed lines mark 100% of sgCTRL. One DNA sample was sequenced per population.

**b,** HG002 T2T-Y ideogram, RNA-seq genes and ddPCR target-gene spans, aligned with the male-specific coordinates in a.

**c,** RNA expression relative to mean sgCTRL expression; points show samples and horizontal bars show group means. {sampling} Blue denotes sgCTRL and red denotes Y-depleted cells. Teal identifies PAR genes, whose expression is not Y-specific. RNA samples derive from one established population per condition.
'''
(OUT/'preview_caption.md').write_text(caption,encoding='utf-8')
table=['| Region | Bins | Original alignment, same bins | New alignment |','|---|---:|---:|---:|']
for r in rows:table.append(f"| {r['region']} | {r['n_bins']} | {r['old_alignment_same_bins_percent']:.2f}% | {r['new_alignment_percent']:.2f}% |")
report=f'''# {cell}: Y-PAR-masked realignment and atlas preview

Both authorised samples are realigned and validated. Preview exports use the line-specific Fig {fig_num} template, with the broader 100-kb Y-bin selection and X-PAR tracks requested after the A549 comparison. Selected exports and draft 2 prose/captions are unchanged.

## Regional comparison

{chr(10).join(table)}

These columns use the same bins and the same broader QDNAseq analysis on the original and masked alignments. They isolate the alignment update from the earlier bin-filter change. The DU145 original baseline is refitted from its two atlas samples, rather than silently compared with a different treatment cohort. {split}

## Validation and interpretation

- QDNAseq GC/mappability correction, normalisation, smoothing, segmentation and segmented normalisation were rerun on the matched pair. The broad reference mask uses non-N sequence >=60% and mappability >=30%, excluding PAR from the correction fit. Displayed MSY bins must lie wholly within the specified regions, have finite estimates in both samples, control corrected CN >0.001, and >=50 primary MAPQ >=37 control reads. Centromere and Yq12 are excluded. No minimum LOY count or effect-direction filter is applied, and no positive numerator floor is imposed. Zero/near-zero bins do not establish absolute DNA absence.
- PAR primary-read depth excludes secondary, supplementary, QC-fail and marked-duplicate records, uses MAPQ >=37, and is normalised by the median depth of the same 21,446 strict-mask autosomal bins in each sample. PAR bins start at the X-PAR boundaries; shorter terminal bins are retained with >=20 control reads and a density equivalent to >=50 reads/100 kb. These tracks do not reuse the old unmasked-reference PAR mappability annotation or apply QDNAseq PAR correction.
- Original BAM checksums are unchanged. Every primary read name, mate, sequence and quality is preserved through recovery and realignment. Recovered FASTQs contain the original BAM's processed reads, not raw sequencer FASTQs.
- Same verified reference/index and BWA-MEM2 2.2.1 algorithm, with fixed `-K 160000000`; five mapping workers per sample allow two concurrent jobs within the available CPU/memory budget (A549 used eight workers). Sorting/indexing uses samtools 1.23.1. No index rebuild was needed.
- Primary-read fingerprints, paired-FASTQ checks, BAM integrity and QDNAseq masks were verified. All {expected_rna} RNA observations and the ideogram rendering code are unchanged. The line's own source tables are used throughout.
- Bin eligibility depends on reference properties and control support, not LOY effect direction. Difficulty-region overlaps remain in annotations; the broader display does not imply every complex sequence is fully resolved.
- PAR tracks report X-mapped combined X/Y coverage, not Y-PAR retention. `combined_PAR_comparison.csv` gives whole-region old/new diagnostics; the plotted lines are medians of bin ratios, a different summary.
- `autosomal_count_comparison.csv`, `read_validation.csv`, the figure source-data validation, and `mask_checkpoints.csv` contain the checks. A reference-threshold or bin-size choice is not a cell-fraction estimator.

## Outputs

PDF/PNG/SVG, caption, code and source snapshots: `{OUT.as_posix()}`.

The four-sample alignment driver is `WGS/reanalysis/Y_PAR_masked_20261027/run_atlas_four.sh`. Fit the original/masked pair with `fit_pairs.R {cell} original|masked`; run `atlas_pair_coverage.py {cell} original|masked` under WSL; then `prepare_inputs.py {cell}`, `build_renderer.py`, `make_atlas_preview.py {cell} --output-dir <preview folder>`, and `compare_and_finalize.py {cell}` from the projects directory.
'''
(HERE/cell/'README.md').write_text(report,encoding='utf-8')
(OUT/'README.md').write_text(f'Preview only.\n\nAnalysis and validation: {HERE/cell/"README.md"}\n\nCaption: preview_caption.md\n',encoding='utf-8')
(HERE/cell/'final_validation.json').write_text(json.dumps({'cell':cell,'read_checks':readchecks,
    'Y_bins':int(new.included.sum()),'RNA_observations':expected_rna,'RNA_and_ideogram_unchanged':True,
    'selected_exports_unchanged':True,'regions':rows,'autosomal_comparison':counts,'display_bin_support':support},indent=2),encoding='utf-8')
print('\n'.join(table));print(par.to_string());print('Validated',cell)
