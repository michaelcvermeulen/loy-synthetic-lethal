"""Publish the local comparison summary only after both atlas previews pass QC."""
from pathlib import Path
import hashlib,json
import pandas as pd

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
PILOT=ROOT.parent/'WGS/reanalysis/Y_PAR_masked_20261027'
OUT=ROOT/'figures/chapter3/drafts/ht1080_du145_ypar_20261027'
samples=['HT1080-sgCTRL-D0','HT1080-LOY-D0','LNCaP-sgCTRL','LNCaP-LOY-2X']
assert all((PILOT/'qc'/f'{s}.REALIGNMENT_COMPLETE.txt').exists() for s in samples)
tables=[];checks=[];reads=[];exports=[]
for cell in ['HT1080','DU145']:
    slug=cell.lower();folder=OUT/slug;stem=f'ch03_main_{slug}_lpwgs_chrY_atlas'
    v=json.loads((HERE/cell/'final_validation.json').read_text())
    assert v['selected_exports_unchanged'] and v['RNA_and_ideogram_unchanged']
    checks.append(v)
    fv=json.loads((folder/'source_data/regional_y_atlases_normalised'/slug/'validation.json').read_text())
    assert fv['no_points_clipped'] and fv['RNA_values_and_means_unchanged']
    for ext,h in fv['export_sha256'].items():
        p=folder/f'{stem}.{ext}'
        assert hashlib.sha256(p.read_bytes()).hexdigest()==h
        exports.append({'cell':cell,'path':str(p),'sha256':h})
    r=pd.read_csv(HERE/cell/'read_validation.csv');r.insert(0,'cell',cell);reads.append(r)
    a=pd.read_csv(HERE/cell/'regional_comparison.csv');a.insert(0,'cell',cell);tables.append(a)
combined=pd.concat(tables,ignore_index=True)
combined.to_csv(HERE/'all_regions_comparison.csv',index=False)
rd=pd.concat(reads,ignore_index=True)
assert set(rd['sample'])==set(samples) and rd.read_content_preserved.all()
rd.to_csv(HERE/'all_reads_validation.csv',index=False)
before=json.loads((HERE/'selected_before.json').read_text(encoding='utf-8-sig'))
assert all(hashlib.sha256(Path(r['path']).read_bytes()).hexdigest().lower()==r['sha256'].lower() for r in before)
table=['| Cell line | Region | 100-kb bins | Original alignment | New alignment |','|---|---|---:|---:|---:|']
for r in combined.itertuples():
    table.append(f'| {r.cell} | {r.region} | {r.n_bins} | {r.old_alignment_same_bins_percent:.2f}% | {r.new_alignment_percent:.2f}% |')
readtable=['| Cell line | Original BAM sample ID | Read pairs |','|---|---|---:|']
for r in rd.itertuples():readtable.append(f'| {r.cell} | {r.sample} | {r.pairs:,} |')
text=f'''# HT1080 and DU145: completed realignments and atlas previews

The four samples authorised on 27 October 2026 are complete. The six additional DU145 treatment samples and A549-BC14 were not processed. The legacy `LNCaP` filenames identify the DU145 samples used by the existing atlas.

{chr(10).join(readtable)}

Every primary read name, mate, sequence and quality was preserved through BAM-to-FASTQ recovery and realignment. Original BAM checksums are unchanged. BAM integrity, reference chromosome lengths, paired FASTQs, QDNAseq mask preservation, source data and figure exports passed their checks.

## Atlas previews

- [HT1080 Fig 3.10 PDF]({(OUT/'ht1080/ch03_main_ht1080_lpwgs_chrY_atlas.pdf').as_posix()}) — [caption]({(OUT/'ht1080/preview_caption.md').as_posix()}) and [analysis]({(HERE/'HT1080/README.md').as_posix()}).
- [DU145 Fig 3.14 PDF]({(OUT/'du145/ch03_main_du145_lpwgs_chrY_atlas.pdf').as_posix()}) — [caption]({(OUT/'du145/preview_caption.md').as_posix()}) and [analysis]({(HERE/'DU145/README.md').as_posix()}).

PNG and SVG exports, numerical tables, code and validation accompany each PDF. Each figure retains its own ideogram and RNA data, including the pre-existing DU145 4.8-Mb Yp split. Both use the broader 100-kb Y-bin display selected after the A549 comparison, with X-mapped PAR measurements directly in panel a. The selected figures, central captions and draft 2 have not been replaced.

## Effect of realignment on male-specific Y

{chr(10).join(table)}

The old/new columns compare the same displayed bins under the same broad-filter, matched-pair analysis. They separate the alignment update from the earlier expansion of the bin display. Zero medians are descriptive at this sequencing depth and do not establish complete DNA absence; the plots label values below 1% as `<1%`.

## PAR interpretation and methods

Only the Y copies of PAR are masked in the alignment reference; the PAR measurements are shown. Teal tracks use X-mapped combined X/Y signal, normalised to the median primary-read depth of the same 21,446 autosomal bins per sample. They cannot measure Y-specific PAR retention. X-PAR positions are projected into labelled PAR spans, not treated as an exact coordinate liftover.

Male-specific Y bins use the rerun QDNAseq correction with reference non-N sequence >=60%, mappability >=30%, valid estimates in both samples and >=50 primary MAPQ >=37 control reads per 100-kb bin. Generic difficulty overlap remains annotated but does not automatically exclude the bin. Centromere, Yq12 and boundary-crossing bins are excluded. No minimum LOY count or effect-direction filter is used. PAR uses direct primary-read depth because the old PAR mappability annotation does not apply after Y-PAR masking.

The four-sample driver is `../../../../WGS/reanalysis/Y_PAR_masked_20261027/run_atlas_four.sh`; it reuses the verified A549 reference and index. Five BWA-MEM2 workers per sample and a fixed 160-million-base batch permit two concurrent alignments. The `fit_pairs.R`, `prepare_inputs.py`, `build_renderer.py`, `compare_and_finalize.py` and `finish_summary.py` scripts record the analysis and figure steps. Per-cell reports contain reproduction commands and detailed support/coverage diagnostics.
'''
(HERE/'README.md').write_text(text,encoding='utf-8')
(HERE/'completion_validation.json').write_text(json.dumps({'status':'complete','authorised_samples':samples,
    'original_BAM_read_content_preserved':True,'selected_exports_unchanged':True,'exports':exports,
    'total_read_pairs':int(rd.pairs.sum()),'cell_checks':checks},indent=2),encoding='utf-8')
print('All four authorised samples and both atlas previews validated.')
print('\n'.join(table))
