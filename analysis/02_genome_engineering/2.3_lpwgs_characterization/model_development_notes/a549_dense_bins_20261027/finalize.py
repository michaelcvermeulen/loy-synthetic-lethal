from pathlib import Path
import hashlib,json
import pandas as pd

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
FIG=ROOT/'figures/chapter3/drafts/a549_dense_bins_20261027'
comparison=pd.read_csv(HERE/'resolution_comparison.csv')
before=json.loads((ROOT/'research/model_development/a549_ypar_realign_20261027/selected_figure_before.json').read_text(encoding='utf-8-sig'))
assert all(hashlib.sha256(Path(r['path']).read_bytes()).hexdigest().lower()==r['sha256'].lower() for r in before)
resolution=json.loads((HERE/'resolution_validation.json').read_text())
assert resolution['all_100kb_raw_counts_reproduced_exactly_by_50kb_sums']
assert resolution['annotation_validation']['n_mismatches']==0
validations=[]
for size,minimum in [(100,50),(50,25)]:
    out=FIG/f'{size}kb'
    snap=out/'source_data/regional_y_atlases_normalised/a549'
    check=json.loads((snap/'validation.json').read_text(encoding='utf-8'))
    z=pd.read_csv(snap/'paired_Y_bins_and_exclusions.csv').query('included == True')
    s=pd.read_csv(snap/'regional_median_paired_ratios.csv').set_index('region')
    p=pd.read_csv(snap/'X_PAR_regional_medians.csv').set_index('region')
    assert check['RNA_values_and_means_unchanged'] and check['no_points_clipped']
    assert z.control_primary_reads.ge(minimum).all()
    assert not z.region.str.startswith('PAR').any()
    assert check['Y_bins_displayed']==len(z)
    caption=f'''**Fig 3.12 {size}-kb preview. Regional Y depletion and X-mapped PAR coverage in A549.**

**a,** lpWGS in Y-depleted versus sgCTRL populations after two successive CP + Q53 nucleofections, using the completed matched Y-PAR-masked alignments. Grey points show {size}-kb male-specific Y bins; red lines show median bin-level percentages of sgCTRL. This broader display retains reference bins with at least 60% non-N sequence and at least 30% mappability, without the earlier blanket exclusion based on the generic difficulty-region annotation. Difficulty overlaps remain in the source annotations. Displayed bins require at least {minimum} primary MAPQ >=37 control reads, a control corrected copy-number estimate >0.001, finite estimates in both samples, and a complete position within male-specific Yp or euchromatic Yq. No minimum LOY read count or effect-direction criterion is imposed. Centromere and Yq12 bins remain excluded; Yq12 is compressed. QDNAseq GC/mappability correction, normalisation, smoothing, segmentation and segmented normalisation were fitted for this bin size and eligibility mask. The retained track contains {int(s.loc['MSY Yp','n_bins'])} Yp and {int(s.loc['euchromatic Yq','n_bins'])} euchromatic Yq bins. Zero or near-zero values do not establish absolute DNA absence.

Teal points show X-mapped PAR1 ({int(p.loc['PAR1_X','n_bins'])} bins) and PAR2 ({int(p.loc['PAR2_X','n_bins'])} bins); teal lines show their median ratios. X-mapped PAR coverage includes homologous X/Y signal and cannot quantify Y-specific PAR retention. These measurements use primary mapped reads with MAPQ >=37, excluding secondary, supplementary, QC-fail and marked-duplicate records; depth is normalised to the same fixed median across 21,446 usable autosomal 100-kb bins in each sample before computing the percentage of sgCTRL. PAR bins start at the respective X-PAR boundaries, and terminal bins may be shorter. They require at least 20 control reads and a density equivalent to at least 50 reads per 100 kb. Relative X-PAR positions are projected within the labelled PAR spans for display; this is not an exact coordinate liftover. The dashed line marks 100% of sgCTRL. One DNA sample was sequenced per population.

**b,** Unchanged HG002 T2T-Y ideogram showing RNA-seq genes above and ddPCR target-gene spans below. Male-specific Y coordinates align with a; X-PAR measurements in a use the labelled schematic projection described above.

**c,** Unchanged RNA-seq expression normalised to mean sgCTRL expression. Points represent four samples per condition; horizontal marks show group means. Blue denotes sgCTRL and red denotes Y-depleted cells. Teal identifies PAR genes, whose expression is not Y-specific. Samples derive from one established population per condition.
'''
    if size==50:
        caption+='\nThe 50-kb version is a resolution sensitivity display. Half of its retained Yp bins contain no primary LOY reads. Its lower bin-median ratio should not be interpreted as additional biological depletion relative to the 100-kb display; bin width, retained positions and correction fitting differ.\n'
    (out/'preview_caption.md').write_text(caption,encoding='utf-8')
    validations.append({'bin_size_kb':size,'Y_bins':len(z),'X_PAR_bins':check['X_PAR_bins'],
        'RNA_observations':check['n_RNA_observations'],'minimum_control_reads':int(z.control_primary_reads.min()),
        'selected_exports_unchanged':True,'RNA_and_ideogram_unchanged':True})
    (out/'README.md').write_text('Preview only.\n\nAnalysis and resolution comparison: '+str(HERE/'README.md')+'\n\nCaption: preview_caption.md\n',encoding='utf-8')
table=['| Bin size / region | Bins | Span (Mb) | Median control reads | Bins with no LOY reads | Median corrected ratio |',
       '|---|---:|---:|---:|---:|---:|']
for r in comparison.itertuples():
    if str(r.size_kb) not in ['100','50']:continue
    table.append(f'| {r.size_kb} kb / {r.region} | {r.n_bins} | {r.genomic_span_Mb:.2f} | {r.control_reads_median:.0f} | {r.LOY_zero_read_percent:.1f}% | {r.median_corrected_percent:.2f}% |')
report=f'''# A549: broader 100-kb track versus 50-kb bins

Recommendation: use the expanded 100-kb preview. It contains 67 male-specific Yp bins instead of ten, with 100 euchromatic Yq bins, and preserves stronger per-bin support than the 50-kb version. Both complete three-panel previews are available; neither replaces selected exports or draft 2 text.

## What was changed

The generic difficulty-overlap annotation is retained as metadata rather than used as a blanket exclusion. Reference mappability >=30%, non-N sequence >=60%, anatomical exclusion of centromere/Yq12, a valid control denominator and observed high-quality control-read support remain required. No selection uses the LOY count or effect direction. The 100-kb display requires at least 50 primary MAPQ >=37 control reads per Y bin; actual minimum is 54 in Yp and 50 in Yq. This is a deliberately broader descriptive display, not a claim that every complex region is fully resolved.

## Resolution comparison

{chr(10).join(table)}

The broader 100-kb track spans 6.70 Mb of the 8.11-Mb male-specific Yp interval (about 83%). The 50-kb version spans 6.60 Mb, so it adds spatial subdivisions rather than more genomic coverage or sequencing information. The minimum control-read threshold is scaled from 50 to 25 to compare similar read densities. At 50 kb, the median control count falls from 123 to 64 on Yp and half the retained Yp bins contain no primary LOY reads. Keeping a 50-read minimum at 50 kb retains 91 Yp bins but covers only 4.55 Mb; that sensitivity is included in `resolution_comparison.csv`.

The 50-kb bin-median Yp estimate falls to 0.76%, compared with 2.10% at 100 kb. This is a resolution/selection/fitting sensitivity in sparse data, not evidence of additional biological Y loss. The expanded 100-kb version provides a denser track without that increase in zero-read bins.

## Validation

- Exact 50-kb GC/non-N, mappability and difficulty annotations were recomputed from the existing aligned reference and original interval sources, not copied from parent bins.
- Recombining 50-kb annotations reproduces the frozen 100-kb annotation in all 31,156 bins outside masked Y-PAR intervals (tolerance 1e-5 percentage points; see `resolution_validation.json`).
- Summing 50-kb QDNAseq read-start counts reproduces the 100-kb counts exactly in all 31,185 genomic bins and both samples. Primary aligned-base coverage on Y is also exactly conserved when reaggregated.
- The same completed BAMs, MAPQ threshold and QDNAseq read-count settings were used. The 50-kb analysis was fitted independently with seed 20261027. PAR measurements retain their fixed autosomal depth normalisers.
- The ideogram and all 88 RNA observations are unchanged. Both figures were visually inspected; selected exports remain unchanged by checksum.

## Reproduction and files

- `WGS/reanalysis/Y_PAR_masked_20261027/build_50kb_comparison.py` builds 50-kb annotation and primary-read coverage under its `a549_50kb_comparison/` folder, using the existing WSL runtime/reference. No realignment is required.
- `fit_50kb.R` counts and fits the matched pair at 50 kb. The broader 100-kb fit is the existing `a549_par_track_relaxed_20261027/map30_no_difficulty` result.
- Run `prepare_previews.py 100` and `prepare_previews.py 50`, then the generated `make_a549_dense_preview.py A549 --output-dir <variant folder>` for each, followed by `compare_resolution.py` and `finalize.py`.
- Complete PDF/PNG/SVG previews, captions and source snapshots: `{FIG.as_posix()}/100kb` and `/50kb`.
'''
(HERE/'README.md').write_text(report,encoding='utf-8')
(HERE/'final_validation.json').write_text(json.dumps({'previews':validations,'counts_conserved':True,
    'recommended_preview':'100kb','selected_exports_unchanged':True},indent=2),encoding='utf-8')
print('\n'.join(table));print('Both previews and captions validated; selected exports unchanged.')
