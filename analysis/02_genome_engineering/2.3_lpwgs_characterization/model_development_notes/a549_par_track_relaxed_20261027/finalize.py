from pathlib import Path
import hashlib,json
import pandas as pd

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
OUT=ROOT/'figures/chapter3/drafts/a549_xpar_relaxed_20261027'
SNAP=OUT/'source_data/regional_y_atlases_normalised/a549'
checks=json.loads((SNAP/'validation.json').read_text(encoding='utf-8'))
assert checks['Y_bins_displayed']==81 and checks['X_PAR_bins']==28
assert checks['RNA_values_and_means_unchanged'] and checks['no_points_clipped']
yp=pd.read_csv(SNAP/'paired_Y_bins_and_exclusions.csv').query('included == True')
assert yp.control_primary_reads.ge(50).all()
assert yp.percent_of_sgCTRL.between(0,156).all()
px=pd.read_csv(SNAP/'X_PAR_bin_ratios.csv')
assert px.included.all() and px.read_count_control.ge(20).all()
assert (px.read_count_control/px.width*100000).ge(50).all()
before=json.loads((ROOT/'research/model_development/a549_ypar_realign_20261027/selected_figure_before.json').read_text(encoding='utf-8-sig'))
assert all(hashlib.sha256(Path(r['path']).read_bytes()).hexdigest().lower()==r['sha256'].lower() for r in before)
comparison=pd.read_csv(HERE/'filter_result_comparison.csv')
table=['| Filter | Yp bins | Yq bins | Median Yp (% control) | Median Yq (% control) |',
       '|---|---:|---:|---:|---:|']
labels={'strict55_20':'Mappability >=55%, difficulty <20%',
        'map30_difficulty50':'Mappability >=30%, difficulty <50% — chosen',
        'map30_difficulty80':'Mappability >=30%, difficulty <80%',
        'map30_no_difficulty':'Mappability >=30%, no difficulty-overlap exclusion'}
for name in labels:
    z=comparison[comparison['filter'].eq(name)].set_index('region')
    a,b=z.loc['MSY Yp'],z.loc['euchromatic Yq']
    table.append(f"| {labels[name]} | {int(a.n_bins)} | {int(b.n_bins)} | {a.median_percent:.2f} | {b.median_percent:.2f} |")
caption='''**Fig 3.12 revised preview. Regional Y depletion and X-mapped PAR coverage in A549.**

**a,** lpWGS in Y-depleted versus sgCTRL populations after two successive CP + Q53 nucleofections. Both libraries use the previously completed matched Y-PAR-masked alignments; PAR measurements are now displayed directly in the main track. Grey points show male-specific Y bins; red lines show regional medians. QDNAseq GC/mappability correction, normalisation, outlier smoothing, segmentation and segmented normalisation were rerun with a less restrictive annotation filter: at least 60% non-N sequence, at least 30% mappability and less than 50% difficulty-region overlap. Displayed Y bins additionally require at least 50 primary MAPQ >=37 control reads per full 100-kb bin, a valid control copy-number denominator and finite estimates in both samples. Bins must lie wholly within male-specific Yp or euchromatic Yq; the centromere and Yq12 remain excluded, and Yq12 is compressed. The retained track contains ten Yp bins and 71 euchromatic Yq bins. No lower bound on the Y-depleted read count is used.

Teal points show 24 PAR1 and four PAR2 bins measured on chromosome X; teal lines show their median ratios. In these alignments, X-mapped PAR reads include homologous X/Y signal and cannot quantify Y-specific PAR retention. PAR coverage is normalised to each sample's median depth across the same 21,446 usable autosomal bins, then expressed as a percentage of sgCTRL. These measurements use primary mapped reads with MAPQ >=37, excluding secondary, supplementary, QC-fail and marked-duplicate records. PAR bins start at their X-PAR boundaries; terminal bins are shorter than 100 kb. They require at least 20 control reads and a control-read density equivalent to at least 50 reads per 100 kb. The original duplicated-reference PAR mappability annotation is not applied to these X-mapped coverage measurements. Relative positions within X-PAR1 and X-PAR2 are projected into the labelled PAR spans for display; this is a schematic placement, not a genomic coordinate liftover. The dashed line marks 100% of sgCTRL. One DNA sample was sequenced per population. Values are descriptive; zero or near-zero bins do not establish absolute DNA absence.

**b,** Unchanged HG002 T2T-Y ideogram showing RNA-seq genes above and ddPCR target-gene spans below. The male-specific Y track in a and the ideogram share Y coordinates; the PAR tracks in a use the explicitly labelled X-derived measurements described above.

**c,** Unchanged RNA-seq expression normalised to mean sgCTRL expression. Points represent four samples per condition; horizontal marks show group means. Blue denotes sgCTRL and red denotes Y-depleted cells. Teal identifies PAR genes, whose expression is not Y-specific. Samples derive from one established population per condition.
'''
(OUT/'preview_caption.md').write_text(caption,encoding='utf-8')
report=f'''# A549 X-PAR tracks and relaxed bin filter

Revised preview only. Selected figure exports and draft 2 text/captions are unchanged; the prior preview is preserved separately.

## What changed

X-PAR1 and X-PAR2 coverage now appears directly in panel a. The hatched PAR exclusions and separate inset were removed. The completed Y-PAR-masked BAMs were reused: the change is to the displayed PAR measurement, not to the reference or alignments. Centromere/Yq12 exclusions remain. All 88 RNA observations and the ideogram are unchanged.

## Filter comparison and choice

The original difficulty-overlap rule, rather than the 55% mappability threshold alone, caused most exclusions. Lowering mappability to 30% while retaining difficulty overlap <20% still leaves 58 male-specific Y bins. The chosen moderate relaxation keeps mappability >=30%, difficulty overlap <50% and non-N bases >=60%, with at least 50 primary MAPQ >=37 control reads per full 100-kb Y bin. All 81 chosen Y bins actually have at least 56 such control reads. Eligibility is not selected on LOY counts or effect direction.

{chr(10).join(table)}

All rows additionally apply the same anatomical exclusions and control-read support requirement. The chosen filter adds 23 Y bins (58 to 81; ten Yp and 71 Yq) while retaining a difficulty-overlap safeguard. More permissive alternatives were examined as sensitivity analyses, not chosen to achieve a particular retention value. The main strong-depletion conclusion persists; the Yp median is sensitive to which small set of bins is retained. This comparison supports a practical preview threshold, not a formally optimised or independently validated filter.

QDNAseq was refitted separately for each alternative on the same matched raw counts, preserving its specified mask throughout and using seed 20261027. The unchanged BAM read-count settings include the historical QDNAseq defaults; the additional displayed-bin support check and PAR diagnostic use primary reads. PAR bins were excluded from the QDNAseq correction fits because their original mappability scores describe the duplicated reference. The autosomal normaliser for the separate PAR depth ratios is held fixed across these displays.

## PAR values and coordinates

All 28 X-PAR bins pass the observed-control-coverage rule. The bin-median ratios are 66.83% for PAR1 and 63.65% for PAR2. The earlier inset used whole-region aggregate coverage ratios (67.41% and 63.26%); the small change reflects using medians of bin ratios for the new track. These are combined X/Y signals. They are displayed at fractional positions within labelled PAR spans and must not be read as exact X-to-Y coordinate liftovers.

## Files and reproduction

- `refit_filters.R`: staged QDNAseq fits; original RDS objects and global annotation were not changed.
- `compare_filters.py`: compares retained bins and results; `filter_result_comparison.csv` and `filter_bin_audit.csv` preserve all alternatives.
- `prepare_preview.py`: derives an isolated renderer, preserving b/c code exactly.
- X-PAR and primary Y control coverage: `{(ROOT.parent/'WGS/reanalysis/Y_PAR_masked_20261027/a549_panel_tracks').as_posix()}`; produced by `a549_panel_tracks.py` without realignment.
- New PDF/PNG/SVG, caption and source snapshots: `{OUT.as_posix()}`.
- Reproduce from the projects directory: `refit_filters.R`, `compare_filters.py`, `prepare_preview.py`, generated `make_a549_xpar_preview.py A549 --output-dir <new preview folder>`, then `finalize.py`.
'''
(HERE/'README.md').write_text(report,encoding='utf-8')
(OUT/'README.md').write_text('Revised preview only.\n\nAnalysis and filter comparison: '+str(HERE/'README.md')+'\n\nCaption: preview_caption.md\n',encoding='utf-8')
(HERE/'final_validation.json').write_text(json.dumps({'Y_bins':81,'X_PAR_bins':28,'RNA_observations':88,
    'minimum_Y_control_primary_reads':int(yp.control_primary_reads.min()),'selected_exports_unchanged':True,
    'PAR_tracks_in_panel_a':True,'PAR_reference_alignment_unchanged':True,'ideogram_and_RNA_unchanged':True},indent=2),encoding='utf-8')
print('\n'.join(table))
print('Verified 81 Y bins, 28 X-PAR bins, unchanged RNA/ideogram and unchanged selected exports.')
