from pathlib import Path
import json
import pandas as pd

out=Path(__file__).resolve().parent
summary=pd.read_csv(out/'regional_comparison.csv')
genome=pd.read_csv(out/'genomewide_comparison.csv')
repro=pd.read_csv(out/'baseline_reproduction_checks.csv')
checks=pd.read_csv(out/'mask_checkpoints.csv')
validation=json.loads((out/'comparison_validation.json').read_text())
assert checks.loc[checks['mode'].eq('corrected'),'usable_bins'].eq(22682).all()
assert checks.loc[checks['mode'].eq('corrected'),'usable_Y_bins'].eq(58).all()
lines=['# WGS filter correction: review before draft-2 adoption','',
'The three WGS analyses were rerun from cached raw 100-kb read counts. **Realignment was not required for this bin-filter correction.** The source notebooks are fixed; this task made no edits to draft 2, selected figures, captions or historical analysis objects. Corrected outputs are staged in this directory. A separate later change to the Results Markdown is recorded in the final file audit below.', '',
'## What was fixed','',
'The notebooks previously called `applyFilters()` with all optional filters disabled. In QDNAseq this replaces the existing `use` mask; it re-enabled 8,503 bins that the custom T2T annotation had rejected. The corrected helper preserves the saved annotation mask exactly: non-N bases ≥60%, mappability ≥55%, difficulty-region/blacklist overlap <20%. It refreshes used-read totals and expected variance and checks annotation identity. No thresholds were tuned against the outcome.', '',
'Usable bins fall from **31,185 to 22,682 genome-wide**, and from **625 to 58 on Y**. The remaining Y bins are seven in male-specific Yp and 51 in euchromatic Yq. All PAR1/PAR2 bins fail the original mask. Missing values remain missing; excluded regions are not assigned zero coverage.', '',
'## Regional retention: current figure versus full rerun','',
'Values are the median of per-bin Y-depleted/sgCTRL corrected-copy-number ratios, expressed as percentages. These are relative DNA measurements, **not percentages of Y-negative cells**. “Mask only” is included to distinguish removing bad bins from refitting the entire analysis.', '',
'| Model | Region | Current % | Mask only % | Full rerun % | Usable bins, current → rerun |',
'|---|---|---:|---:|---:|---:|']
for cell in ('HT1080','A549','DU145'):
    sub=summary[summary.cell.eq(cell)]
    for region in sub.region.unique():
        d=sub[sub.region.eq(region)].set_index('analysis')
        def fmt(mode):
            x=d.loc[mode,'median_percent']
            return 'Not quantifiable' if pd.isna(x) else f'{x:.2f}'
        lines.append(f"| {cell} | {region} | {fmt('current')} | {fmt('mask_only')} | {fmt('full_rerun')} | {int(d.loc['current','n_bins'])} → {int(d.loc['full_rerun','n_bins'])} |")
lines += ['',
'The comparison retains the current figure’s 0.001 copy-number floor, control-denominator eligibility rule, and exclusion of bins crossing regional boundaries. The full table also reports ratios without the depleted-sample floor and counts of floor-limited bins. **HT1080 Yq is floor-limited: 33/51 corrected bins are at or below this floor, and the median without the floor is zero. Its 0.23% is a plotting/comparison value, not a measured positive residual fraction.** PAR2 is omitted from the zoomed comparison plot, and neither PAR region receives a corrected retention estimate.', '',
'The broad pattern persists: strong Y depletion on both measured arms in A549; appreciable residual male-specific Yp signal with very low euchromatic Yq signal in HT1080; and residual distal Yp with much lower proximal Yp/Yq signal in DU145. Most regional median changes arise from restricting the bin set; refitting shifts the mask-only medians much less.', '',
'![Before and after, regional Y](../../../figures/chapter3/drafts/wgs_filter_rerun_20261027/ch03_supp_wgs_filter_before_after.png)', '',
'Grey shows the current analysis; purple diamonds show the corrected rerun. Left: all plotted usable Y bins in the 0–27.5-Mb view. Right: regional medians; n denotes bins, not biological replicates. Both ends of each connecting line are shown, with identical values appearing superimposed. The old PAR points remain visible for comparison, but are excluded from the corrected analysis.', '',
'## Genome-wide comparison','',
'This check uses common usable autosomal bins, so differences reflect refitting and normalization rather than simply dropping different locations. CN differences are depleted minus control on the QDNAseq median-normalized linear scale; they are not absolute copies per cell.', '',
'| Model | Common autosomal bins | Correlation of old/new bin CN differences | Median absolute change in CN difference | 95th percentile absolute change |',
'|---|---:|---:|---:|---:|']
for row in genome[genome.scope.eq('autosomes')].itertuples():
    lines.append(f'| {row.cell} | {row.common_bins:,} | {row.old_new_bin_delta_pearson:.4f} | {row.median_abs_change_in_delta_cn:.4f} | {row.p95_abs_change_in_delta_cn:.4f} |')
lines += ['', 'Segmented differences require a separate assessment:', '',
    '| Model | Old/new segmented log2-ratio correlation | Median absolute change | 95th percentile absolute change | Bins changing by >0.1 log2 |',
    '|---|---:|---:|---:|---:|']
for row in genome[genome.scope.eq('autosomes')].itertuples():
    lines.append(f'| {row.cell} | {row.old_new_segment_log2_ratio_pearson:.3f} | {row.median_abs_segment_log2_ratio_change:.3f} | {row.p95_abs_segment_log2_ratio_change:.3f} | {100*row.fraction_segment_log2_ratio_changes_gt_0_1:.2f}% |')
lines += ['',
'In particular, A549’s segmented differences have low correlation despite mostly small absolute changes. Its median absolute change is about 0.025 log2 and 95% of bins change by ≤0.061 log2, but approximately 2.67% change by >0.1 log2. Do not describe the segmented genome as unchanged or reuse old focal/copy-number-class annotations without regeneration. Some changes also occur in an unmodified baseline rerun because CBS is stochastic; the fixed-seed baseline is retained for that distinction.', '',
'The genome-wide diagnostic retains the full observed range of each model. Grey includes the old unfiltered bins; purple contains only corrected usable bins. Chromosome-level summaries are saved separately in `chromosome_comparison.csv`; neither chromosome-averaged stability nor a high correlation guarantees unchanged focal calls.', '',
'![Genome-wide comparison](../../../figures/chapter3/drafts/wgs_filter_rerun_20261027/ch03_supp_wgs_filter_genomewide.png)', '',
'## Reproduction and validation','',
'- Used the original R 4.5.2 and QDNAseq 1.46.0, with DNAcopy 1.84.0 and CGHcall 2.72.0. Full session details are saved.',
'- Started with cached raw counts whose feature metadata exactly matched the original annotation. Ran GC/mappability correction, median normalization, outlier smoothing, square-root CBS segmentation, segmented normalization and the original optional CGHcall step.',
'- Preserved the original cohorts: two HT1080 samples, two A549 samples and eight historical `LNCaP-*` samples (DU145). All passed the original five-million-read threshold. The main DU145 comparison remains `LNCaP-sgCTRL` versus `LNCaP-LOY-2X`.',
'- Repeated the old all-enabled calculation from the same raw counts and software, with a fixed seed (20261027). The table below distinguishes deterministic bin-value reproduction from segmentation variability; the original historical seed was not recorded.', '',
'| Model | Assay | Maximum absolute difference, reproduced versus saved | Correlation |',
'|---|---|---:|---:|']
for row in repro.itertuples():
    lines.append(f'| {row.cohort} | {row.assay} | {row.max_absolute_difference:.6g} | {row.correlation:.6f} |')
lines += ['',
'The raw counts, fitted GC/mappability surfaces and smoothed copy-number values **before segmentation reproduce exactly for every cohort** (`presegmentation_reproduction_checks.csv`). CBS and its subsequent normalization explain the small baseline differences after that point. For the six main samples, these are uniform rescalings: +0.390% in A549 Y-depleted, +0.102% in DU145 control, and none in the other four. After accounting for those scale factors, the bin shapes reproduce to floating-point precision. The regional table includes the fixed-seed baseline as an additional comparison in its CSV.', '',
f'- Verified the intended mask at all {len(checks[checks["mode"].eq("corrected")])} corrected pipeline checkpoints.',
'- Verified the old figure’s regional medians directly from the historical RDS exports, not only from a previously computed summary.',
'- Tested the helper on deliberately reset all-TRUE masks: it restores the original mask and read-count metadata. All three edited notebooks parse, and rejected CN/segment bins are NA in all corrected outputs.',
'- '+(out/'protected_files_check.txt').read_text().strip()+' **at rerun completion**. Hashes cover draft-2 files, selected regional figures, captions, the selected index, historical CN objects and the original bin annotation.',
'- A final audit subsequently detected a separate change to `draftv2/chapter3/chapter3_results.md`; this task did not write to that file and left the newer contents intact. The other 42 protected files, including the selected figure exports and captions, remain unchanged. See `final_protected_file_audit.json` for hashes and timestamps.', '',
'## Does this solve PAR counting?','',
'It solves the lost bin-quality mask. It does not resolve the assignment of homologous PAR reads to X versus Y in the original alignment. With the original quality criteria restored, these data provide no usable PAR bins for a Y-specific retention estimate. Therefore the current PAR DNA percentages should be withdrawn when the figures are updated. The ideogram can still mark PAR locations, and PAR RNA should remain labelled as shared X/Y expression.', '',
'A separately validated PAR-aware analysis would be needed if total PAR dosage is a specific objective; a Y-PAR-masked reference is one possible remapping approach. Masking alone would quantify combined homologous signal, not identify whether each read came from X or Y. It is unnecessary to realign merely to adopt this corrected non-PAR WGS analysis. See the preceding [PAR audit](../par_counting_audit_20261027/README.md) and the [T2T reference definitions](https://github.com/marbl/CHM13#analysis-set).', '',
'## Scope and adoption','',
'This is a staged reanalysis, not a manuscript update. RNA-seq, ATAC-seq and ddPCR have not been rerun. If the corrected WGS is adopted, regenerate WGS-dependent regional/genome-wide figures and DNA–RNA copy-number annotations/integrations, then revise their captions and affected quantitative prose. Do not infer that those downstream numbers are unchanged from these WGS checks alone.', '',
'The main regional patterns can be assessed from the corrected table above, but Yp is sparsely represented: seven bins overall and only two in DU145 distal Yp. These data do not independently pinpoint a breakpoint at 4.8 Mb or provide biological replicate uncertainty.', '',
'## Files and reproduction','',
'- `run_rerun.R`: runs all three cohorts and both corrected/baseline pipelines without replacing historical objects.',
'- `../../../../WGS/code/preserve_t2t_bin_mask.R`: shared correction helper; all three source notebooks now use it and invalidate their correction cache when it changes.',
'- `inputs/`: cached raw-count snapshots and pre-edit notebook backups. Original caches and annotation are read-only inputs.',
'- `HT1080/`, `A549/`, `DU145/`: corrected and baseline CN objects, correction-fit objects, bin tables, sample QC and CGHcall status.',
'- `regional_comparison.csv`, `genomewide_comparison.csv`, `chromosome_comparison.csv`: quantitative comparisons.',
'- `baseline_reproduction_checks.csv`, `mask_checkpoints.csv`, `fix_validation.csv`, `comparison_validation.json`: validation outputs.',
'- `input_manifest.csv`, `protected_files_before.csv`, `protected_files_check.txt`, `sessionInfo.txt`, `run.log`: provenance and run details.',
'- `compare_results.py`: creates the review figures under `figures/chapter3/drafts/wgs_filter_rerun_20261027/`, outside the selected figure folders.', '',
'Run the R scripts from `A:/Dropbox/LOY_paralogs/projects` using `C:/Program Files/R/R-4.5.2/bin/Rscript.exe --vanilla`. For the Python comparison, use the bundled Python with `SEPT_11_2026_THESIS/.tmp/qi_figure_python` on `PYTHONPATH`. The source notebooks will produce their normal output files if explicitly rendered later; this staged rerun does not render them.', '']
(out/'README.md').write_text('\n'.join(lines),encoding='utf-8')
print(out/'README.md')
