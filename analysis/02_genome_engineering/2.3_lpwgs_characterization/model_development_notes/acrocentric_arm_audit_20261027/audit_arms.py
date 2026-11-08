"""Audit the saved arm heatmap against original bin quality and cached read counts."""
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd

OUT=Path(__file__).resolve().parent
THESIS=OUT.parents[2]
PROJECTS=THESIS.parent
RERUN=OUT.parent/'wgs_filter_rerun_20261027'
FIGDATA=THESIS/'figures/chapter3/source_data/du145_arm_recurrence'
targets=['chr13p','chr14p','chr15p','chr21p','chr22p']
comparators=['chr5p','chr8p','chr8q','chr20p']
paths={
    'annotation':RERUN/'bin_annotation.csv',
    'cytobands':PROJECTS/'WGS/chm13v2.0_cytobands_allchrs.bed',
    'original_arm_summary':FIGDATA/'arm_DNA_vs_RNA.csv',
    'original_bins':PROJECTS/'WGS/figures/DU145_recurrence_2026-09-30/wgs_bins.csv',
    'corrected_bins':RERUN/'DU145/corrected_bins.csv',
    'raw_count_object':RERUN/'inputs/DU145_cached_raw_counts.rds',
    'counts_export':OUT/'DU145_raw_read_counts.csv',
}
a=pd.read_csv(paths['annotation'])
a['chromosome']='chr'+a.chromosome.astype(str)
assert a.use.dtype==bool
assert (a.use==(a.bases.ge(60)&a.mappability.ge(55)&a.blacklist.lt(20))).all()
c=pd.read_csv(paths['cytobands'],sep='\t',header=None,names=['chrom','start','end','band','stain'])
arms=['chr'+str(ch)+arm for ch in list(range(1,23))+['X','Y'] for arm in ['p','q']]
a['arm']=None
for chrom in a.chromosome.unique():
    z=c[c.chrom.eq(chrom)]
    mid=(a.start+a.end)/2
    a.loc[a.chromosome.eq(chrom)&mid.le(z.loc[z.band.str.startswith('p'),'end'].max()),'arm']=chrom+'p'
    a.loc[a.chromosome.eq(chrom)&mid.ge(z.loc[z.band.str.startswith('q'),'start'].min()),'arm']=chrom+'q'
counts=pd.read_csv(paths['counts_export'])
samples=list(counts.columns[1:])
a=a.merge(counts,on='bin_id',validate='one_to_one')
assert len(a)==31185 and not a[samples].isna().any().any()
source=pd.read_csv(paths['original_arm_summary'])
original=pd.read_csv(paths['original_bins'],usecols=['sample','bin_id','use'])
assert original.use.all()
fixed=pd.read_csv(paths['corrected_bins'])
bad_ids=set(a.loc[a.arm.isin(targets),'bin_id'])
bad_fixed=fixed[fixed.bin_id.isin(bad_ids)]
assert len(bad_fixed)==len(bad_ids)*len(samples)
assert not bad_fixed.use.any() and bad_fixed[['copynumber','segmented']].isna().all().all()
rows=[];read_rows=[]
normal=a[a.use&a.chromosome.isin(['chr'+str(i) for i in range(1,23)])][samples].median()
for arm in arms:
    z=a[a.arm.eq(arm)]
    src=source[source.arm.eq(arm)]
    assert len(src)==7 and src.bins.eq(len(z)).all(),arm
    rows.append(dict(arm=arm,total_100kb_bins=len(z),bins_used_in_old_heatmap=int(src.bins.iloc[0]),
        bins_passing_full_quality_mask=int(z.use.sum()),
        bins_passing_mappability_55=int(z.mappability.ge(55).sum()),
        bins_passing_difficulty_overlap_below_20=int(z.blacklist.lt(20).sum()),
        median_mappability_percent=float(z.mappability.median()),
        minimum_difficulty_overlap_percent=float(z.blacklist.min()),
        median_difficulty_overlap_percent=float(z.blacklist.median()),
        control_total_cached_reads=int(z['LNCaP-sgCTRL'].sum()),
        control_bins_with_reads=int(z['LNCaP-sgCTRL'].gt(0).sum()),
        control_median_reads_per_bin=float(z['LNCaP-sgCTRL'].median()),
        control_median_counts_relative_to_usable_autosome=float(z['LNCaP-sgCTRL'].median()/normal['LNCaP-sgCTRL'])))
    for sample in samples:
        read_rows.append(dict(arm=arm,sample=sample,total_cached_reads=int(z[sample].sum()),
            bins_with_reads=int(z[sample].gt(0).sum()),median_reads_per_bin=float(z[sample].median()),
            minimum_reads_per_bin=int(z[sample].min()),maximum_reads_per_bin=int(z[sample].max()),
            median_counts_relative_to_usable_autosome=float(z[sample].median()/normal[sample])))
summary=pd.DataFrame(rows).set_index('arm').loc[arms].reset_index()
summary.to_csv(OUT/'arm_bin_support.csv',index=False)
pd.DataFrame(read_rows).to_csv(OUT/'arm_read_support_all_samples.csv',index=False)
a[a.arm.isin(targets+comparators)].to_csv(OUT/'target_and_comparator_bin_details.csv',index=False)
manifest={k:{'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for k,p in paths.items()}
(OUT/'input_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
lines=['# Supp Fig 3.11: acrocentric short-arm support audit','',
'The five masked arms have **zero bins passing the existing full quality mask**. The original heatmap had nevertheless included every nominal bin on these arms. This supports retaining the display mask, but the reason is unreliable arm-level measurement, not absence of sequencing reads.', '',
'| Arm | Bins included originally | Mappability >=55% | Full quality mask | Median annotated mappability |',
'|---|---:|---:|---:|---:|']
for arm in targets+comparators:
    r=summary.set_index('arm').loc[arm]
    lines.append(f"| {arm[3:]} | {int(r.total_100kb_bins)} | {int(r.bins_passing_mappability_55)} | {int(r.bins_passing_full_quality_mask)} | {r.median_mappability_percent:.2f}% |")
lines += ['',
'The saved mask requires non-N sequence >=60%, annotated mappability >=55%, and difficulty-region overlap <20%. Every bin on the five acrocentric short arms fails the difficulty-overlap criterion; only 0-7 bins per arm meet the mappability criterion alone. These thresholds predate the heatmap masking and were not chosen using the observed arm differences.', '',
'## Read counts are present','',
'| Arm | Control bins with counted reads | Control total cached read counts | Median counts per 100-kb bin | Median / usable-autosomal median |',
'|---|---:|---:|---:|---:|']
for arm in targets+comparators:
    r=summary.set_index('arm').loc[arm]
    lines.append(f"| {arm[3:]} | {int(r.control_bins_with_reads)}/{int(r.total_100kb_bins)} | {int(r.control_total_cached_reads):,} | {r.control_median_reads_per_bin:g} | {r.control_median_counts_relative_to_usable_autosome:.3f} |")
lines += ['',
'These are cached QDNAseq alignment counts before GC/mappability correction and bin-quality exclusion, using the original read-counting settings. They are not sequencing depth, absolute copy number, or independent observations of an arm-level change. Ratios here only compare median bin counts with the same sample\'s usable autosomal-bin median. All eight libraries are tabulated in `arm_read_support_all_samples.csv`.', '',
'## Why the old heatmap still contained strong colours','',
'The historical `wgs_bins.csv` has `use=TRUE` for all bins, including these arms. Its arm summarizer therefore included 165, 114, 172, 111 and 142 bins for 13p, 14p, 15p, 21p and 22p respectively. This is consistent with the previously diagnosed QDNAseq filter-reset problem. The corrected staged DU145 output marks all these bins excluded and their CN/segment values missing in every library.', '',
'The existing heatmap uses medians of saved segment differences. Having many nominal bins (often sharing segment values) does not establish many independent reliable measurements. This audit does not prove that every signal on these arms is artefactual, or that biological changes are absent; it shows that these data do not support interpreting them as well-resolved arm-level CN changes.', '',
'The grey RNA cells use a separate availability rule (at least 20 genes per arm) and should not be used as evidence of WGS coverage or mappability. No inference of zero genes is made from missing RNA summaries.', '',
'## Scope and reproduction','',
'- No figure, caption, manuscript, source arm table or historical analysis was changed in this audit.',
'- Reproduced original arm bin counts for every arm and all seven heatmap DNA rows; checked the corrected exclusion in all eight DU145 libraries.',
'- Arm assignment follows the existing CHM13v2.0 cytoband/midpoint rule used by the figure source.',
'- Run `export_counts.R` with the existing R 4.5.2/QDNAseq environment, then `audit_arms.py` with the existing pandas environment.',
'- See `arm_bin_support.csv`, `arm_read_support_all_samples.csv`, `target_and_comparator_bin_details.csv`, and `input_manifest.json` for the numerical evidence and input checksums.',
'- The staged corrected WGS results remain separate from adopted draft-2 numerical figures; this check does not adopt that reanalysis.', '']
(OUT/'README.md').write_text('\n'.join(lines),encoding='utf-8')
print(summary[summary.arm.isin(targets+comparators)].to_string(index=False))
print('All original arm counts reproduced; corrected exclusion verified for all eight libraries.')
