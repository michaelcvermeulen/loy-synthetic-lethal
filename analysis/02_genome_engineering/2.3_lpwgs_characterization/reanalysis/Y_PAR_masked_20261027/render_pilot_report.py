"""Render the completed pilot using the Windows plotting runtime."""
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'comparison'
assert (ROOT/'qc/REALIGNMENT_COMPLETE.txt').exists()
source_checksum=json.loads((ROOT/'qc/source_bam_sha256.json').read_text(encoding='utf-8-sig'))
source_bam=ROOT.parents[1]/'t2t_bam/A549-CTRL2.bam'
digest=hashlib.sha256()
with source_bam.open('rb') as source:
    for block in iter(lambda:source.read(8*1024*1024),b''):digest.update(block)
assert digest.hexdigest().upper()==source_checksum['Hash'],'Original BAM checksum changed'
(ROOT/'qc/source_bam_unchanged.json').write_text(json.dumps({'path':str(source_bam),
    'sha256':digest.hexdigest(),'matches_before_pilot':True},indent=2))
meta=json.loads((OUT/'analysis_parameters.json').read_text())
cov=pd.read_csv(OUT/'regional_coverage.csv')
hist=pd.read_csv(OUT/'PAR_MAPQ_distribution.csv')
norm=pd.read_csv(OUT/'autosomal_normalization.csv')
combo=pd.read_csv(OUT/'combined_PAR_coverage.csv')
assert set(cov.alignment)=={'original','Y_PAR_masked'}
old=json.loads((ROOT/'qc/A549-CTRL2.original_read_qc.json').read_text())
new=json.loads((ROOT/'qc/A549-CTRL2.masked_read_qc.json').read_text())
fq=json.loads((ROOT/'qc/A549-CTRL2.fastq_validation.json').read_text())
assert old['read_fingerprint']==new['read_fingerprint']==fq['read_fingerprint']
move=json.loads((OUT/'read_reassignment_summary.json').read_text())

rows=[]
for par in ('PAR1','PAR2'):
    for label in ('original','Y_PAR_masked'):
        h=hist[hist.region.str.endswith('_'+par)&hist.alignment.eq(label)]
        total=int(h.read_count.sum());hq=int(h.loc[h.MAPQ.ge(37),'read_count'].sum())
        rel=combo[combo.region.eq(par)&combo.alignment.eq(label)&combo.minimum_MAPQ.eq(37)].iloc[0]
        rows.append({'region':par,'alignment':label,'all_MAPQ_read_ends':total,'MAPQ37_read_ends':hq,
                     'MAPQ37_percent':100*hq/total,'autosome_normalized_PAR_depth_Q37':rel.relative_to_median_autosomal_depth})
summary=pd.DataFrame(rows);summary.to_csv(OUT/'pilot_summary.csv',index=False)
fonts=ROOT.parents[2]/'SEPT_11_2026_THESIS/skills/thesis-figure-style/assets/fonts'
for p in fonts.glob('*.ttf'):font_manager.fontManager.addfont(str(p))
plt.rcParams.update({'font.family':'Open Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,
                     'pdf.fonttype':42,'svg.fonttype':'none'})
fig,axes=plt.subplots(1,2,figsize=(7.5,3.4),layout='constrained')
for label,color,offset,legend in [('original','#999999',-.18,'Original reference'),('Y_PAR_masked','#7E3F8F',.18,'Y-PAR masked')]:
    d=summary[summary.alignment.eq(label)].set_index('region').loc[['PAR1','PAR2']]
    x=np.arange(2)+offset
    bars=axes[0].bar(x,d.MAPQ37_percent,width=.34,color=color,label=legend)
    axes[0].bar_label(bars,labels=[f'{v:.1f}%' for v in d.MAPQ37_percent],padding=3,fontsize=8)
    bars=axes[1].bar(x,d.autosome_normalized_PAR_depth_Q37,width=.34,color=color)
    axes[1].bar_label(bars,labels=[f'{v:.2f}' for v in d.autosome_normalized_PAR_depth_Q37],padding=3,fontsize=8)
axes[0].set(xticks=[0,1],xticklabels=['PAR1','PAR2'],ylim=(0,115),ylabel='PAR read ends with MAPQ ≥37 (%)')
fig.legend(*axes[0].get_legend_handles_labels(),frameon=False,loc='outside upper center',ncol=2,fontsize=8)
axes[1].set(xticks=[0,1],xticklabels=['PAR1','PAR2'],ylabel='Combined PAR / autosomal depth\n(MAPQ ≥37)')
axes[1].set_ylim(0,summary.autosome_normalized_PAR_depth_Q37.max()*1.23)
for i,ax in enumerate(axes):ax.text(-.13,1.04,'ab'[i],transform=ax.transAxes,fontsize=12,fontweight='semibold')
for ext in ('png','pdf','svg'):
    fig.savefig(OUT/f'A549_CTRL2_PAR_masking_pilot.{ext}',dpi=600)
plt.close(fig)

lines=['# A549-CTRL2 Y-PAR-masked alignment pilot','',
'One Y-retaining A549 control sample was recovered from its existing BAM and realigned with BWA-MEM2 2.2.1. The comparison below tests PAR read assignment; it does **not** estimate Y-specific retention or compare a depleted population with its control. No other sample or draft-2 figure was updated.', '',
f"Recovered and validated **{fq['pairs']:,} read pairs** ({old['counts']['primary']:,} primary read ends). Read-name, mate, forward-sequence and quality fingerprints agree between the original BAM, paired FASTQs and new BAM. No singleton reads were lost. The input is the previously processed/trimmed sequence stored in the BAM, not newly recovered raw sequencer output.", '',
'## PAR mapping and coverage','',
'| Region | Original high-MAPQ / all PAR read ends | Masked high-MAPQ / all PAR read ends | Original → masked PAR/autosomal depth |',
'|---|---:|---:|---:|']
for par in ('PAR1','PAR2'):
    d=summary[summary.region.eq(par)].set_index('alignment');a=d.loc['original'];b=d.loc['Y_PAR_masked']
    lines.append(f"| {par} | {int(a.MAPQ37_read_ends):,} / {int(a.all_MAPQ_read_ends):,} ({a.MAPQ37_percent:.1f}%) | {int(b.MAPQ37_read_ends):,} / {int(b.all_MAPQ_read_ends):,} ({b.MAPQ37_percent:.1f}%) | {a.autosome_normalized_PAR_depth_Q37:.3f} → {b.autosome_normalized_PAR_depth_Q37:.3f} |")
lines += ['', 'High MAPQ means ≥37, matching the original QDNAseq read-counting threshold. Counts combine X- and Y-PAR primary alignments; after masking, homologous reads are expected to align to X. The full tables also provide MAPQ 0 and 30 sensitivity results.', '',
'Masking substantially restores high-MAPQ PAR signal in this control, supporting ambiguous X/Y-PAR mapping as the main reason for its low filtered PAR coverage on the original reference. This pilot establishes a workable BAM-to-FASTQ realignment route. A matched Y-depleted sample is still needed to assess changes in combined PAR dosage; the masked alignment cannot identify the X or Y origin of homologous PAR reads.', '',
f"![PAR masking pilot]({(OUT/'A549_CTRL2_PAR_masking_pilot.png').as_posix()})", '',
'The combined coverage diagnostic sums aligned bases across X and Y PAR intervals, divides by the X-PAR span, then normalizes to the median depth of 21,446 usable autosomal 100-kb bins. X and Y PAR spans differ in length; this is a descriptive combined-signal measure, not an absolute copy-number call. Overlapping paired reads can both contribute. No GC correction or masked-reference mappability model has been fitted for this pilot, and the old unmasked-reference PAR exclusion mask is not reused to assess newly mapped PAR reads.', '',
'## Tracking the original Y-PAR reads','',
'| Original location | Read ends followed | Remapped to corresponding X-PAR | Previously MAPQ <37, now X-PAR MAPQ ≥37 | Unmapped after realignment |',
'|---|---:|---:|---:|---:|']
for key in ('chrY_PAR1','chrY_PAR2'):
    x=move[key]
    lines.append(f"| {key} | {x['original_PAR_read_ends']:,} | {x['new_X_PAR_read_ends']:,} | {x['old_below37_new_X_PAR_ge37']:,} | {x['new_unmapped']:,} |")
lines += ['', 'All original PAR primary read ends were followed through the complete new alignment. Detailed destinations and before/after MAPQs are in `comparison/PAR_read_reassignment.csv`.', '',
'## Whole-library checks','',
'| Metric | Original | Y-PAR masked |','|---|---:|---:|']
for field,label in [('primary','Primary read ends'),('primary_unmapped','Unmapped primary read ends'),('primary_properly_paired','Properly paired primary read ends')]:
    lines.append(f"| {label} | {old['counts'].get(field,0):,} | {new['counts'].get(field,0):,} |")
for q in (30,37):
    d=norm[norm.minimum_MAPQ.eq(q)].set_index('alignment')
    lines.append(f"| Median autosomal depth, MAPQ ≥{q} | {d.loc['original','median_autosomal_depth']:.5f} | {d.loc['Y_PAR_masked','median_autosomal_depth']:.5f} |")
auto=pd.read_csv(OUT/'original.autosomal_bedcov.Q37.tsv',sep='\t',header=None)
after=pd.read_csv(OUT/'Y_PAR_masked.autosomal_bedcov.Q37.tsv',sep='\t',header=None)
assert auto.iloc[:,:4].equals(after.iloc[:,:4])
a=auto.iloc[:,4]/(auto.iloc[:,2]-auto.iloc[:,1]);b=after.iloc[:,4]/(after.iloc[:,2]-after.iloc[:,1])
auto_stats={'common_autosomal_bins':len(a),'coverage_pearson':a.corr(b),'median_absolute_depth_change':(b-a).abs().median(),'p95_absolute_depth_change':(b-a).abs().quantile(.95),
    'bins_with_changed_coverage':int(a.ne(b).sum()),'percent_identical_bins':100*float(a.eq(b).mean()),'maximum_absolute_depth_change':float((b-a).abs().max())}
(OUT/'autosomal_stability.json').write_text(json.dumps(auto_stats,indent=2))
lines += ['',f"Across the common autosomal bins, old/new MAPQ ≥37 depth has Pearson r = {auto_stats['coverage_pearson']:.6f}; {auto_stats['percent_identical_bins']:.2f}% of bins have exactly unchanged aligned-base totals. The median absolute depth change is {auto_stats['median_absolute_depth_change']:.6f}, its 95th percentile is {auto_stats['p95_absolute_depth_change']:.6f}, and its maximum is {auto_stats['maximum_absolute_depth_change']:.6f}.", '',
'## Reference and alignment provenance','',
'- Original BAM: `../../t2t_bam/A549-CTRL2.bam`. It includes paired primary reads with sequences and qualities; no primary hard clipping, QC-failure flags or duplicate flags were found.',
'- The original BAM SHA256 is unchanged after the pilot; the check is recorded in `qc/source_bam_unchanged.json`.',
'- Local source reference: `H:/reference/chm13v2.0.fa.gz`. Only chrY PAR1 [0, 2,458,320) and PAR2 [62,122,809, 62,460,029) were replaced with Ns; 2,795,540 reference positions were masked.',
'- The prepared reference matches every contig of the official `chm13v2.0_maskedY.fa.gz` analysis set by sequence MD5. The downloaded official compressed-file MD5 is `bd90ddb80c86af7fcfeefe7a0909b175`. This uses the original mitochondrial sequence, not the separate rCRS replacement reference.',
'- All non-Y bases remain unchanged relative to the local source FASTA. Reference chromosome names and lengths match the original BAM; original BAM sequence dictionary tags do not independently provide a checksum of the historical HPC FASTA.',
'- BWA-MEM2 2.2.1 and samtools 1.23.1. Mapping uses eight workers and fixed `-K 160000000`, the original 16-worker default batch size. Name collation changes the input order; minor placement differences outside PAR can therefore occur. Autosomal stability is checked above.',
'- Reads are not trimmed again. All primary reads, including unmapped and MAPQ-0 reads, are recovered; secondary/supplementary records are omitted to avoid duplicating sequence records.',
'- The index is stored in the dedicated WSL scratch folder `/home/micha/codex-work/wgs-ypar-pilot-20261027/reference/`. An extra 16-GiB task-specific swap file was used for index construction and is removed by the index script afterward.', '',
'## Files','',
'- `fastq/A549-CTRL2.R1.fastq.gz` and `fastq/A549-CTRL2.R2.fastq.gz`: recovered pairs.',
'- `alignment/A549-CTRL2.YPARmasked.bam` and `.bam.bai`: new alignment and index.',
'- `qc/`: original/new read QC, fingerprints, paired-FASTQ verification, headers and flagstats.',
'- `reference/`: official PAR BED, prepared reference, checksum evidence and masking provenance.',
'- `comparison/`: quantitative tables, read reassignment, coverage diagnostics and review plots.',
'- `runtime/environment_explicit.txt`: exact software package build list.',
'- `logs/`: extraction, reference preparation, indexing and alignment logs.', '',
'Reproduction: `setup_runtime.sh`, `prepare_reference.py`, `verify_official_reference.py`, `build_index.sh`, `recover_pilot.sh`, `align_pilot.sh`, then `compare_pilot.py` and `render_pilot_report.py`. The pilot extraction used `recover_from_collated.sh` to reuse a completed name-collated BAM; the recovered-read fingerprint validates its equivalence. Scripts use separate output files and refuse to overwrite completed FASTQs/BAMs.', '',
'Primary documentation: [T2T analysis-set definitions](https://github.com/marbl/CHM13#analysis-set), [samtools FASTQ recovery](https://www.htslib.org/doc/samtools-fasta.html), [BWA-MEM2](https://github.com/bwa-mem2/bwa-mem2).', '']
(ROOT/'README.md').write_text('\n'.join(lines),encoding='utf-8')
print(summary.to_string(index=False))
print(json.dumps(auto_stats,indent=2))
