"""Tabulate regional CPM and explain its relationship to atlas bin medians."""
from pathlib import Path
import csv,json,statistics

ROOT=Path(__file__).resolve().parent
THESIS=ROOT.parents[2]/'SEPT_11_2026_THESIS'
OUT=ROOT/'regional_cpm_20261028'
PAIRS={'A549':['A549-CTRL2','A549-LOY'],'HT1080':['HT1080-sgCTRL-D0','HT1080-LOY-D0'],
       'DU145':['LNCaP-sgCTRL','LNCaP-LOY-2X']}
REGIONS=['MSY_Yp','PAR1_X','euchromatic_Yq','PAR2_X']
LABELS={'MSY_Yp':'Male-specific Yp','PAR1_X':'PAR1 (X-mapped)','euchromatic_Yq':'Euchromatic Yq',
        'PAR2_X':'PAR2 (X-mapped)','distal_MSY_Yp':'Distal MSY Yp','proximal_MSY_Yp':'Proximal MSY Yp'}
def read(path):
    with path.open(newline='',encoding='utf-8-sig') as f:return list(csv.DictReader(f))
def save(path,data):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)

counts=read(OUT/'regional_reads_and_cpm_all_definitions.csv')
ratios=read(OUT/'CPM_ratios_vs_figure_medians.csv')
denom={r['sample']:int(r['eligible_library_read_ends']) for r in counts}
bin_rows=[];comparison=[]
for cell,(control,loy) in PAIRS.items():
    base=THESIS/'figures/chapter3/drafts'
    source=(base/'a549_dense_bins_20261027/100kb' if cell=='A549' else base/'ht1080_du145_ypar_20261027'/cell.lower())/'source_data/regional_y_atlases_normalised'/cell.lower()
    coverage=read(source/'binned_coverage.csv')
    lookup={(r['sample'],r['chromosome'],int(r['start0']),int(r['end'])):int(r['read_count']) for r in coverage}
    y=[r for r in read(source/'paired_Y_bins_and_exclusions.csv') if r['included'].lower()=='true']
    par=[r for r in read(source/'X_PAR_bin_ratios.csv') if r['included'].lower()=='true']
    groups={r:[] for r in REGIONS}
    if cell=='DU145':groups.update({'distal_MSY_Yp':[],'proximal_MSY_Yp':[]})
    name_map={'MSY Yp':'MSY_Yp','euchromatic Yq':'euchromatic_Yq',
              'distal MSY Yp':'distal_MSY_Yp','proximal MSY Yp':'proximal_MSY_Yp'}
    for r in y:
        interval=('chrY',int(float(r['start']))-1,int(float(r['end'])))
        region=name_map[r['region']]
        groups[region].append(interval)
        if region in ['distal_MSY_Yp','proximal_MSY_Yp']:groups['MSY_Yp'].append(interval)
    for r in par:groups[r['track']].append(('chrX',int(r['start0']),int(r['end'])))
    for region,intervals in groups.items():
        entries=[]
        for chrom,start,end in intervals:
            a=lookup[(control,chrom,start,end)];b=lookup[(loy,chrom,start,end)]
            assert a>0
            entries.append({'cell':cell,'region':region,'chromosome':chrom,'start0':start,'end':end,
                'control_read_ends':a,'LOY_read_ends':b,
                'raw_bin_CPM_percent_of_control':100*(b/denom[loy])/(a/denom[control])})
        bin_rows+=entries
        whole=next(r for r in ratios if r['cell']==cell and r['scope']=='whole_region' and r['region']==region)
        plotted=next(r for r in ratios if r['cell']==cell and r['scope']=='plotted_bins' and r['region']==region)
        higher=[r for r in entries if r['raw_bin_CPM_percent_of_control']>20]
        comparison.append({'cell':cell,'region':region,'n_bins':len(entries),
            'whole_region_CPM_percent_of_control':float(whole['LOY_CPM_percent_of_control']),
            'plotted_union_CPM_percent_of_control':float(plotted['LOY_CPM_percent_of_control']),
            'median_raw_bin_CPM_percent_of_control':statistics.median(r['raw_bin_CPM_percent_of_control'] for r in entries),
            'figure_median_percent_of_control':plotted['figure_median_percent_of_control'],
            'bins_with_zero_LOY_reads':sum(r['LOY_read_ends']==0 for r in entries),
            'bins_above_20_percent_raw_retention':len(higher),
            'fraction_of_summed_LOY_bin_counts_in_bins_above20':sum(r['LOY_read_ends'] for r in higher)/sum(r['LOY_read_ends'] for r in entries)})
save(OUT/'per_bin_raw_CPM_ratios.csv',bin_rows)
save(OUT/'figure_comparison_summary.csv',comparison)

table=['| Sample | Male-specific Yp | PAR1 (X-mapped) | Euchromatic Yq | PAR2 (X-mapped) |',
       '|---|---:|---:|---:|---:|']
cpmtable=[table[0],table[1]];readtable=[table[0],table[1]]
wide=[]
for cell,pair in PAIRS.items():
    for condition,sample in zip(['sgCTRL','Y-depleted'],pair):
        z={r['region']:r for r in counts if r['sample']==sample and r['scope']=='whole_region'}
        label=f'{cell} {condition}'
        table.append('| '+label+' | '+' | '.join(f"{int(z[r]['primary_read_ends']):,} ({float(z[r]['CPM']):,.2f})" for r in REGIONS)+' |')
        cpmtable.append('| '+label+' | '+' | '.join(f"{float(z[r]['CPM']):,.2f}" for r in REGIONS)+' |')
        readtable.append('| '+label+' | '+' | '.join(f"{int(z[r]['primary_read_ends']):,}" for r in REGIONS)+' |')
        w={'cell':cell,'condition':condition,'sample':sample,'eligible_library_read_ends':denom[sample]}
        for region in REGIONS:
            w[region+'_reads']=int(z[region]['primary_read_ends']);w[region+'_CPM']=float(z[region]['CPM'])
        wide.append(w)
save(OUT/'regional_reads_CPM_wide.csv',wide)
rat=['| Cell line | Male-specific Yp | PAR1 | Euchromatic Yq | PAR2 |','|---|---:|---:|---:|---:|']
for cell in PAIRS:
    v={r['region']:r for r in ratios if r['cell']==cell and r['scope']=='whole_region'}
    rat.append('| '+cell+' | '+' | '.join(f"{float(v[r]['LOY_CPM_percent_of_control']):.2f}%" for r in REGIONS)+' |')
detail=['| Cell / region | Whole-region CPM ratio | Plotted-bin CPM ratio | Median raw bin CPM ratio | Figure median |','|---|---:|---:|---:|---:|']
for r in comparison:
    figure=f"{float(r['figure_median_percent_of_control']):.2f}%" if r['figure_median_percent_of_control'] else 'split in figure'
    detail.append(f"| {r['cell']} / {LABELS[r['region']]} | {r['whole_region_CPM_percent_of_control']:.2f}% | {r['plotted_union_CPM_percent_of_control']:.2f}% | {r['median_raw_bin_CPM_percent_of_control']:.2f}% | {figure} |")
du=next(r for r in comparison if r['cell']=='DU145' and r['region']=='proximal_MSY_Yp')
text=f'''# Regional read totals and CPM — six realigned atlas samples

Values below are **primary read-end counts (CPM)** from the completed Y-PAR-masked BAMs. Both numerator and denominator require MAPQ >=37 and exclude unmapped, secondary, supplementary, QC-failed and marked-duplicate records (`-F 0xF04`). R1 and R2 are counted separately. The denominator is all eligible mapped primary read ends across all reference contigs, including chrM; it is not the number of pairs or QDNAseq's used-read total.

CPM = regional eligible read ends / total eligible mapped read ends x 1,000,000.

## Whole-region counts and CPM

{chr(10).join(table)}

These are totals over the full stated intervals, without the figure's bin mappability or control-count exclusions. Male-specific Yp is chrY:2,458,321–10,565,750. Euchromatic Yq is chrY:10,883,086–27,449,937, matching the figure's Yq region. PAR1 and PAR2 are counted on chrX:1–2,394,410 and chrX:153,925,835–154,259,566, respectively. Coordinates here are one-based inclusive. PAR values are combined X/Y signal after Y-PAR masking, not Y-specific PAR counts.

The complete CSV also contains Yq including Yq12 (excluding PAR2), Yq12 alone, and DU145's distal/proximal Yp subdivisions. Reads overlapping a subdivision boundary may count in each subdivision; the combined region counts each eligible alignment once.

## CPM in Y-depleted cells as a percentage of sgCTRL

{chr(10).join(rat)}

## Comparison with the figures

{chr(10).join(detail)}

Whole-region totals, totals restricted to plotted bins and medians across bins answer different questions. The Y figure lines are medians of QDNAseq-corrected bin ratios; PAR lines are medians of autosome-normalised bin-depth ratios. CPM uses a library-wide eligible-read denominator, without GC/mappability correction. Therefore the percentages need not equal the plotted labels.

DU145 proximal Yp illustrates the distinction: its whole-region CPM ratio is {du['whole_region_CPM_percent_of_control']:.2f}%, its plotted-bin union ratio is {du['plotted_union_CPM_percent_of_control']:.2f}%, and its median raw bin ratio is {du['median_raw_bin_CPM_percent_of_control']:.2f}%, compared with the figure median of {float(du['figure_median_percent_of_control']):.2f}%. Of the {du['n_bins']} plotted proximal bins, {du['bins_above_20_percent_raw_retention']} have raw retention above 20%; these account for {100*du['fraction_of_summed_LOY_bin_counts_in_bins_above20']:.1f}% of the summed per-bin Y-depleted counts. Their higher counts raise the pooled total relative to the median. Per-bin counts are supplied as a diagnostic and can double-count reads crossing bin boundaries; regional and plotted-union totals above use the multi-region iterator to count each alignment once per union.

HT1080 has nonzero Yq read totals even though its corrected median bin value is zero (displayed as <1%); the plot does not claim complete absence of Yq reads.

## Denominators

'''
text+='| Sample ID | Eligible mapped primary read ends |\n|---|---:|\n'
for sample,n in denom.items():text+=f'| {sample} | {n:,} |\n'
text+='''
## Files and validation

- `regional_reads_CPM_wide.csv`: compact six-sample table with both counts and CPM.
- `whole_region_reads_and_cpm.csv`: the four requested regions in long format.
- `regional_reads_and_cpm_all_definitions.csv`: whole regions, additional Yq/Yp definitions and exact plotted-bin unions.
- `CPM_ratios_vs_figure_medians.csv`, `figure_comparison_summary.csv` and `per_bin_raw_CPM_ratios.csv`: comparison diagnostics.
- `provenance.json` and `beds/`: filters, reference intervals, figure-input checksums and counted interval unions.

All six BAMs passed samtools quickcheck. X-PAR counts exactly reproduce the earlier coverage calculations; high-quality Y-mapped PAR counts are zero as expected from the reference mask. Plotted-union counts do not exceed their containing whole regions. BAM sizes and modification times are unchanged. These calculations do not modify the figures or draft 2.

Reproduction: run `count_regional_cpm.py` under the existing WSL runtime, then `render_cpm_summary.py` using Python. Legacy LNCaP sample IDs correspond to DU145.
'''
(OUT/'README.md').write_text(text,encoding='utf-8')
(OUT/'CPM_table.md').write_text('\n'.join(cpmtable)+'\n',encoding='utf-8')
(OUT/'read_counts_table.md').write_text('\n'.join(readtable)+'\n',encoding='utf-8')
print('\n'.join(table));print('\n'.join(rat));print(json.dumps(du,indent=2))
