"""Construct a reproducible common peak universe and count whole fragments."""
from bisect import bisect_left, insort
import importlib.util
import json
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('prep', ROOT/'code/02_prepare.py')
prep = importlib.util.module_from_spec(spec); spec.loader.exec_module(prep)


def interval_counts(fragment_starts, fragment_ends, starts, ends):
    """Count each fragment once per overlapping feature (BED half-open)."""
    return (np.searchsorted(np.sort(fragment_starts), ends, side='left')
            - np.searchsorted(np.sort(fragment_ends), starts, side='right'))


def union_intervals(frame):
    rows = []
    for chrom, group in frame.groupby('chrom', sort=False):
        merged = []
        for start, end in group[['start','end']].sort_values('start').itertuples(index=False, name=None):
            if merged and start <= merged[-1][1]:
                merged[-1][1] = max(end, merged[-1][1])
            else:
                merged.append([start, end])
        rows.extend((chrom, start, end) for start, end in merged)
    return pd.DataFrame(rows, columns=['chrom','start','end'])


def fragments_in_union(fragment_starts, fragment_ends, starts, ends):
    """Count a fragment once across the entire union, including long fragments."""
    starts, ends = np.asarray(starts), np.asarray(ends)
    if len(starts) == 0:
        return 0
    order = np.argsort(starts); starts, ends = starts[order], ends[order]
    index = np.searchsorted(starts, fragment_ends, side='left') - 1
    valid = index >= 0
    return int(np.sum(valid & (ends[np.maximum(index, 0)] > fragment_starts)))


def main():
    params = json.loads((ROOT/'config/params.json').read_text())
    samples = pd.read_csv(ROOT/'config/samples.tsv', sep='\t')
    sizes = dict(pd.read_csv(ROOT/'config/bam.chrom.sizes', sep='\t', header=None).itertuples(index=False, name=None))
    blacklist = prep.load_intervals(ROOT/'config/hg38-blacklist.v2.bed')
    par = prep.load_intervals(ROOT/'config/hg38_PAR.bed')
    columns = ['chrom','start','end','name','score','strand','signal','pvalue','qvalue','summit_offset']
    peaksets = {}; candidates = []
    for sample in samples['sample']:
        peaks = pd.read_csv(ROOT/f'peaks/{sample}/{sample}_peaks.narrowPeak', sep='\t', names=columns)
        peaks = peaks[peaks.chrom.isin(params['canonical_chromosomes'])].copy()
        peaks = peaks[[not prep.overlaps(blacklist, r.chrom, r.start, r.end) for r in peaks.itertuples()]]
        if peaks.empty:
            raise ValueError(f'No peaks survive blacklist filtering for {sample}')
        peaks[columns].to_csv(ROOT/f'peaks/{sample}/{sample}.blacklist_filtered.narrowPeak',
                              sep='\t',header=False,index=False)
        peaks['summit'] = peaks.start + peaks.summit_offset
        peaks['rank_percentile'] = peaks.pvalue.rank(method='average', pct=True)
        peaksets[sample] = peaks
        for row in peaks.itertuples():
            start, end = int(row.summit)-250, int(row.summit)+251
            if row.summit_offset < 0 or start < 0 or end > sizes[row.chrom]:
                continue
            if not prep.overlaps(blacklist, row.chrom, start, end):
                candidates.append((row.chrom,start,end,sample,float(row.rank_percentile)))
    consensus = pd.DataFrame(candidates, columns=['chrom','start','end','summit_sample','rank_percentile'])
    # Reproducibility: a called summit lies inside this candidate's 501 bp window.
    for sample, peaks in peaksets.items():
        support = np.zeros(len(consensus), dtype=bool)
        for chrom, group in consensus.groupby('chrom', sort=False):
            summits = np.sort(peaks.loc[peaks.chrom == chrom, 'summit'].to_numpy())
            support[group.index] = (np.searchsorted(summits, group.end, side='left')
                                    > np.searchsorted(summits, group.start, side='left'))
        consensus[sample] = support.astype(int)
    ctrl = samples.loc[samples.condition == 'CTRL', 'sample'].tolist()
    loy = samples.loc[samples.condition == 'LOY', 'sample'].tolist()
    consensus['support_CTRL'] = consensus[ctrl].sum(axis=1)
    consensus['support_LOY'] = consensus[loy].sum(axis=1)
    consensus = consensus[(consensus.support_CTRL >= 2) | (consensus.support_LOY >= 2)]
    consensus = consensus.sort_values(['rank_percentile','chrom','start','summit_sample'],
                                      ascending=[False,True,True,True], kind='stable')
    occupied = {}; kept = []
    for row in consensus.itertuples():
        selected = occupied.setdefault(row.chrom, [])
        index = bisect_left(selected, row.start)
        if index > 0 and selected[index-1] + 501 > row.start:
            continue
        if index < len(selected) and selected[index] < row.end:
            continue
        insort(selected, row.start); kept.append(row.Index)
    consensus = consensus.loc[kept].copy()
    order = {c:i for i,c in enumerate(params['canonical_chromosomes'])}
    consensus['chrom_order'] = consensus.chrom.map(order)
    consensus = consensus.sort_values(['chrom_order','start']).drop(columns='chrom_order').reset_index(drop=True)
    consensus.insert(0, 'peak_id', consensus.chrom+':'+consensus.start.astype(str)+'-'+consensus.end.astype(str))
    consensus['is_PAR'] = [prep.overlaps(par,r.chrom,r.start,r.end) for r in consensus.itertuples()]
    consensus['primary_test_region'] = (consensus.chrom != 'chrY') & ~consensus.is_PAR
    assert consensus.peak_id.is_unique and (consensus.end-consensus.start).eq(501).all()
    consensus.to_csv(ROOT/'peaks/consensus_membership.tsv', sep='\t', index=False)
    consensus[['chrom','start','end','peak_id']].to_csv(ROOT/'peaks/consensus.bed', sep='\t', header=False,index=False)
    all_peaks = union_intervals(pd.concat(list(peaksets.values())))
    background = []
    for chrom in params['canonical_chromosomes']:
        if chrom in ['chrX','chrY']:
            continue
        chrom_peaks = all_peaks[all_peaks.chrom == chrom]
        peak_intervals = {chrom: (chrom_peaks.start.tolist(), chrom_peaks.end.tolist())}
        for start in range(0, sizes[chrom]-9999, 10000):
            end = start+10000
            if not prep.overlaps(blacklist,chrom,start,end) and not prep.overlaps(peak_intervals,chrom,start,end):
                background.append((f'{chrom}:{start}-{end}',chrom,start,end))
    background = pd.DataFrame(background,columns=['bin_id','chrom','start','end'])
    consensus_union = union_intervals(consensus)
    count_matrix = pd.DataFrame({'peak_id':consensus.peak_id})
    bin_matrix = background.copy()
    qc_rows = []
    for sample in samples['sample']:
        print('Counting', sample, flush=True)
        fragments = np.load(ROOT/f'counts/{sample}_fragments.npz')
        counts = np.zeros(len(consensus),dtype=np.int64)
        bin_counts = np.zeros(len(background),dtype=np.int64)
        own_union = union_intervals(peaksets[sample])
        own_rip = shared_rip = 0
        for chrom in params['canonical_chromosomes']:
            start,end = fragments[chrom+'_start'],fragments[chrom+'_end']
            subset = consensus[consensus.chrom == chrom]
            counts[subset.index] = interval_counts(start,end,subset.start,subset.end)
            bins = background[background.chrom == chrom]
            bin_counts[bins.index] = interval_counts(start,end,bins.start,bins.end)
            own = own_union[own_union.chrom == chrom]
            shared = consensus_union[consensus_union.chrom == chrom]
            own_rip += fragments_in_union(start,end,own.start,own.end)
            shared_rip += fragments_in_union(start,end,shared.start,shared.end)
        count_matrix[sample] = counts; bin_matrix[sample] = bin_counts
        qc = json.loads((ROOT/f'qc/{sample}_preparation.json').read_text())
        qc_rows.append({'sample':sample,'n_peaks':len(peaksets[sample]),
                       'accepted_fragments':qc['accepted_fragments'],
                       'fragments_in_own_peaks':own_rip,'frip_own':own_rip/qc['accepted_fragments'],
                       'fragments_in_consensus':shared_rip,'frip_consensus':shared_rip/qc['accepted_fragments'],
                       'extra_consensus_assignments':int(counts.sum()-shared_rip)})
    count_matrix.to_csv(ROOT/'counts/peak_counts.tsv',sep='\t',index=False)
    bin_matrix.to_csv(ROOT/'counts/background_bin_counts.tsv',sep='\t',index=False)
    pd.DataFrame(qc_rows).to_csv(ROOT/'qc/peak_qc.tsv',sep='\t',index=False)
    (ROOT/'counts/COUNTING_METHOD.txt').write_text(
      'Integer whole-fragment overlaps on the same BED intervals for every sample.\n'
      'A fragment is counted once per overlapping peak; multiple peak assignments are allowed.\n'
      'FRiP counts each fragment once across the union, regardless of how many peaks it overlaps.\n'
      'Extra consensus assignments report the difference between the count sum and unique union overlaps.\n'
      'Background uses 10 kb autosomal bins not intersecting any sample peak or blacklist interval.\n'
      'Primary differential testing excludes chrY and X/Y PAR intervals.\n')
    (ROOT/'counts/COMPLETE.json').write_text(json.dumps({
        'samples':samples['sample'].tolist(), 'consensus_peaks':len(consensus),
        'background_bins':len(background)}, indent=2)+'\n')
    print(len(consensus),'consensus peaks;',len(background),'background bins',flush=True)


if __name__ == '__main__':
    main()
