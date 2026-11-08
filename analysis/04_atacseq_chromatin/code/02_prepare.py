"""Validate complete paired fragments, prepare BAMs, and compute fragment/TSS QC."""
from __future__ import annotations
import argparse
from array import array
from bisect import bisect_left
from collections import Counter, defaultdict
import csv
import json
from pathlib import Path
import tempfile
import time

import numpy as np
import pysam
import pyBigWig

ROOT = Path(__file__).resolve().parents[1]


def load_intervals(path):
    groups = defaultdict(list)
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith('#'):
                continue
            chrom, start, end, *_ = line.split()
            groups[chrom].append((int(start), int(end)))
    result = {}
    for chrom, intervals in groups.items():
        merged = []
        for start, end in sorted(intervals):
            if merged and start <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(end, merged[-1][1]))
            else:
                merged.append((start, end))
        result[chrom] = ([x[0] for x in merged], [x[1] for x in merged])
    return result


def overlaps(intervals, chrom, start, end):
    starts, ends = intervals.get(chrom, ([], []))
    index = bisect_left(starts, end) - 1
    return index >= 0 and ends[index] > start


def eligible_pairs(bam, canonical, blacklist, metrics, min_mapq=30):
    """Coordinate stream; retain both mates only when the complete pair is eligible.

    Pending records are flushed at chromosome boundaries. Mate keys include RG;
    multiple primary records for a mate are a fatal provenance ambiguity.
    """
    pending = {}
    previous = (-1, -1)
    previous_chrom = None
    for read in bam.fetch(until_eof=True):
        metrics['source_records'] += 1
        if read.is_unmapped:
            metrics['unmapped_records'] += 1
            continue
        coordinate = (read.reference_id, read.reference_start)
        if coordinate < previous:
            raise ValueError('Source BAM is not coordinate sorted')
        previous = coordinate
        chrom = read.reference_name
        if chrom != previous_chrom:
            metrics['unmatched_primary_records'] += len(pending)
            pending.clear()
            previous_chrom = chrom
        if chrom not in canonical:
            metrics['noncanonical_records'] += 1
            continue
        if read.is_secondary or read.is_supplementary:
            metrics['nonprimary_records'] += 1
            continue
        if (not read.is_paired or not read.is_proper_pair or read.mate_is_unmapped
                or read.next_reference_id != read.reference_id):
            metrics['improper_records'] += 1
            continue
        key = (read.query_name, read.get_tag('RG') if read.has_tag('RG') else '')
        mate = pending.pop(key, None)
        if mate is None:
            pending[key] = read
            continue
        metrics['complete_primary_pairs'] += 1
        if read.is_read1 == mate.is_read1 or read.is_read2 == mate.is_read2:
            raise ValueError(f'Multiple primary alignments for one mate: {key}')
        if (read.next_reference_start != mate.reference_start
                or mate.next_reference_start != read.reference_start):
            raise ValueError(f'Nonreciprocal mate coordinates: {key}')
        if read.is_duplicate or mate.is_duplicate:
            metrics['duplicate_pairs'] += 1
            continue
        if read.is_qcfail or mate.is_qcfail:
            metrics['qc_failed_pairs'] += 1
            continue
        if min(read.mapping_quality, mate.mapping_quality) < min_mapq:
            metrics['low_mapq_pairs'] += 1
            continue
        start = min(read.reference_start, mate.reference_start)
        end = max(read.reference_end, mate.reference_end)
        if end <= start or read.is_reverse == mate.is_reverse:
            metrics['invalid_geometry_pairs'] += 1
            continue
        if abs(read.template_length) != end - start or read.template_length != -mate.template_length:
            metrics['inconsistent_template_length_pairs'] += 1
            continue
        if overlaps(blacklist, chrom, start, end):
            metrics['blacklist_pairs'] += 1
            continue
        metrics['accepted_fragments'] += 1
        yield read, mate, chrom, start, end
    metrics['unmatched_primary_records'] += len(pending)


def tss_profile(cuts, sizes, blacklist, path):
    profile = np.zeros(4000, dtype=np.int64)
    used = set()
    with open(path) as handle:
        for chrom, start, end, gene, score, strand in csv.reader(handle, delimiter='\t'):
            center = int(start)
            key = (chrom, center, strand)
            if chrom == 'chrY' or chrom not in cuts or key in used:
                continue
            if center < 2000 or center + 2000 > sizes[chrom]:
                continue
            if overlaps(blacklist, chrom, center - 2000, center + 2000):
                continue
            positions = cuts[chrom]
            left, right = np.searchsorted(positions, [center - 2000, center + 2000])
            offsets = positions[left:right].astype(np.int64) - center + 2000
            counts = np.bincount(offsets, minlength=4000)
            profile += counts if strand == '+' else counts[::-1]
            used.add(key)
    flank = np.r_[profile[:100], profile[-100:]].mean()
    score = float(profile[1950:2050].mean() / flank) if flank > 0 else None
    return profile, len(used), score


def write_coverage(path, fragments, sizes, denominator, bin_size=25):
    with pyBigWig.open(str(path), 'w') as bw:
        bw.addHeader(list(sizes.items()))
        for chrom, length in sizes.items():
            start, end = fragments[chrom]
            n = (length + bin_size - 1) // bin_size
            sb, eb = start // bin_size, end // bin_size
            delta = (np.bincount(sb, minlength=n+1) - np.bincount(eb, minlength=n+1))
            integral = np.cumsum(delta)[:n].astype(float) * bin_size
            integral -= np.bincount(sb, weights=start % bin_size, minlength=n+1)[:n]
            integral += np.bincount(eb, weights=end % bin_size, minlength=n+1)[:n]
            widths = np.full(n, bin_size); widths[-1] = length - (n-1)*bin_size
            values = integral / widths * (1e6 / denominator)
            nz = np.flatnonzero(values > 0)
            for offset in range(0, len(nz), 250000):
                chunk = nz[offset:offset+250000]
                starts = chunk * bin_size
                bw.addEntries([chrom]*len(chunk), starts.tolist(),
                              ends=np.minimum(starts+bin_size, length).tolist(),
                              values=values[chunk].tolist())


def prepare(sample, input_bam, output_root=ROOT):
    begin = time.time()
    params = json.loads((output_root/'config/params.json').read_text())
    blacklist = load_intervals(output_root/'config/hg38-blacklist.v2.bed')
    canonical = params['canonical_chromosomes']
    metrics = Counter()
    fragments = {chrom: (array('I'), array('I')) for chrom in canonical}
    cuts = {chrom: array('I') for chrom in canonical}
    histogram = Counter()
    output_bam = output_root/f'bam/{sample}.bam'
    if output_bam.exists() or (output_root/f'qc/{sample}_preparation.json').exists():
        raise FileExistsError(f'Existing output for {sample}; inspect before rerunning')
    with tempfile.TemporaryDirectory(prefix=f'atac-sept08-{sample}-') as scratch:
        unsorted = Path(scratch)/'accepted.bam'
        with pysam.AlignmentFile(str(input_bam), 'rb', threads=4) as src:
            sizes = {chrom: src.get_reference_length(chrom) for chrom in canonical}
            header = src.header.to_dict()
            header['HD']['SO'] = 'unsorted'
            header.setdefault('PG', []).append({'ID': 'SEPT08_pair_filter', 'PN': '02_prepare.py',
              'VN': '1', 'CL': 'MAPQ>=30 both mates; complete primary proper pairs; no duplicates/QCfail; canonical; fragment blacklist'})
            with pysam.AlignmentFile(str(unsorted), 'wb', header=header, threads=4) as dst:
                for read, mate, chrom, start, end in eligible_pairs(
                        src, set(canonical), blacklist, metrics, params['min_mapq_both_mates']):
                    dst.write(read); dst.write(mate)
                    fragments[chrom][0].append(start); fragments[chrom][1].append(end)
                    histogram[end-start] += 1
                    for alignment in (read, mate):
                        cut = alignment.reference_end-5 if alignment.is_reverse else alignment.reference_start+4
                        if 0 <= cut < sizes[chrom]:
                            cuts[chrom].append(cut)
        if metrics['accepted_fragments'] == 0:
            raise ValueError('No eligible fragments')
        print(sample, 'full BAM stream read;', dict(metrics), flush=True)
        partial = str(output_bam)+'.partial'
        pysam.sort('-@', str(params['threads']), '-m', '512M', '-O', 'BAM',
                   '-T', str(Path(scratch)/'sort'), '-o', partial, str(unsorted))
        pysam.quickcheck(partial)
        Path(partial).replace(output_bam)
        pysam.index('-@', str(params['threads']), str(output_bam))
    payload = {}
    for chrom in canonical:
        start = np.asarray(fragments[chrom][0], dtype=np.uint32)
        end = np.asarray(fragments[chrom][1], dtype=np.uint32)
        order = np.argsort(start, kind='stable')
        fragments[chrom] = (start[order], end[order])
        payload[chrom+'_start'], payload[chrom+'_end'] = fragments[chrom]
        cuts[chrom] = np.sort(np.asarray(cuts[chrom], dtype=np.uint32))
    np.savez_compressed(output_root/f'counts/{sample}_fragments.npz', **payload)
    profile, n_tss, enrichment = tss_profile(cuts, sizes, blacklist, output_root/'config/tss_knownGene.bed')
    np.savetxt(output_root/f'qc/{sample}_tss_profile.tsv',
               np.c_[np.arange(-2000, 2000), profile], delimiter='\t', fmt='%d',
               header='offset\tinsertions', comments='')
    np.savetxt(output_root/f'qc/{sample}_fragment_lengths.tsv',
               np.array(sorted(histogram.items())), delimiter='\t', fmt='%d',
               header='length\tfragments', comments='')
    denominator = sum(len(fragments[c][0]) for c in canonical if c != 'chrY')
    write_coverage(output_root/f'tracks/{sample}.fragment_cpm.bw', fragments, sizes, denominator)
    result = dict(metrics)
    result.update(sample=sample, tss_enrichment=enrichment, n_tss=n_tss,
      nfr_fraction=sum(n for length, n in histogram.items() if length < 100)/metrics['accepted_fragments'],
      chrY_fragments=len(fragments['chrY'][0]), nonY_fragments=denominator,
      chromosome_fragments={c: len(fragments[c][0]) for c in canonical},
      tss_method=params['tss_score'], original_mitochondrial_fraction=None,
      track_method='25 bp mean fragment coverage per million eligible non-Y fragments; uncovered bins omitted',
      elapsed_seconds=round(time.time()-begin, 1))
    (output_root/f'qc/{sample}_preparation.json').write_text(json.dumps(result, indent=2)+'\n')
    print(sample, 'preparation complete', result['elapsed_seconds'], 'seconds', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('sample')
    args = parser.parse_args()
    with open(ROOT/'config/samples.tsv') as handle:
        sample_rows = {r['sample']: r for r in csv.DictReader(handle, delimiter='\t')}
    row = sample_rows[args.sample]
    prepare(args.sample, (ROOT/row['bam']).resolve())
