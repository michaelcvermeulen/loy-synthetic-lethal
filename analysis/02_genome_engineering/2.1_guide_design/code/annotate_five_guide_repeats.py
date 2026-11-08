"""Annotate exact 20-nt Y matches using published T2T repeat BED tracks.

Run from projects: python sgRNA_guide_development/code/annotate_five_guide_repeats.py
Input coordinates are 1-based; BED intervals are zero-based, half-open.
"""
import csv
import hashlib
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT.parents[1] / 'data'
OUT = ROOT / 'five_guide_deep_dive'


def read_csv(path):
    with path.open(newline='', encoding='utf-8-sig') as handle:
        return list(csv.DictReader(handle))


def write_csv(name, rows):
    with (OUT / name).open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def read_bed(name):
    rows = []
    for line in (OUT / name).read_text().splitlines():
        fields = line.split('\t')
        if fields[0] == 'chrY':
            rows.append((int(fields[1]), int(fields[2]), fields[3]))
    return rows


def overlap(start, end, bed):
    return ';'.join(sorted({name for lo, hi, name in bed
                            if lo < end and start < hi})) or 'unannotated'


def rc(sequence):
    return sequence.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


labels = read_csv(ROOT / 'guide_labels.csv')
hits = {row['guide']: row for row in read_csv(DATA / 'top_hits.csv')}
fasta = DATA / 'chm13v2.0_chrY.fa'
with fasta.open() as handle:
    header = next(handle).strip()
    genome = ''.join(line.strip() for line in handle).upper()
assert 'CP086569.2' in header and len(genome) == 62460029
censat = read_bed('chm13v2.0_censat_v2.1.bed')
composite = read_bed('chm13v2.0_composite-repeats_2022DEC.bed')
site_rows, summary = [], []
for label in labels:
    guide = label['guide']
    row = hits[guide]
    positions = sorted(int(x.strip()) for x in row['location_of_exact_hits_Y'].split(',') if x.strip())
    assert len(positions) == int(row['chrY_exact_hits'])
    records = []
    # Independently rescan both orientations to ensure coordinate lists are complete.
    scanned = []
    for motif in (guide, rc(guide)):
        pos = genome.find(motif)
        while pos >= 0:
            scanned.append(pos + 1)
            pos = genome.find(motif, pos + 1)
    assert sorted(scanned) == positions, label['guide_id']
    for pos in positions:
        start, end = pos - 1, pos - 1 + 20
        sequence = genome[start:end]
        assert sequence in (guide, rc(guide))
        strand = '+' if sequence == guide else '-'
        pam = genome[end:end+3] if strand == '+' else rc(genome[start-3:start])
        record = dict(guide_id=label['guide_id'], guide=guide, start_1based=pos,
                      end_1based=end, strand=strand, pam=pam, canonical_NGG=pam.endswith('GG'),
                      censat=overlap(start, end, censat),
                      composite_repeat=overlap(start, end, composite))
        records.append(record)
    site_rows.extend(records)
    counts = Counter((r['censat'], r['composite_repeat']) for r in records)
    for (sat, comp), count in counts.items():
        group = [r for r in records if (r['censat'], r['composite_repeat']) == (sat, comp)]
        summary.append(dict(guide_id=label['guide_id'], guide=guide, exact_matches=count,
                            NGG_matches=sum(r['canonical_NGG'] for r in group),
                            start_min=min(r['start_1based'] for r in group),
                            start_max=max(r['start_1based'] for r in group),
                            censat=sat, composite_repeat=comp))
    print(label['guide_id'], len(records), 'exact;', sum(r['canonical_NGG'] for r in records), 'NGG')
    print('  CenSat:', dict(Counter(r['censat'] for r in records)))
    print('  Composite:', dict(Counter(r['composite_repeat'] for r in records)))
write_csv('exact_Y_repeat_sites.csv', site_rows)
write_csv('repeat_annotation_summary.csv', summary)
paths = [ROOT / 'guide_labels.csv', DATA / 'top_hits.csv', fasta,
         OUT / 'chm13v2.0_censat_v2.1.bed', OUT / 'chm13v2.0_composite-repeats_2022DEC.bed']
write_csv('repeat_annotation_sources.csv', [dict(path=str(p), sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in paths])
