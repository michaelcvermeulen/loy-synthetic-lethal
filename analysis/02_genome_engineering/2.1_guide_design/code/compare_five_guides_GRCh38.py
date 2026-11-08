"""Compare exact guide matches on assembled chrY; excludes unplaced scaffolds."""
import csv
import gzip
import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'five_guide_deep_dive'
REFERENCE = OUT / 'GRCh38_UCSC_chrY.fa.gz'
URL = 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chrY.fa.gz'
with gzip.open(REFERENCE, 'rt') as handle:
    sequence = ''.join(line.strip() for line in handle if not line.startswith('>')).upper()
assert len(sequence) == 57227415
with (ROOT / 'guide_labels.csv').open() as handle:
    labels = list(csv.DictReader(handle))

def rc(s):
    return s.translate(str.maketrans('ACGT', 'TGCA'))[::-1]

sites, summary = [], []
for row in labels:
    hits = []
    for strand, motif in [('+', row['guide']), ('-', rc(row['guide']))]:
        pos = sequence.find(motif)
        while pos >= 0:
            pam = sequence[pos+20:pos+23] if strand == '+' else rc(sequence[pos-3:pos])
            hits.append(dict(guide_id=row['guide_id'], start_1based=pos+1,
                             strand=strand, pam=pam, NGG=pam.endswith('GG')))
            pos = sequence.find(motif, pos+1)
    sites.extend(hits)
    summary.append(dict(guide_id=row['guide_id'], guide=row['guide'],
                        exact_chrY_matches=len(hits), NGG_chrY_matches=sum(h['NGG'] for h in hits),
                        reference_url=URL, reference_sha256=hashlib.sha256(REFERENCE.read_bytes()).hexdigest()))
    print(row['guide_id'], 'exact', len(hits), 'NGG', sum(h['NGG'] for h in hits), flush=True)
for name, rows in [('GRCh38_exact_Y_sites.csv', sites), ('GRCh38_exact_Y_summary.csv', summary)]:
    with (OUT / name).open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
