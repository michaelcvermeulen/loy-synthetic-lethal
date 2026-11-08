"""Retrieve public GSE193885 counts and source documents; keep raw downloads immutable."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import argparse
import gzip
import hashlib
import json
import re
import requests
import pandas as pd
from bs4 import BeautifulSoup

ROOT = Path(__file__).resolve().parent
CACHE = ROOT / 'resources/hdac12_external/GSE193885'
BASE = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE193nnn/GSE193885/'
SOURCES = [
    ('GSE193885_family.soft.gz', BASE + 'soft/GSE193885_family.soft.gz'),
    ('GSE193885_counts_raw.tsv.gz', BASE + 'suppl/GSE193885_rnaseq_deseq_global_counts_raw.tsv.gz'),
    ('hess_2022.html', 'https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1010376'),
    ('hess_2022.pdf', 'https://journals.plos.org/plosgenetics/article/file?id=10.1371/journal.pgen.1010376&type=printable'),
    ('hess_2022_S3_RNAseq.xlsb', 'https://journals.plos.org/plosgenetics/article/file?id=10.1371/journal.pgen.1010376.s008&type=supplementary'),
]
SOURCES += [(accession + '.xml', 'https://www.ebi.ac.uk/ena/browser/api/xml/' + accession)
            for accession in ['SRX13821339', 'SRX13821329', 'SRX13821332',
                              'SAMN25064695', 'SAMN25064719', 'SAMN25064716']]

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def fetch(item):
    name, url = item
    path = CACHE / name
    meta = CACHE / (name + '.retrieval.json')
    if path.exists() and meta.exists():
        old = json.loads(meta.read_text())
        if old.get('status') == 'retrieved':
            assert sha(path) == old['sha256'], f'Cached resource changed: {name}'
            return old
    row = dict(path=path.relative_to(ROOT).as_posix(), url=url,
               retrieved_at=datetime.now(timezone.utc).isoformat())
    try:
        r = requests.get(url, timeout=(15, 55), headers={'User-Agent': 'HT1080-HDAC-reference-analysis/1.0'})
        row.update(http_status=r.status_code, final_url=r.url, content_type=r.headers.get('Content-Type', ''))
        r.raise_for_status()
        assert not any(t in r.content[:10000] for t in [b'reCAPTCHA', b'Checking your browser']), 'Browser challenge'
        if name.endswith('.gz'):
            assert r.content.startswith(b'\x1f\x8b'), 'Expected gzip data'
            gzip.decompress(r.content)
        if name.endswith('.pdf'):
            assert r.content.startswith(b'%PDF'), 'Expected PDF'
        if name.endswith('.xlsb'):
            assert r.content.startswith(b'PK'), 'Expected workbook archive'
        path.write_bytes(r.content)
        row.update(status='retrieved', bytes=path.stat().st_size, sha256=sha(path))
    except Exception as exc:
        row.update(status='unavailable', error=str(exc))
    meta.write_text(json.dumps(row, indent=2), encoding='utf-8')
    print(name, row['status'], row.get('bytes', row.get('error', '')), flush=True)
    return row

def parse_soft():
    records = []
    current = None
    with gzip.open(CACHE / 'GSE193885_family.soft.gz', 'rt', encoding='utf-8') as stream:
        for line in stream:
            if line.startswith('^SAMPLE = '):
                if current:
                    records.append(current)
                current = {'GSM': line.strip().split(' = ', 1)[1]}
            elif current is not None and line.startswith('!Sample_') and ' = ' in line:
                key, value = line.rstrip().split(' = ', 1)
                key = key.removeprefix('!Sample_')
                current.setdefault(key, []).append(value)
    if current:
        records.append(current)
    samples = pd.DataFrame([{k: '; '.join(v) if isinstance(v, list) else v for k, v in r.items()} for r in records])
    assert len(samples) == 69 and samples.GSM.nunique() == 69
    samples.to_csv(ROOT / 'outputs/tables/HDACEXT_GEO_sample_metadata.csv', index=False)
    (CACHE / 'sample_metadata.json').write_text(json.dumps(records, indent=2), encoding='utf-8')
    counts = pd.read_csv(CACHE / 'GSE193885_counts_raw.tsv.gz', sep='\t', nrows=5)
    print('Count columns:', list(counts.columns), flush=True)
    print(samples[['GSM', 'title', 'characteristics_ch1']].to_string(index=False), flush=True)
    html = CACHE / 'hess_2022.html'
    if html.exists():
        soup = BeautifulSoup(html.read_bytes(), 'html.parser')
        (CACHE / 'hess_2022_text.txt').write_text(soup.get_text('\n', strip=True), encoding='utf-8')
    files = [p for p in CACHE.iterdir() if p.is_file() and p.name != 'manifest.csv']
    pd.DataFrame([dict(path=p.relative_to(ROOT).as_posix(), sha256=sha(p), bytes=p.stat().st_size)
                  for p in sorted(files)]).to_csv(CACHE / 'manifest.csv', index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--fetch', action='store_true')
    args = parser.parse_args()
    CACHE.mkdir(parents=True, exist_ok=True)
    if args.fetch:
        with ThreadPoolExecutor(max_workers=4) as pool:
            list(pool.map(fetch, SOURCES))
    parse_soft()
