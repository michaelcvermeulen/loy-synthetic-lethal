"""Download small public metadata records, never local experimental data."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import json
import urllib.request

OUT = Path(__file__).resolve().parent/'provenance/chipseq_resources/public_metadata'
OUT.mkdir(parents=True, exist_ok=True)
URLS = {
    'ENCFF461DZD.json': 'https://www.encodeproject.org/files/ENCFF461DZD/?format=json',
    'PMC13081975.xml': 'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13081975/fullTextXML',
    'PMC12311792.xml': 'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12311792/fullTextXML',
    'ENCSR295VRM.json': 'https://www.encodeproject.org/annotations/ENCSR295VRM/?format=json',
}
for acc in ['ENCFF688MFO', 'ENCFF090WXN', 'ENCFF633FBI', 'ENCFF408PEC',
            'ENCFF827ACX', 'ENCFF544VBI']:
    URLS[acc+'.json'] = f'https://www.encodeproject.org/files/{acc}/?format=json'

def fetch(item):
    name, url = item
    try:
        if (OUT/name).exists():
            return {'file': name, 'url': url, 'bytes': (OUT/name).stat().st_size, 'status': 'cached'}
        request = urllib.request.Request(url, headers={'User-Agent': 'HT1080-chromatin-resource-audit/1.0'})
        with urllib.request.urlopen(request, timeout=25) as handle:
            data = handle.read(10000000)
        (OUT/name).write_bytes(data)
        return {'file': name, 'url': url, 'bytes': len(data), 'status': 'downloaded'}
    except Exception as exc:
        return {'file': name, 'url': url, 'status': 'unavailable', 'error': str(exc)}

if __name__ == '__main__':
    with ThreadPoolExecutor(max_workers=3) as pool:
        rows = list(pool.map(fetch, URLS.items()))
    (OUT/'fetch_status.json').write_text(json.dumps(rows, indent=2), encoding='utf-8')
    print(json.dumps(rows, indent=2))
