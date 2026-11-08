"""Freeze primary bibliographic records and available full text for source interpretation."""
from pathlib import Path
import csv
import hashlib
import json
import requests
from concurrent.futures import ThreadPoolExecutor

root = Path(__file__).resolve().parent
dest = root / 'resources/epigenetic_deep_dive'
pmids = ['16618801', '17603471', '18509334', '18600261', '17339329',
         '18413731', '18443585', '18488029', '17452456', '11137995',
         '20713442', '16135788', '17709396', '18172295']

def retrieve(pmid):
    path = dest / f'PMID_{pmid}.json'
    url = 'https://www.ebi.ac.uk/europepmc/webservices/rest/search'
    if not path.exists():
        response = requests.get(url, params={'query': f'EXT_ID:{pmid} AND SRC:MED',
                                'format': 'json', 'resultType': 'core'}, timeout=45)
        response.raise_for_status()
        result = response.json()['resultList']['result']
        assert len(result) == 1 and result[0]['id'] == pmid
        path.write_text(json.dumps(result[0], indent=2), encoding='utf-8')
    record = json.loads(path.read_text(encoding='utf-8'))
    pmc = record.get('pmcid', '')
    xml = dest / f'{pmc}.xml'
    status = 'not_available'
    if pmc:
        if not xml.exists():
            try:
                response = requests.get(f'https://www.ebi.ac.uk/europepmc/webservices/rest/{pmc}/fullTextXML', timeout=45)
                if response.ok and '<article' in response.text:
                    xml.write_bytes(response.content)
            except requests.RequestException:
                pass
        if xml.exists():
            status = 'full_text_cached'
    return {'pmid': pmid, 'title': record['title'], 'year': record.get('pubYear'),
            'doi': record.get('doi', ''), 'pmcid': pmc, 'full_text_status': status,
            'primary_url': f'https://pubmed.ncbi.nlm.nih.gov/{pmid}/'}

if __name__ == '__main__':
    with ThreadPoolExecutor(max_workers=4) as pool:
        rows = list(pool.map(retrieve, pmids))
    with (dest / 'primary_source_index.csv').open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    paths = sorted(p for p in dest.iterdir() if p.is_file() and p.name != 'manifest.csv')
    with (dest / 'manifest.csv').open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=['path', 'sha256'])
        writer.writeheader()
        writer.writerows({'path': p.relative_to(root).as_posix(),
                          'sha256': hashlib.sha256(p.read_bytes()).hexdigest()} for p in paths)
    print('\n'.join(f"{r['pmid']} {r['full_text_status']} {r['title']}" for r in rows))
