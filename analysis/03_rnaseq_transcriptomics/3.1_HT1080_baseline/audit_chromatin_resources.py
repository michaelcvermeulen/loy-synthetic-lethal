"""Read-only inventory of earlier local chromatin analyses and their provenance."""
from pathlib import Path
from datetime import datetime, timezone
from collections import Counter
import hashlib
import json
import re
import csv
import gzip
import pandas as pd
from bs4 import BeautifulSoup

ROOT=Path(__file__).resolve().parent
SOURCE=Path('A:/multiome/results/ChIPseq')
OUT=ROOT/'provenance/chipseq_resources'


def main():
    textdir=OUT/'legacy_html_text';textdir.mkdir(parents=True,exist_ok=True)
    files=sorted(p for p in SOURCE.rglob('*') if p.is_file())
    rows=[{'relative_path':p.relative_to(SOURCE).as_posix(),'absolute_path':str(p),
           'bytes':p.stat().st_size,'modified_utc':datetime.fromtimestamp(p.stat().st_mtime,timezone.utc).isoformat(),
           'extension':p.suffix.lower()} for p in files]
    pd.DataFrame(rows).to_csv(OUT/'local_file_inventory.csv',index=False)
    selected=list(SOURCE.glob('*.html'))
    selected += [SOURCE/'ACTIVE_MOTIF_ATAC'/p for p in [
        'atacseq_html/01_consensus_peaks.html','atacseq_html/02_diffbind.html',
        'atacseq_html/04_atac_rna_integration.html','atacseq_html/14_combined_master_table.html',
        'nf_core/01_consensus_peaks.html','nf_core/02_diffbind.html','nf_core/highlights.html']]
    reports=[];accessions=[];paths=[]
    for path in selected:
        if not path.exists():continue
        raw=path.read_bytes(); markup=raw.decode('utf-8',errors='replace')
        markup=re.sub(r'<script\b[^>]*>.*?</script>','',markup,flags=re.I|re.S)
        markup=re.sub(r'<style\b[^>]*>.*?</style>','',markup,flags=re.I|re.S)
        markup=re.sub(r'<img\b[^>]*>','',markup,flags=re.I|re.S)
        soup=BeautifulSoup(markup,'html.parser')
        text=soup.get_text('\n',strip=True)
        relative=path.relative_to(SOURCE).as_posix()
        outname=relative.replace('/','__').replace('.html','.txt')
        (textdir/outname).write_text(text,encoding='utf-8')
        heads=[h.get_text(' ',strip=True) for h in soup.find_all(['h1','h2','h3'])]
        reports.append({'report':relative,'sha256':hashlib.sha256(raw).hexdigest(),
            'text_file':outname,'characters':len(text),'headings':' | '.join(heads)})
        for acc in sorted(set(re.findall(r'\b(?:(?:GSE|GSM|SRR|SRX)\d+|(?:ENCSR|ENCFF)[A-Z0-9]+)\b',text))):
            accessions.append({'report':relative,'accession':acc})
        for match in re.finditer(r'["\']([A-Za-z]:[/\\][^"\'\n]+)["\']',text):
            value=match.group(1);resolved=Path(value.replace('\\','/'))
            paths.append({'report':relative,'referenced_path':value,'exists':resolved.exists(),'is_directory':resolved.is_dir()})
    pd.DataFrame(reports).to_csv(OUT/'legacy_report_inventory.csv',index=False)
    pd.DataFrame(accessions).to_csv(OUT/'legacy_report_accessions.csv',index=False)
    pd.DataFrame(paths).drop_duplicates().to_csv(OUT/'legacy_referenced_paths.csv',index=False)
    result={'audited_utc':datetime.now(timezone.utc).isoformat(),'source':str(SOURCE),'files':len(files),
        'bytes':sum(r['bytes'] for r in rows),'extensions':dict(Counter(r['extension'] for r in rows)),
        'reports_extracted':len(reports),'source_modified':False}
    (OUT/'local_inventory_summary.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))
    print(pd.DataFrame(reports)[['report','characters']].to_string(index=False))


def upstream_inventory():
    source = Path('A:/multiome/data/ChIPseq')
    rows = []
    samples = []
    for p in sorted(source.rglob('*')):
        if not p.is_file():
            continue
        rows.append({'relative_path': p.relative_to(source).as_posix(),
                     'absolute_path': str(p), 'bytes': p.stat().st_size,
                     'extension': p.suffix.lower()})
        if 'series_matrix' not in p.name or not p.name.endswith('.gz'):
            continue
        fields = {}
        with gzip.open(p, 'rt', encoding='utf-8') as handle:
            for line in handle:
                if line.startswith('!Sample_'):
                    values = next(csv.reader([line], delimiter='\t'))
                    fields.setdefault(values[0], []).append(values[1:])
        accessions = fields.get('!Sample_geo_accession', [[]])[0]
        for i, accession in enumerate(accessions):
            row = {'accession': accession, 'metadata_source': str(p)}
            for key, vals in fields.items():
                row[key.removeprefix('!Sample_')] = ' | '.join(v[i] for v in vals if len(v) > i)
            samples.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(OUT/'upstream_file_inventory.csv', index=False)
    pd.DataFrame(samples).to_csv(OUT/'geo_local_sample_metadata.csv', index=False)
    df['dataset'] = df.relative_path.str.split('/').str[0]
    print(df.groupby('dataset').agg(files=('bytes','size'), bytes=('bytes','sum')).to_string())
    print('GEO samples:', len(samples))


if __name__=='__main__':
    import sys
    if '--upstream-only' not in sys.argv:
        main()
    upstream_inventory()
