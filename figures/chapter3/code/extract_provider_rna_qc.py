"""Extract provider RNA QC using the author-accepted May sample mapping."""
from pathlib import Path
import csv
import io
import re
import zipfile
import hashlib
import json
from pypdf import PdfReader
import openpyxl

PROJECTS = Path(__file__).resolve().parents[4]
OUT = PROJECTS/'SEPT_11_2026_THESIS/figures/chapter3/source_data/provider_rna_qc'
OUT.mkdir(parents=True,exist_ok=True)
rows=[]
sources=[]
def record(p):
    sources.append({'path':str(p.relative_to(PROJECTS)), 'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
p=PROJECTS/'RNAseq/analysis/A549/QC_reprocess_2026/_rin/Sample QC report of_30-1192550804_250619054238.pdf'
record(p)
s=' '.join(' '.join(page.extract_text().split()) for page in PdfReader(p).pages)
hits=list(re.finditer(r'(A549|HT1080)-\s*(LOY|CTRL|PARENTAL)-\s*(\d)\s+(SHE\d+)',s))
assert len(hits)==20
for i,m in enumerate(hits):
    segment=s[m.end():hits[i+1].start() if i+1<len(hits) else s.index('QC Note Legend')]
    # Values after the final NA are the integrity score and DV200.
    vals=re.findall(r'NA\s+([\d.]+)\s+([\d.]+)',segment)
    integrity,size=map(float,vals[-1])
    rows.append(dict(qc_id=f'{m[1]}-{m[2]}-{m[3]}',provider_id=m[4],group=m[1],integrity_type='RINe',integrity=integrity,size_metric='DV200',size_percent=size,mapping_status='sample name in provider summary',source=p.name))
p=PROJECTS/'RNAseq/analysis/A549/QC_reprocess_2026/_rin/SampleQC_1-20_30-1192550804.pdf'
record(p)
traces=' '.join(page.extract_text() for page in PdfReader(p).pages)
assert 'RINe' in traces
p=PROJECTS/'thesis/Methods/provider_RNA_QC/30-1301012049_qc_report.xlsx'
record(p)
w=openpyxl.load_workbook(p,data_only=True)
for sample,species,conc,integrity,size,status,action in list(w['QC Metrics'].values)[1:]:
    rows.append(dict(qc_id=sample.replace('LNCaP','DU145'),provider_id='',group='DU145',integrity_type='RIN (reported)',integrity=float(integrity),size_metric='DV200',size_percent=float(size),mapping_status='author-confirmed LNCaP to DU145 naming',source=p.name))
p=PROJECTS/'thesis/Methods/provider_RNA_QC/30-1328169507-09-25-2026_164628.zip'
record(p)
with zipfile.ZipFile(p) as z:
    name='30-1328169507_FA_Traces_18MAY2026.pdf'
    doc=PdfReader(io.BytesIO(z.read(name)))
    assert len(doc.pages)==18
    for page in doc.pages:
        t=page.extract_text()
        sid=re.search(r'Sample:\s*(\S+)',t)[1]
        integrity=float(re.search(r'RQN\s*([\d.]+)',t)[1])
        size=float(re.search(r'Smear Analysis 200 nt to 10000 nt\s+[\d.]+ ng/ul\s+([\d.]+) %Total',t)[1])
        rows.append(dict(qc_id='',provider_id=sid,group='HT1080 MAY 2026',integrity_type='RQN',integrity=integrity,size_metric='200-10000 nt fraction',size_percent=size,mapping_status='pending provider-to-library mapping',source=name))
numbered=PROJECTS/'RNAseq/MAY_26_2026_sample_sheet.csv'
relabel=PROJECTS/'RNAseq/counts/STAR_v2_PARfix/MAY_2026_HT1080_ADDBACK/sample_relabeling_2026-05-26.csv'
submission=PROJECTS/'thesis/Methods/provider_RNA_QC/20260925160351_NGSSSF.zip'
for p in [numbered,relabel,submission]:record(p)
with numbered.open(encoding='utf-8-sig',newline='') as f:
    original_to_arm={r['old']:r['new'] for r in csv.DictReader(f)}
with relabel.open(encoding='utf-8-sig',newline='') as f:
    arm_to_final={r['old_filename'].removesuffix('.txt'):r['new_filename'].removesuffix('.txt') for r in csv.DictReader(f)}
for r in rows:
    if r['group']=='HT1080 MAY 2026':
        ordinal=int(r['provider_id'][3:])-795+1
        old=f'HT1080-sgCTRL-{ordinal}'
        arm=original_to_arm[old]
        r['qc_id']='AB-'+arm_to_final.get(arm,arm)
        r['mapping_status']='Author accepted 2026-09-25: assume consecutive provider ID to submission row; numbered-to-arm mapping and expression relabel documented'
assert len(rows)==44
assert len({r['qc_id'] for r in rows})==44
assert all(0<=r['integrity']<=10 and 0<=r['size_percent']<=100 for r in rows)
with (OUT/'provider_qc.tsv').open('w',encoding='utf-8',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');writer.writeheader();writer.writerows(rows)
(OUT/'provenance.json').write_text(json.dumps(sources,indent=2),encoding='utf-8')
print('44 provider records extracted; 26 sample-name matches, 18 author-accepted provider-order matches.')
