"""Check local links, package the portable report, and verify archive extraction."""
from pathlib import Path
from html.parser import HTMLParser
from urllib.parse import urlsplit, unquote
import argparse
import hashlib
import json
import re
import shutil
import zipfile

ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT/'exports/HT1080_ATAC_RNA_portable'
AUDIT = ROOT/'integration/portable_atac_rna'
ARCHIVE = ROOT/'exports/HT1080_ATAC_RNA_portable.zip'
EXTRACTED = AUDIT/'extracted_check/Shared Report/HT1080_ATAC_RNA_portable'


class Links(HTMLParser):
    def __init__(self):
        super().__init__(); self.links=[]
    def handle_starttag(self,tag,attrs):
        self.links += [v for k,v in attrs if k in ('href','src') and v]


def dump(path,value):
    path.write_text(json.dumps(value,indent=2),encoding='utf-8')


def prepare():
    for filename in ('METHODS.txt','selection_summary.json','verification.json'):
        shutil.copy2(ROOT/'integration/gene_atlas'/filename,PACKAGE/'integration/gene_atlas'/filename)
    atlas = PACKAGE/'report/gene_atlas.html'
    text = atlas.read_text(encoding='utf-8')
    text = re.sub(r'<a href="../integration/gene_atlas/gene_diagrams_PDFs.zip">.*?</a>',
                  'All 251 PDF books are included in the integration/gene_atlas/figures folder',text)
    text = text.replace('>Main ATAC report</a>','>Interactive report</a>')
    atlas.write_text(text,encoding='utf-8')
    problems=[]; count=0
    for page in PACKAGE.rglob('*.html'):
        parser=Links();parser.feed(page.read_text(encoding='utf-8'))
        for link in parser.links:
            url=urlsplit(link)
            if url.scheme or url.netloc or not url.path:continue
            target=(page.parent/unquote(url.path)).resolve()
            if not target.is_relative_to(PACKAGE.resolve()) or not target.is_file():
                problems.append({'page':str(page.relative_to(PACKAGE)),'link':link})
            count+=1
    assert not problems,json.dumps(problems,indent=2)
    assert (PACKAGE/'report/interactive_assets/data.js').read_bytes()==(ROOT/'report/interactive_assets/data.js').read_bytes()
    assert not list(PACKAGE.rglob('*.bam')) and not list(PACKAGE.rglob('*fragments.npz'))
    dump(AUDIT/'package_verification.json',{'status':'passed','local_HTML_links_checked':count,'broken_links':0,'scatter_RNA_data_byte_identical':True,'BAMs_or_individual_fragments_included':False})
    verification=PACKAGE/'verification';verification.mkdir(exist_ok=True)
    for filename in ('build_manifest.json','coverage_integral_verification.json','package_verification.json'):
        shutil.copy2(AUDIT/filename,verification/filename)
    files=sorted(p for p in PACKAGE.rglob('*') if p.is_file())
    with zipfile.ZipFile(ARCHIVE,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as archive:
        for p in files:archive.write(p,arcname=str(Path(PACKAGE.name)/p.relative_to(PACKAGE)))
    with zipfile.ZipFile(ARCHIVE) as archive:
        assert archive.testzip() is None
        # Only locally generated, path-validated archive entries are extracted.
        extraction_root=EXTRACTED.parent.resolve()
        for entry in archive.infolist():
            assert (extraction_root/entry.filename).resolve().is_relative_to(extraction_root)
        archive.extractall(extraction_root)
    for p in files:
        q=EXTRACTED/p.relative_to(PACKAGE)
        assert hashlib.sha256(p.read_bytes()).digest()==hashlib.sha256(q.read_bytes()).digest()
    result={'status':'passed','files':len(files),'zip_bytes':ARCHIVE.stat().st_size,'extracted_bytes':sum(p.stat().st_size for p in files),'zip_CRC':'passed','extracted_checksums':'all_match','extracted_path':str(EXTRACTED)}
    dump(AUDIT/'archive_verification.json',result)
    print(json.dumps(result,indent=2),flush=True)


def finalize():
    browser=json.loads((AUDIT/'browser_verification.json').read_text(encoding='utf-8'))
    assert browser['status']=='passed'
    assert browser['portable_validation']['HTTP_requests']==0
    assert Path(browser['portable_validation']['package_root']).resolve()==EXTRACTED.resolve()
    manifest=json.loads((AUDIT/'build_manifest.json').read_text())
    manifest['status']='complete'
    dump(AUDIT/'build_manifest.json',manifest)
    # Only verification documents are added/updated after testing the extracted report.
    with zipfile.ZipFile(ARCHIVE) as archive:
        keep=[(i,archive.read(i.filename)) for i in archive.infolist() if i.filename==PACKAGE.name+'/verification/build_manifest.json']
    assert len(keep)==1
    # Avoid duplicate ZIP members: rebuild from the verified package with new verification records.
    for filename in ('build_manifest.json','browser_verification.json','archive_verification.json'):
        shutil.copy2(AUDIT/filename,PACKAGE/'verification'/filename)
    files=sorted(p for p in PACKAGE.rglob('*') if p.is_file())
    with zipfile.ZipFile(ARCHIVE,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as archive:
        for p in files:archive.write(p,arcname=str(Path(PACKAGE.name)/p.relative_to(PACKAGE)))
    with zipfile.ZipFile(ARCHIVE) as archive:
        assert archive.testzip() is None
        for p in files:
            assert hashlib.sha256(archive.read(str(Path(PACKAGE.name)/p.relative_to(PACKAGE)).replace('\\','/'))).digest()==hashlib.sha256(p.read_bytes()).digest()
    sha=hashlib.sha256(ARCHIVE.read_bytes()).hexdigest()
    ARCHIVE.with_suffix('.zip.sha256').write_text(sha+'  '+ARCHIVE.name+'\n')
    result={'status':'complete','archive':str(ARCHIVE),'bytes':ARCHIVE.stat().st_size,'files':len(files),'sha256':sha,'genes':manifest['genes'],'panels':manifest['panels'],'offline_browser_verification':'passed','tested_extracted_package':str(EXTRACTED)}
    dump(AUDIT/'ANALYSIS_COMPLETE.json',result)
    print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--finalize',action='store_true');args=parser.parse_args()
    finalize() if args.finalize else prepare()
