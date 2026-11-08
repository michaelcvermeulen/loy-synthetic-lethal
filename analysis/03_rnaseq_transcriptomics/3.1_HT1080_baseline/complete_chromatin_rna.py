"""Check publication artifacts and record completion after scientific and visual review."""
from pathlib import Path
import argparse
import datetime
import hashlib
import json
import platform
import matplotlib
import numpy as np
import pandas as pd
from PIL import Image
from pypdf import PdfReader
from bs4 import BeautifulSoup

ROOT=Path(__file__).resolve().parent
PROV=ROOT/'provenance/chromatin_rna'


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--visual-reviewed',action='store_true')
    args=parser.parse_args()
    scientific=json.loads((PROV/'validation.json').read_text())
    camera=json.loads((PROV/'camera_validation.json').read_text())
    assert scientific['status']==camera['status']=='passed'
    manifest=pd.read_csv(ROOT/'outputs/figures/chromatin_rna/manifest.csv')
    artifacts=[]
    for r in manifest.itertuples():
        for ext in ['pdf','svg','png']:
            p=ROOT/getattr(r,ext)
            sha=hashlib.sha256(p.read_bytes()).hexdigest()
            assert sha==getattr(r,f'{ext}_sha256')
            artifacts.append({'path':str(p.relative_to(ROOT)),'sha256':sha})
        image=Image.open(ROOT/r.png)
        assert all(abs(v-600)<.2 for v in image.info['dpi'])
        assert abs(image.width-r.width_mm/25.4*600)<2
        pdf=PdfReader(ROOT/r.pdf)
        assert len(pdf.pages)==1
        assert abs(float(pdf.pages[0].mediabox.width)-r.width_mm/25.4*72)<1
    report=ROOT/'outputs/report/chromatin_rna_enrichment.html'
    soup=BeautifulSoup(report.read_text(encoding='utf-8'),'html.parser')
    links=0
    for a in soup.find_all('a',href=True):
        href=a['href']
        if href.startswith(('http://','https://','#','mailto:')):
            continue
        assert (report.parent/href.split('#')[0]).resolve().exists(),href
        links+=1
    assert len(soup.find_all('img'))==4
    for r in pd.read_csv(PROV/'input_manifest.csv').itertuples():
        assert hashlib.sha256((ROOT/r.path).read_bytes()).hexdigest()==r.sha256,r.path
    for p in [ROOT/'CHROMATIN_RNA_ENRICHMENT.md',report,ROOT/'outputs/captions/chromatin_rna.md']:
        artifacts.append({'path':str(p.relative_to(ROOT)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
    sources=['run_chromatin_rna.R','run_chromatin_rna.py','summarize_chromatin_rna.py','plot_chromatin_rna.py',
             'validate_chromatin_rna.R','complete_chromatin_rna.py']
    sources=[{'path':p,'sha256':hashlib.sha256((ROOT/p).read_bytes()).hexdigest()} for p in sources]
    prior=json.loads((PROV/'completion.json').read_text()) if (PROV/'completion.json').exists() else {}
    unchanged=prior.get('artifacts')==artifacts and prior.get('visual_review')=='passed'
    visual=args.visual_reviewed or unchanged
    result={'status':'complete' if visual else 'awaiting_visual_review',
            'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
            'scientific_validation':scientific,'camera_validation':camera,
            'publication_figures':len(manifest),'local_report_links_checked':links,
            'visual_review':'passed' if visual else 'pending','artifacts':artifacts,'scripts':sources}
    result['python_environment']={'python':platform.python_version(),'pandas':pd.__version__,
                                  'numpy':np.__version__,'matplotlib':matplotlib.__version__}
    (PROV/'completion.json').write_text(json.dumps(result,indent=2))
    print(json.dumps({k:v for k,v in result.items() if k not in ['artifacts','scripts']}))


if __name__=='__main__':
    main()
