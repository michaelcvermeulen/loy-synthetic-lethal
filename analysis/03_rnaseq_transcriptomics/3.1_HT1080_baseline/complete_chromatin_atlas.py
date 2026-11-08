"""Finalize output hashes after scientific validation and figure review."""
from pathlib import Path
from datetime import datetime,timezone
import argparse
import ast
import hashlib
import json
import re
import platform
import pandas as pd
from PIL import Image
from pypdf import PdfReader

ROOT=Path(__file__).resolve().parent
PROV=ROOT/'provenance/chromatin_atlas'
FIG=ROOT/'outputs/figures/chromatin_atlas'
OUT=ROOT/'outputs/tables/chromatin_atlas'


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--visual-reviewed',action='store_true')
    args=parser.parse_args()
    for p in [PROV/'validation.json',PROV/'signal_validation.json']:
        assert json.loads(p.read_text())['status']=='passed'
    png=FIG/'HT1080_measured_chromatin_reference.png'
    im=Image.open(png)
    assert min(im.size)>4000 and im.info['dpi'][0]>599
    pdf=PdfReader(FIG/'HT1080_measured_chromatin_reference.pdf')
    assert len(pdf.pages)==1
    for p in [ROOT/'CHROMATIN_REFERENCE_ATLAS.md',ROOT/'outputs/report/chromatin_reference_atlas.html']:
        txt=p.read_text(encoding='utf-8')
        links=re.findall(r'\]\(([^)]+)\)',txt) if p.suffix=='.md' else re.findall(r'href="([^"]+)"',txt)
        for target in links:
            if not target.startswith(('http:','https:','data:','#')):
                assert (p.parent/target).exists(),target
    scripts=[ROOT/f for f in ['prepare_chromatin_atlas.py','summarize_chromatin_signal.py',
              'build_chromatin_peak_features.R','finalize_chromatin_atlas.py','validate_chromatin_atlas.py',
              'plot_chromatin_atlas.py','run_chromatin_atlas.py','complete_chromatin_atlas.py']]
    for p in scripts:
        if p.suffix=='.py':ast.parse(p.read_text(encoding='utf-8'))
    old={}
    if (PROV/'completion.json').exists():old=json.loads((PROV/'completion.json').read_text())
    png_hash=sha(png)
    reviewed=args.visual_reviewed or (old.get('figure_png_sha256')==png_hash and old.get('visual_review')=='passed')
    artifacts=[p for d in [OUT,FIG,PROV] for p in d.rglob('*') if p.is_file() and p.name!='completion.json']
    artifacts+=scripts+[ROOT/'config/chromatin_atlas.json',ROOT/'CHROMATIN_REFERENCE_ATLAS.md',
                         ROOT/'outputs/report/chromatin_reference_atlas.html',
                         ROOT/'outputs/captions/chromatin_reference_atlas.md']
    result={'completed_utc':datetime.now(timezone.utc).isoformat(),
            'status':'complete' if reviewed else 'scientific_checks_passed_visual_review_pending',
            'scientific_validation':'passed','signal_validation':'passed',
            'visual_review':'passed' if reviewed else 'pending','figure_png_sha256':png_hash,
            'png_pixels':im.size,'png_dpi':im.info['dpi'],
            'python_version':platform.python_version(),'pandas_version':pd.__version__,
            'new_enrichment_run':False,'original_inputs_modified':False,
            'artifact_hashes':[{'path':str(p.relative_to(ROOT)),'bytes':p.stat().st_size,'sha256':sha(p)} for p in artifacts]}
    (PROV/'completion.json').write_text(json.dumps(result,indent=2))
    print(json.dumps({k:v for k,v in result.items() if k!='artifact_hashes'},indent=2))


if __name__=='__main__':main()
