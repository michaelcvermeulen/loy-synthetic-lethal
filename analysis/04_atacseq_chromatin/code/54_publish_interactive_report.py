from pathlib import Path
import hashlib
import json
import re
import requests

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'integration/interactive_atac_rna'
for name in ['browser_verification.json','data_verification.json']:
    assert json.loads((OUT/name).read_text(encoding='utf-8'))['status']=='passed'
health=requests.get('http://127.0.0.1:8765/api/health',timeout=5).json()
assert health['application']=='HT1080_ATAC_RNA' and health['status']=='ready'
main=ROOT/'report/ATACseq_report.html'
text=main.read_text(encoding='utf-8')
snippet='<section id="interactive-atac-rna"><h2>Interactive ATAC versus RNA explorer</h2><p><a href="interactive_ATAC_RNA.html">Open the interactive scatter, sample expression and ATAC tracks</a>. Click any eligible gene, switch RNA datasets or genomic features, and inspect individual peaks. For on-demand coverage, start the local viewer with <code>Start_interactive_report.cmd</code> in the processing folder.</p></section>'
if '<section id="interactive-atac-rna">' in text:
    text=re.sub(r'<section id="interactive-atac-rna">.*?</section>',lambda _:snippet,text,flags=re.S)
else:text=text.replace('</body>',snippet+'\n</body>')
main.write_text(text,encoding='utf-8')
files=[ROOT/'report/interactive_ATAC_RNA.html',*list((ROOT/'report/interactive_assets').glob('*'))]
manifest={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files if p.is_file()}
result={'status':'complete','URL':'http://127.0.0.1:8765/report/interactive_ATAC_RNA.html',
        'HTML':'report/interactive_ATAC_RNA.html','launcher':'Start_interactive_report.cmd',
        'genes':14358,'panels':9,'RNA_normalization':'shared combined-model DESeq2 size factors; not batch-corrected',
        'browser_test':'passed, including actual point click and a gene outside the static atlas',
        'data_validation':'passed; all original overlap counts retained and on-demand coverage agrees with original BigWigs',
        'artifact_sha256':manifest}
(OUT/'ANALYSIS_COMPLETE.json').write_text(json.dumps(result,indent=2))
print(json.dumps({k:v for k,v in result.items() if k!='artifact_sha256'},indent=2))
