from pathlib import Path
import hashlib,json,urllib.request
ROOT=Path(__file__).resolve().parents[1]
out=ROOT/'integration/RNA_interaction_CN/sources';out.mkdir(parents=True,exist_ok=True)
url='https://raw.githubusercontent.com/asntech/QDNAseq.hg38/v1.2.0/data/hg38.100kbp.SR50.rda'
path=out/'hg38.100kbp.SR50.rda'
if not path.exists():
    with urllib.request.urlopen(url,timeout=90) as response:path.write_bytes(response.read())
manifest={'source_url':url,'package':'QDNAseq.hg38','version':'1.2.0','doi':'10.5281/zenodo.17948766','sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'bytes':path.stat().st_size}
(out/'hg38_annotation_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
print(json.dumps(manifest,indent=2))
