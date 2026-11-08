from pathlib import Path
import hashlib,json,csv
from datetime import datetime,timezone
out=Path(__file__).resolve().parent
root=out.parents[3]
rows=[]
with (out/'protected_files_before.csv').open(encoding='utf-8-sig',newline='') as handle:
    for row in csv.DictReader(handle):
        p=root/row['path']
        digest=hashlib.md5(p.read_bytes()).hexdigest()
        rows.append({'path':row['path'],'before_md5':row['md5'],'final_md5':digest,
            'unchanged':digest==row['md5'],'last_modified_utc':datetime.fromtimestamp(p.stat().st_mtime,timezone.utc).isoformat()})
changed=[r for r in rows if not r['unchanged']]
assert all(r['path']=='SEPT_11_2026_THESIS/draftv2/chapter3/chapter3_results.md' for r in changed)
(out/'final_protected_file_audit.json').write_text(json.dumps({'checked_utc':datetime.now(timezone.utc).isoformat(),'files':rows,'changes_after_rerun':changed},indent=2))
if changed:
    print('Results Markdown changed separately after the rerun; all other protected files unchanged.')
else:
    print('All protected files unchanged.')
