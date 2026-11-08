from pathlib import Path
import json,re
ROOT=Path(__file__).resolve().parent
for sample in ['HT1080-sgCTRL-D0','HT1080-LOY-D0','LNCaP-sgCTRL','LNCaP-LOY-2X']:
    qc=ROOT/'qc';log=ROOT/'logs'/f'{sample}.pipeline.log';align=ROOT/'logs'/f'{sample}.alignment.log'
    stage='queued';total=None;done=None;last=''
    if log.exists():
        lines=log.read_text(errors='replace').splitlines();last=lines[-1] if lines else ''
        stage='read recovery/QC'
    original=qc/f'{sample}.original_read_qc.json'
    if original.exists():total=json.loads(original.read_text())['counts']['primary']
    if align.exists():
        stage='aligning';done=sum(map(int,re.findall(r'Processed (\d+) reads',align.read_text(errors='replace'))))
    if (qc/f'{sample}.ALIGNMENT_FILE_COMPLETE.txt').exists():stage='validating aligned reads'
    if (qc/f'{sample}.REALIGNMENT_COMPLETE.txt').exists():stage='complete'
    print(json.dumps(dict(sample=sample,stage=stage,processed_read_ends=done,total_primary_reads=total,last_pipeline_message=last)))
