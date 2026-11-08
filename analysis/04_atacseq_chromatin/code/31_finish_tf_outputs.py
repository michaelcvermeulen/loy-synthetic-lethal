"""Finish figures, cross-method tables, and report after verified TOBIAS completion."""
from pathlib import Path
import subprocess,time,sys,json
ROOT=Path(__file__).resolve().parents[1]
marker=ROOT/'motifs/TF_extension/TOBIAS/ANALYSIS_COMPLETE.json'
deadline=time.monotonic()+6*3600
print('Waiting for the complete TOBIAS comparison',flush=True)
while not marker.exists():
    if time.monotonic()>deadline:raise TimeoutError(marker)
    time.sleep(30)
original_python='/home/micha/.local/share/atacseq-sept08-2026-venv/bin/python'
for script,program in [('28_summarize_tobias.py',sys.executable),
                       ('29_footprint_aggregate.py',sys.executable),
                       ('32_verify_tf_extension.py',sys.executable),
                       ('09_report.py',original_python),('11_verify_outputs.py',original_python)]:
    print('Finalizing',script,flush=True)
    logpath=ROOT/f'logs/31_{Path(script).stem}.log'
    with logpath.open('w') as log:
        result=subprocess.run([program,str(ROOT/'code'/script)],stdout=log,stderr=subprocess.STDOUT)
    if result.returncode:
        print(logpath.read_text()[-6000:],flush=True)
        raise RuntimeError(f'{script} failed: {logpath}')
(ROOT/'motifs/TF_extension/OUTPUTS_COMPLETE.json').write_text(json.dumps(dict(
    calculations_complete=True,report_rebuilt=True,verification_passed=True,
    technical_batch_assignments='unconfirmed',visual_review='separate final review'),indent=2)+'\n')
print('TF figures, integrated tables, report and verification complete',flush=True)
