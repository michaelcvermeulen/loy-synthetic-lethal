"""One-time restart of the early LOY stage to improve total worker utilization."""
from pathlib import Path
import os, signal, subprocess, time, json, sys
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'motifs/TF_extension/TOBIAS'
assert not (OUT/'23_ATACorrect_LOY.complete.json').exists()
lines=subprocess.check_output(['ps','-eo','pid=,ppid=,args='],text=True).splitlines()
procs={}
for line in lines:
    pid,parent,args=line.strip().split(None,2);procs[int(pid)]=(int(parent),args)
targets={p for p,(_,args) in procs.items() if str(ROOT/'code/30_start_loy_footprinting.py') in args}
targets.discard(os.getpid());assert len(targets)==1,targets
while True:
    children={p for p,(parent,_) in procs.items() if parent in targets}
    if children.issubset(targets):break
    targets|=children
assert any('ATACorrect' in procs[p][1] and '--prefix LOY' in procs[p][1] for p in targets)
for p in sorted(targets,reverse=True):
    try:os.kill(p,signal.SIGTERM)
    except ProcessLookupError:pass
time.sleep(3)
for p in targets:
    path=Path(f'/proc/{p}/stat')
    assert not path.exists() or path.read_text().split()[2]=='Z',f'Process {p} remains active'
src=(OUT/'LOY').resolve();dst=(OUT/'LOY_interrupted_worker_reallocation').resolve()
assert src.parent==OUT.resolve() and dst.parent==OUT.resolve() and not dst.exists()
src.rename(dst)
for p in [ROOT/'logs/23_ATACorrect_LOY.log',OUT/'23_ATACorrect_LOY.command.json']:
    if p.exists():p.rename(dst/p.name)
(dst/'INTERRUPTED.json').write_text(json.dumps(dict(reason='Restart before read correction with 12 workers instead of 4; scientific settings unchanged',stopped_pids=sorted(targets)),indent=2)+'\n')
print('Restarting LOY with 12 workers',flush=True)
subprocess.run([sys.executable,str(ROOT/'code/23_run_tobias.py'),'LOY'],check=True,
               env={**os.environ,'TF_TOBIAS_CORES':'12'})
