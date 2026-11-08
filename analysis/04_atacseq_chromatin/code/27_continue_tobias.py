"""Start CTRL correction as soon as its verified pool is ready, then finish both conditions."""
from pathlib import Path
import subprocess,time,sys
ROOT=Path(__file__).resolve().parents[1];out=ROOT/'motifs/TF_extension/TOBIAS'
deadline=time.monotonic()+4*3600
def wait_for(path):
    print('Waiting for',path.name,flush=True)
    while not path.exists():
        if time.monotonic()>deadline:raise TimeoutError(path)
        time.sleep(15)
wait_for(out/'CTRL_matched_pool.complete.json')
subprocess.run([sys.executable,str(ROOT/'code/23_run_tobias.py'),'CTRL'],check=True)
wait_for(out/'POOL_COMPLETE.json')
subprocess.run([sys.executable,str(ROOT/'code/23_run_tobias.py')],check=True)
