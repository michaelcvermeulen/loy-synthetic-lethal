"""Allow independent condition correction to overlap once the LOY pool is verified."""
from pathlib import Path
import time, subprocess, sys, os
ROOT=Path(__file__).resolve().parents[1]
marker=ROOT/'motifs/TF_extension/TOBIAS/LOY_matched_pool.complete.json'
deadline=time.monotonic()+4*3600
while not marker.exists():
    if time.monotonic()>deadline:raise TimeoutError(marker)
    time.sleep(15)
subprocess.run([sys.executable,str(ROOT/'code/23_run_tobias.py'),'LOY'],check=True,
               env={**os.environ,'TF_TOBIAS_CORES':'4'})
