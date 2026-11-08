"""Run four deterministic partitions of the frozen analysis, merge, and validate."""
from pathlib import Path
import argparse
import hashlib
import json
import subprocess
import time
import pandas as pd

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'outputs/tables/chromatin_rna'
PROV=ROOT/'provenance/chromatin_rna'
RSCRIPT='C:/Program Files/R/R-4.5.2/bin/Rscript.exe'


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--statistics-only',action='store_true')
    args=parser.parse_args()
    PROV.mkdir(parents=True,exist_ok=True)
    if not (PROV/'input_manifest.csv').exists():
        subprocess.run([RSCRIPT,'--vanilla',str(ROOT/'run_chromatin_rna.R'),'--freeze-only'],cwd=ROOT,check=True)
    frozen=pd.read_csv(PROV/'input_manifest.csv')
    for r in frozen.itertuples():
        assert hashlib.sha256((ROOT/r.path).read_bytes()).hexdigest()==r.sha256,r.path
    workers=[]
    for i in range(1,5):
        done=OUT/f'partitions/part_{i}/completion.json'
        if done.exists() and json.loads(done.read_text())['status']=='partition_complete':
            print(f'Partition {i}: previously completed',flush=True)
            continue
        log=(PROV/f'partition_{i}.log').open('w',encoding='utf-8')
        p=subprocess.Popen([RSCRIPT,'--vanilla',str(ROOT/'run_chromatin_rna.R'),f'--partition={i}'],
                           cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,
                           creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0))
        workers.append((i,p,log))
        print(f'Started partition {i}, PID {p.pid}',flush=True)
    while workers:
        for i,p,log in workers[:]:
            code=p.poll()
            if code is not None:
                log.close();workers.remove((i,p,log))
                print(f'Partition {i} finished with exit {code}',flush=True)
                if code:
                    raise RuntimeError(f'Partition {i} failed; inspect its log. Other authorized workers may still finish.')
        if workers:
            time.sleep(5)
    for name in ['fgsea_results','camera_results']:
        pieces=[pd.read_csv(OUT/f'partitions/part_{i}/{name}.csv') for i in range(1,5)]
        combined=pd.concat(pieces,ignore_index=True)
        keys=['scenario','contrast','pathway']+(['correlation_mode'] if name.startswith('camera') else [])
        assert not combined.duplicated(keys).any()
        combined.to_csv(OUT/f'{name}.csv',index=False)
    # Serial finalization reuses every completed result, rechecks frozen inputs and writes pooled sensitivity q values.
    subprocess.run([RSCRIPT,'--vanilla',str(ROOT/'run_chromatin_rna.R')],cwd=ROOT,check=True)
    subprocess.run(['python','-X','utf8',str(ROOT/'summarize_chromatin_rna.py')],cwd=ROOT,check=True)
    subprocess.run([RSCRIPT,'--vanilla',str(ROOT/'validate_chromatin_rna.R')],cwd=ROOT,check=True)
    if not args.statistics_only:
        subprocess.run(['python','-X','utf8',str(ROOT/'plot_chromatin_rna.py')],cwd=ROOT,check=True)
        subprocess.run(['python','-X','utf8',str(ROOT/'complete_chromatin_rna.py')],cwd=ROOT,check=True)
    print('Chromatin RNA analysis and checks complete.',flush=True)


if __name__=='__main__':
    main()
