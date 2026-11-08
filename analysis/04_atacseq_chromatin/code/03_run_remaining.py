"""Run disjoint sample jobs with bounded concurrency after inspecting the pilot."""
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def run_sample(sample):
    with open(ROOT/f'logs/{sample}_preparation.log', 'w') as log:
        subprocess.run([sys.executable,'-u',str(ROOT/'code/02_prepare.py'),sample],
                       check=True,stdout=log,stderr=subprocess.STDOUT)
    subprocess.run([sys.executable,'-u',str(ROOT/'code/03_call_peaks.py'),sample],check=True)
    return sample


if __name__ == '__main__':
    qc = json.loads((ROOT/'qc/CTRL_11_preparation.json').read_text())
    assert qc['accepted_fragments'] > 0 and qc['tss_enrichment'] is not None
    with open(ROOT/'config/samples.tsv') as handle:
        samples = [row['sample'] for row in csv.DictReader(handle,delimiter='\t') if row['sample'] != 'CTRL_11']
    subprocess.run([sys.executable,str(ROOT/'code/test_counts.py')],check=True)
    with ThreadPoolExecutor(max_workers=2) as pool:
        jobs = {pool.submit(run_sample,s):s for s in samples}
        for job in as_completed(jobs):
            print('Completed',job.result(),flush=True)
    subprocess.run([sys.executable,'-u',str(ROOT/'code/04_consensus_counts.py')],check=True)
