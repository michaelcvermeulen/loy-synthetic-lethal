"""Call fresh sample peaks with pinned MACS3 and preserve the exact command."""
import argparse
import json
from pathlib import Path
import shlex
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser()
parser.add_argument('sample')
args = parser.parse_args()
params = json.loads((ROOT/'config/params.json').read_text())
qc = json.loads((ROOT/f'qc/{args.sample}_preparation.json').read_text())
if qc['accepted_fragments'] <= 0 or qc['tss_enrichment'] is None:
    raise ValueError('Preparation/QC is incomplete')
output = ROOT/'peaks'/args.sample
output.mkdir(exist_ok=True)
name = args.sample
if (output/f'{name}_peaks.narrowPeak').exists():
    raise FileExistsError('Peak output already exists; inspect before rerunning')
command = [str(Path(sys.executable).parent/'macs3'), 'callpeak',
  '-t', str(ROOT/f'bam/{name}.bam'), '-f', 'BAMPE',
  '-g', str(params['effective_genome_size']), '-n', name,
  '--outdir', str(output), '-q', str(params['peak_qvalue']),
  '--keep-dup', 'all', '--call-summits']
(ROOT/f'logs/{name}_macs3_command.txt').write_text(shlex.join(command)+'\n')
with open(ROOT/f'logs/{name}_macs3.log', 'w') as log:
    subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)
with open(output/f'{name}_peaks.narrowPeak') as handle:
    count = sum(1 for line in handle if line.strip() and not line.startswith('#'))
if count == 0:
    raise ValueError('MACS3 returned no peaks')
print(name, count, 'fresh peaks', flush=True)
