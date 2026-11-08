"""Compare peak discovery at exactly the same eligible fragment depth."""
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def call_matched(sample, depth, seed):
    config = json.loads((ROOT/'config/params.json').read_text())
    fragments = np.load(ROOT/f'counts/{sample}_fragments.npz')
    chroms = config['canonical_chromosomes']
    sizes = np.array([len(fragments[c+'_start']) for c in chroms], dtype=np.int64)
    rng = np.random.default_rng(seed)
    allocations = rng.multivariate_hypergeometric(sizes, depth)
    outdir = ROOT/'peaks'/'matched_depth'/sample
    outdir.mkdir(parents=True, exist_ok=True)
    outfile = outdir/f'{sample}_peaks.narrowPeak'
    if outfile.exists():
        raise FileExistsError(outfile)
    with tempfile.TemporaryDirectory(prefix=f'atac-depth-{sample}-') as tmp:
        bed = Path(tmp)/'fragments.bed'
        with open(bed,'w',buffering=1024*1024) as handle:
            for chrom, n in zip(chroms, allocations):
                start,end = fragments[chrom+'_start'],fragments[chrom+'_end']
                indices = rng.choice(len(start), size=int(n), replace=False)
                for offset in range(0,len(indices),100000):
                    selected = indices[offset:offset+100000]
                    handle.writelines(f'{chrom}\t{a}\t{b}\n' for a,b in zip(start[selected],end[selected]))
        command = [str(Path(sys.executable).parent/'macs3'),'callpeak','-t',str(bed),
                   '-f','BEDPE','-g',str(config['effective_genome_size']),'-n',sample,
                   '--outdir',str(outdir),'-q',str(config['peak_qvalue']),
                   '--keep-dup','all','--call-summits']
        with open(ROOT/f'logs/{sample}_matched_depth_macs3.log','w') as log:
            log.write('Exact fragment depth: '+str(depth)+'; seed: '+str(seed)+'\n')
            log.flush()
            subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)
    with open(outfile) as handle:
        peaks = sum(bool(line.strip()) and not line.startswith('#') for line in handle)
    print(sample, 'matched-depth peaks:',peaks,flush=True)
    return {'sample':sample,'matched_fragments':depth,'seed':seed,'matched_depth_peaks':peaks}


if __name__ == '__main__':
    samples = pd.read_csv(ROOT/'config/samples.tsv',sep='\t')['sample'].tolist()
    depth = min(json.loads((ROOT/f'qc/{s}_preparation.json').read_text())['accepted_fragments'] for s in samples)
    with ThreadPoolExecutor(max_workers=2) as pool:
        jobs = [pool.submit(call_matched,s,depth,20260908+i) for i,s in enumerate(samples)]
        rows = [job.result() for job in jobs]
    pd.DataFrame(rows).to_csv(ROOT/'qc/matched_depth_peak_counts.tsv',sep='\t',index=False)
