"""Read existing 25 bp fragment-CPM BigWigs; cache exact window averages for plotting."""
from pathlib import Path
from contextlib import ExitStack
import hashlib
import json
import math
import shutil
import numpy as np
import pandas as pd
import pyBigWig

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'integration/gene_atlas'
windows = pd.read_csv(OUT / 'data/plot_windows.csv')
samples = pd.read_csv(ROOT / 'config/samples.tsv', sep='\t')['sample'].tolist()
cache = {}
audit = []
manifest = []
with ExitStack() as stack:
    tracks = {}
    for sample in samples:
        path = ROOT / f'tracks/{sample}.fragment_cpm.bw'
        tracks[sample] = stack.enter_context(pyBigWig.open(str(path)))
        with path.open('rb') as handle:
            sha = hashlib.file_digest(handle, 'sha256').hexdigest()
        manifest.append({'sample': sample, 'path': str(path.relative_to(ROOT)), 'sha256': sha})
    for i, row in enumerate(windows.itertuples(index=False)):
        path = OUT / f'coverage/{row.window_id}.npz'
        key = (row.chrom, row.start, row.end)
        width = max(25, int(math.ceil((row.end - row.start) / 1800 / 25)) * 25)
        starts = np.arange(row.start, row.end, width, dtype=np.int64)
        ends = np.minimum(starts + width, row.end)
        if key in cache:
            shutil.copyfile(cache[key], path)
        else:
            coverage = []
            for sample in samples:
                values = tracks[sample].values(row.chrom, row.start, row.end, numpy=True)
                assert len(values) == row.end - row.start
                values = np.nan_to_num(values, nan=0., posinf=0., neginf=0.)
                assert np.all(values >= 0)
                integral = np.add.reduceat(values.astype(np.float64), starts - row.start)
                means = integral / (ends - starts)
                assert np.isclose(np.sum(means * (ends - starts)), np.sum(values, dtype=np.float64), rtol=1e-10)
                coverage.append(means)
            np.savez_compressed(path, starts=starts, ends=ends, values=np.array(coverage), samples=np.array(samples), bin_width=width)
            cache[key] = path
        audit.append({'window_id': row.window_id, 'chrom': row.chrom, 'start': row.start, 'end': row.end,
                      'plot_bin_bp': width, 'bins': len(starts), 'samples': len(samples)})
        if (i + 1) % 25 == 0:
            print(f'Coverage {i + 1}/{len(windows)} windows', flush=True)
pd.DataFrame(audit).to_csv(OUT / 'data/coverage_audit.csv', index=False)
pd.DataFrame(manifest).to_csv(OUT / 'data/coverage_source_manifest.csv', index=False)
(OUT / 'coverage/COMPLETE.json').write_text(json.dumps({'windows': len(windows), 'unique_coordinate_windows': len(cache),
    'samples': samples, 'pyBigWig': pyBigWig.__version__, 'numpy': np.__version__,
    'source_normalization': 'Mean fragment-span coverage per million accepted non-Y fragments, original 25 bp BigWigs',
    'plot_binning': 'Exact mean including zero uncovered bases; at least 25 bp, up to approximately 1800 bins per window',
    'validation': 'All means preserve the full basewise signal integral'}, indent=2))
print(f'Saved coverage for {len(windows)} windows', flush=True)
