"""Native TOBIAS aggregate profiles for three TFAP2C models on shared sites."""
from pathlib import Path
import json
import subprocess
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'motifs/TF_extension/TFAP2C_review'
FP = ROOT/'motifs/TF_extension/TOBIAS'
ids = ['MA0524.3', 'MA0814.3', 'MA0815.1']
signals, labels = [], []
for kind in ['uncorrected', 'corrected']:
    for condition in ['CTRL', 'LOY']:
        source = FP/f'{condition}/{condition}_{kind}.bw'
        cache = Path('/tmp/atacseq-sept08-2026-tf-tools/aggregate')/source.name
        chosen = cache if cache.exists() and cache.stat().st_size == source.stat().st_size else source
        signals.append(str(chosen)); labels.append(f'{condition}_{kind}')
fig, axes = plt.subplots(3, 2, figsize=(11, 10), sharex=True)
all_curves = []
for row, motif in enumerate(ids):
    bed = FP/f'BINDetect/{motif}/beds/{motif}_all.bed'
    table = OUT/f'{motif}_aggregate_native.txt'
    command = ['/tmp/atacseq-sept08-2026-tf-tools/tobias_venv/bin/TOBIAS', 'PlotAggregate',
               '--TFBS', str(bed), '--signals', *signals, '--TFBS-labels', motif,
               '--signal-labels', *labels, '--output', str(OUT/f'{motif}_aggregate_native.pdf'),
               '--output-txt', str(table), '--flank', '60', '--remove-outliers', '1',
               '--smooth', '1', '--share-y', 'signals', '--plot-boundaries',
               '--title', f'TFAP2C {motif}: all shared motif sites']
    (OUT/f'{motif}_aggregate_command.json').write_text(json.dumps(command, indent=2)+'\n')
    with (OUT/f'{motif}_aggregate.log').open('w') as log:
        subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
    curves = {}
    for line in table.read_text().splitlines():
        if not line.strip() or line.startswith('#'):
            continue
        label, regions, values = line.split('\t')
        curves[label] = np.array([float(v) for v in values.split(',')])
    assert set(curves) == set(labels)
    assert all(len(v) == 120 and np.isfinite(v).all() for v in curves.values())
    x = np.r_[np.arange(-60, 0), np.arange(1, 61)]
    all_curves.append(pd.DataFrame({'motif_id': motif, 'bp_from_center': x, **curves}))
    for col, kind in enumerate(['uncorrected', 'corrected']):
        ax = axes[row, col]
        for condition, color in [('CTRL', '#0072B2'), ('LOY', '#D55E00')]:
            ax.plot(x, curves[f'{condition}_{kind}'], color=color, label=condition, lw=1.2)
        ax.axvline(0, color='grey', ls='--', lw=.6)
        ax.set_title(f'{motif}: {kind}')
        ax.set_ylabel('Mean Tn5 signal')
        ax.legend(frameon=False)
    lo = min(axes[row, 0].get_ylim()[0], axes[row, 1].get_ylim()[0])
    hi = max(axes[row, 0].get_ylim()[1], axes[row, 1].get_ylim()[1])
    for ax in axes[row]: ax.set_ylim(lo, hi)
    print(f'Completed {motif}', flush=True)
for ax in axes[-1]: ax.set_xlabel('bp from motif center (native TOBIAS coordinates)')
fig.suptitle('TFAP2C: aggregate profiles over identical sites in CTRL and LOY')
fig.text(.5, .01, 'All motif sites; no condition-specific bound-site selection, smoothing, outlier removal or rescaling.\n'
         'Condition pools provide descriptive profiles, not biological-replicate confidence intervals.', ha='center', fontsize=9)
fig.tight_layout(rect=(0, .055, 1, .96))
for ext in ['png', 'pdf']: fig.savefig(OUT/f'TFAP2C_aggregate_footprints.{ext}', dpi=180)
pd.concat(all_curves, ignore_index=True).to_csv(OUT/'TFAP2C_aggregate_curves.csv', index=False)
