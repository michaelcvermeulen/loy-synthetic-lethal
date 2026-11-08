"""CTCF QC using shared motif sites, without selecting condition-specific bound sites."""
from pathlib import Path
import json, subprocess, shutil
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'motifs/TF_extension/TOBIAS'
assert (OUT/'ANALYSIS_COMPLETE.json').exists()
cv = pd.read_csv(ROOT/'motifs/TF_extension/chromVAR/motifs_independent_limma.csv')
motif_id = cv.loc[cv.motif_name.eq('CTCF'), 'id'].sort_values().iloc[0]
native = pd.read_csv(OUT/'BINDetect/bindetect_results.txt',sep='\t')
ctcf = native[native.motif_id.eq(motif_id)]
assert len(ctcf)==1
prefix = ctcf.output_prefix.iloc[0]
bed = OUT/f'BINDetect/{prefix}/beds/{prefix}_all.bed'
assert bed.exists(), bed
site_count = sum(1 for line in bed.open() if line.strip() and not line.startswith('#'))
assert site_count > 0
cache = Path('/tmp/atacseq-sept08-2026-tf-tools/aggregate')
cache.mkdir(exist_ok=True)
signals, labels = [], []
for kind in ['uncorrected', 'corrected']:
    for condition in ['CTRL', 'LOY']:
        src = OUT/f'{condition}/{condition}_{kind}.bw'
        dst = cache/src.name
        assert src.exists()
        if not dst.exists() or dst.stat().st_size != src.stat().st_size:
            shutil.copyfile(src, dst)
        signals.append(str(dst)); labels.append(f'{condition}_{kind}')
native_pdf = OUT/'CTCF_aggregate_native.pdf'
table = OUT/'CTCF_aggregate_native.txt'
command = ['/tmp/atacseq-sept08-2026-tf-tools/tobias_venv/bin/TOBIAS',
    'PlotAggregate', '--TFBS', str(bed), '--signals', *signals,
    '--TFBS-labels', f'CTCF_{motif_id}', '--signal-labels', *labels,
    '--output', str(native_pdf), '--output-txt', str(table),
    '--flank', '60', '--remove-outliers', '1', '--smooth', '1',
    '--share-y', 'signals', '--plot-boundaries', '--title', 'CTCF: shared motif-site QC']
(OUT/'29_PlotAggregate.command.json').write_text(json.dumps(command, indent=2)+'\n')
with (ROOT/'logs/29_PlotAggregate.log').open('w') as log:
    subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
curves = {}
for line in table.read_text().splitlines():
    if line.startswith('#') or not line.strip(): continue
    label, regions, values = line.split('\t')
    curves[label] = np.array([float(v) for v in values.split(',')])
assert set(curves) == set(labels)
assert all(len(x) == 120 and np.isfinite(x).all() for x in curves.values())
# Preserve the TOBIAS native coordinate convention (no zero position).
x = np.r_[np.arange(-60,0), np.arange(1,61)]
pd.DataFrame({'bp_from_center':x, **curves}).to_csv(OUT/'CTCF_aggregate.csv', index=False)
fig, axes = plt.subplots(1,2,figsize=(11,4.5),sharey=True)
for ax, kind in zip(axes, ['uncorrected','corrected']):
    for condition, color in [('CTRL','#0072b2'),('LOY','#d55e00')]:
        ax.plot(x, curves[f'{condition}_{kind}'], label=condition, color=color, lw=1.3)
    ax.axvline(0,color='grey',ls='--',lw=.7)
    ax.set(title=f'Tn5 signal: {kind}', xlabel='bp from motif center')
    ax.legend(frameon=False)
axes[0].set_ylabel('Mean signal (native ATACorrect output scale)')
fig.suptitle(f'CTCF {motif_id}: {site_count:,} shared motif sites')
fig.text(.5,.01,'All motif sites in tested peaks; no condition-specific binding filter, outlier removal, smoothing, or rescaling.',ha='center',fontsize=8)
fig.tight_layout(rect=(0,.04,1,.94))
for ext in ['png','pdf']: fig.savefig(OUT/f'CTCF_aggregate.{ext}',dpi=220)
(OUT/'29_PlotAggregate.complete.json').write_text(json.dumps(dict(motif_id=motif_id,output_prefix=prefix,site_count=site_count,shared_sites=True,command=command),indent=2)+'\n')
print(f'CTCF aggregate complete: {site_count} sites',flush=True)
