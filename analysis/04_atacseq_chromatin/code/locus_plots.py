"""Display all replicate coverage tracks at prespecified top differential loci."""
from contextlib import ExitStack
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyBigWig


def create_loci(root, annotated, samples):
    selected = []
    for direction in ['gained', 'lost']:
        count = 0
        for row in annotated[annotated.direction == direction].sort_values(
                ['FDR', 'peak_id']).itertuples(index=False):
            if any(row.chrom == old.chrom and abs(row.start-old.start) < 10000
                   for old in selected):
                continue
            selected.append(row)
            count += 1
            if count == 2:
                break
    if not selected:
        return
    annotated.set_index('peak_id').loc[[row.peak_id for row in selected]].reset_index().to_csv(
        root/'annotation/example_loci.csv', index=False)
    fig, axes = plt.subplots(len(selected), 1, figsize=(12, 3.4*len(selected)),
                             squeeze=False)
    with ExitStack() as stack:
        tracks = {sample: stack.enter_context(pyBigWig.open(str(
            root/f'tracks/{sample}.fragment_cpm.bw'))) for sample in samples['sample']}
        sizes = next(iter(tracks.values())).chroms()
        for ax, row in zip(axes[:, 0], selected):
            left = max(0, ((row.start+250-5000)//25)*25)
            right = min((sizes[row.chrom]//25)*25, left+10000)
            coordinates = np.arange(left, right, 25)+12.5
            for i, sample in samples.iterrows():
                values = tracks[sample['sample']].values(row.chrom, left, right,
                                                        numpy=True)
                coverage = np.nan_to_num(values, nan=0).reshape(-1, 25).mean(axis=1)
                ax.plot(coordinates, coverage, label=sample['sample'],
                        color='#0072B2' if sample.condition == 'CTRL' else '#D55E00',
                        linestyle=['-', '--', ':'][i % 3], linewidth=1, alpha=.85)
            ax.axvspan(row.start, row.end, color='#777777', alpha=.15)
            symbol = row.nearest_symbol if pd.notna(row.nearest_symbol) else 'unassigned'
            ax.set_title(f'{row.peak_id} | nearest gene: {symbol} | '
                         f'LOY/CTRL log2FC={row.logFC:.2f}, FDR={row.FDR:.2g}', fontsize=10)
            ax.set_xlabel(f'{row.chrom} position (hg38, 0-based)')
            ax.set_ylabel('Fragment coverage\nper million non-Y fragments')
            ax.ticklabel_format(axis='x', style='plain', useOffset=False)
            ax.spines[['top', 'right']].set_visible(False)
            ax.legend(ncol=6, fontsize=8, loc='upper right')
    fig.tight_layout()
    fig.savefig(root/'annotation/example_loci.png', dpi=180)
    plt.close(fig)
