"""Selected A549 S01: QC, saved model comparisons and Y-gene sample expression."""
from pathlib import Path
import hashlib
import json
import shutil
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / 'figures/chapter3'
DATA = OUT / 'source_data/a549_qc_sensitivity_combined'
sys.path.insert(0, str(ROOT / 'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style

BLUE, RED, TEAL, BROWN = '#3C5488', '#E64B35', '#00A087', '#7E6148'
STEM = 'ch03_supp_a549_rna_qc'
GENES = ['KDM5D', 'CALB2', 'PDE3A', 'FGA']


def render():
    font = apply_style(8)
    DATA.mkdir(parents=True, exist_ok=True)
    inputs = {
        'primary_PCA_QC.csv': OUT / 'source_data/a549_qc_expanded/primary_PCA_QC.csv',
        'PCA_associations.csv': OUT / 'source_data/a549_qc_expanded/PCA_associations.csv',
        'DE_nuclear_unadjusted.csv': OUT / 'source_data/a549_primary_nuclear_RINe/DE_nuclear_unadjusted.csv',
        'DE_nuclear_RINe_adjusted.csv': OUT / 'source_data/a549_primary_nuclear_RINe/DE_nuclear_RINe_adjusted.csv',
        'Y_gene_sample_logCPM.csv': OUT / 'source_data/a549_Y_inclusive_volcano/Y_gene_sample_logCPM.csv',
        'Y_gene_statistics.csv': OUT / 'source_data/a549_Y_inclusive_volcano/Y_gene_statistics.csv',
    }
    for name, path in inputs.items():
        shutil.copy2(path, DATA / name)
    qc = pd.read_csv(inputs['primary_PCA_QC.csv'])
    associations = pd.read_csv(inputs['PCA_associations.csv']).set_index('metric')
    un = pd.read_csv(inputs['DE_nuclear_unadjusted.csv'])
    de = pd.read_csv(inputs['DE_nuclear_RINe_adjusted.csv'])
    assert len(qc) == 8 and qc.sample_key.is_unique
    assert qc.groupby('treatment').size().to_dict() == {'CTRL': 4, 'LOY': 4}
    assert len(de) == len(un) == 22843 and set(de.gene_id) == set(un.gene_id)
    assert (de.FDR < .05).sum() == 27 and (un.FDR < .05).sum() == 16
    fig = plt.figure(figsize=(7.2, 5.8))
    positions = [.085, .395, .705]
    top_axes = [fig.add_axes([x, .615, .265, .29]) for x in positions]
    variables = [('integrity', 'RINe', 1),
                 ('qm_exonic_pct', 'Exonic reads\n(% of classified)', 1),
                 ('mitochondrial_fraction', 'Mitochondrial reads\n(% of assigned)', 100)]
    qc_display = []
    for i, (ax, (key, label, scale)) in enumerate(zip(top_axes, variables)):
        rho = float(associations.loc[key, 'spearman_rho'])
        assert int(associations.loc[key, 'n']) == 8
        for row in qc.itertuples():
            x = getattr(row, key) * scale
            colour = BLUE if row.treatment == 'CTRL' else RED
            ax.scatter(x, row.PC1, c=colour, s=25)
            ax.annotate(str(int(row.sample_index)), (x, row.PC1), xytext=(4, 3),
                        textcoords='offset points', fontsize=6.5, color=colour)
            qc_display.append({'panel': 'abc'[i], 'sample_key': row.sample_key,
                               'metric': key, 'x': x, 'PC1': row.PC1, 'rho': rho})
        ax.set(xlabel=label, ylabel='PC1' if i == 0 else '')
        ax.margins(.17, .22)
        ax.tick_params(labelsize=7)
        if i:
            ax.tick_params(labelleft=False)
        ax.text(.97, .025, f'ρ = {rho:.2f}; n = 8', transform=ax.transAxes,
                ha='right', fontsize=7)
        fig.text(positions[i] - .055, .927, 'abc'[i], fontsize=11, weight='semibold')
    fig.legend(handles=[Line2D([], [], marker='o', ls='', color=c, label=l)
                        for c, l in [(BLUE, 'sgCTRL'), (RED, 'Y-depleted')]],
               loc='upper center', bbox_to_anchor=(.53, .993), ncol=2, frameon=False,
               fontsize=7.5, handletextpad=.4, columnspacing=1.2)

    xlim = (min(-10, np.floor(min(de.logFC.min(), un.logFC.min())) - 1),
            max(5, np.ceil(max(de.logFC.max(), un.logFC.max())) + 1))
    ylim = (-.15, max(-np.log10(de.PValue.min()), -np.log10(un.PValue.min())) + 1.1)
    bottom_axes = [fig.add_axes([.085, .085, .385, .355]),
                   fig.add_axes([.585, .085, .385, .355])]
    positions_by_model = [
        {'KDM5D': (-10.7, 5.5), 'CALB2': (-4.2, 2.0), 'PDE3A': (2.6, 3.0), 'FGA': (3.1, 1.1)},
        {'KDM5D': (-10.1, 5.7), 'CALB2': (-4.0, 6.0), 'PDE3A': (2.5, 5.85), 'FGA': (3.2, 4.8)},
    ]
    saved_labels = []
    cutoffs = {}
    for i, (ax, values, model) in enumerate(zip(bottom_axes, [un, de], ['Condition only', 'Condition + RINe'])):
        y_gene = values.chrom.eq('chrY') & ~values.PAR
        other = ~y_gene & ~values.PAR
        significant = values.FDR < .05
        for mask, colour, size in [(other & ~significant, '#C8C8C8', 3),
                                   (other & significant & (values.logFC < 0), BLUE, 13),
                                   (other & significant & (values.logFC > 0), RED, 13),
                                   (values.PAR, TEAL, 17), (y_gene, BROWN, 18)]:
            points = values[mask]
            ax.scatter(points.logFC, -np.log10(points.PValue), s=size, c=colour,
                       alpha=.6 if size == 3 else 1, linewidths=.25,
                       edgecolors='white' if size > 3 else 'none', rasterized=size == 3,
                       zorder=2 if size == 3 else 4)
        ax.set(xlim=xlim, ylim=ylim, xlabel='log₂FC (Y-depleted / sgCTRL)',
               ylabel='−log₁₀(nominal P)')
        for x in [-1, 1]:
            ax.axvline(x, c='#AAAAAA', ls=':', lw=.5)
        cutoff = float(values.loc[significant, 'PValue'].max())
        cutoffs[model] = cutoff
        ax.axhline(-np.log10(cutoff), c='#666666', ls='--', lw=.6)
        for gene in GENES:
            row = values[values.gene_name.eq(gene)].iloc[0]
            ax.annotate(gene, (row.logFC, -np.log10(row.PValue)),
                        xytext=positions_by_model[i][gene], textcoords='data', fontsize=6.5,
                        arrowprops=dict(arrowstyle='-', color='#888888', lw=.45), zorder=6)
            saved_labels.append({'model': model, 'gene_id': row.gene_id, 'gene_name': gene,
                                 'logFC': row.logFC, 'PValue': row.PValue, 'FDR': row.FDR})
        ax.text(.025, .975, model, transform=ax.transAxes, fontsize=7.3, va='top')
        ax.text(.98, .90, f'{int(significant.sum())} significant\n{int((significant & other).sum())} non-Y/non-PAR',
                transform=ax.transAxes, ha='right', va='top', fontsize=7)
        fig.text([.03, .53][i], .454, 'de'[i], fontsize=11, weight='semibold')
    fig.legend(handles=[Line2D([], [], marker='o', ls='', color=c, label=l)
                        for c, l in [(BROWN, 'Y non-PAR'), (TEAL, 'PAR'),
                                     (BLUE, 'Other: lower'), (RED, 'Other: higher')]],
               loc='center', bbox_to_anchor=(.53, .487), ncol=4, fontsize=7,
               frameon=False, handletextpad=.4, columnspacing=1.3)

    # Extend below a-e, retaining their physical sizes and relative spacing.
    original_height, extra_height = 5.8, 2.8
    total_height = original_height + extra_height
    axes_positions = [(ax, ax.get_position().bounds) for ax in fig.axes]
    text_positions = [(text, text.get_position()) for text in fig.texts]
    legend_positions = [(legend, legend.get_bbox_to_anchor().transformed(fig.transFigure.inverted()).bounds)
                        for legend in fig.legends]
    def translated_y(y):
        return (y * original_height + extra_height) / total_height
    fig.set_size_inches(7.2, total_height)
    for ax, (x, y, width, height) in axes_positions:
        ax.set_position([x, translated_y(y), width, height * original_height / total_height])
    for text, (x, y) in text_positions:
        text.set_position((x, translated_y(y)))
    for legend, (x, y, width, height) in legend_positions:
        legend.set_bbox_to_anchor((x, translated_y(y), width, height * original_height / total_height),
                                  transform=fig.transFigure)

    # Frozen S04b matrix: same genes, order, raw logCPM scale and statistical support.
    y_stats = pd.read_csv(inputs['Y_gene_statistics.csv'])
    sample_matrix = pd.read_csv(inputs['Y_gene_sample_logCPM.csv'], index_col=0)
    samples = [f'A549-{group}-{i}' for group in ['CTRL', 'LOY'] for i in range(1, 5)]
    assert sample_matrix.shape == (13, 8)
    assert list(sample_matrix.index) == y_stats.gene_id.tolist()
    assert list(sample_matrix.columns) == samples
    assert np.isfinite(sample_matrix.to_numpy()).all()
    pd.testing.assert_frame_equal(y_stats[['gene_id', 'gene_name', 'PValue', 'FDR']].reset_index(drop=True),
                                  de.set_index('gene_id').loc[y_stats.gene_id].reset_index()[['gene_id', 'gene_name', 'PValue', 'FDR']],
                                  check_exact=False, rtol=1e-12, atol=0)
    heat = fig.add_axes([.235, .085, .44, .218])
    vmin = float(np.floor(sample_matrix.to_numpy().min()))
    vmax = float(np.ceil(sample_matrix.to_numpy().max()))
    im = heat.imshow(sample_matrix.to_numpy(), cmap='magma', aspect='auto',
                     vmin=vmin, vmax=vmax, interpolation='nearest')
    heat.set_yticks(range(len(y_stats)), y_stats.gene_name, fontsize=7)
    for tick, gene in zip(heat.get_yticklabels(), y_stats.gene_name):
        if gene == 'KDM5D':
            tick.set_fontweight('bold')
    heat.set_xticks(range(8), [str(i) for i in range(1, 5)] * 2, fontsize=7)
    heat.tick_params(length=0)
    heat.axvline(3.5, color='white', lw=2)
    for start, colour, label in [(-.5, BLUE, 'sgCTRL'), (3.5, RED, 'Y-depleted')]:
        heat.plot([start, start + 4], [-.85, -.85], color=colour, lw=3, clip_on=False)
        heat.text(start + 2, -1.25, label, ha='center', va='bottom', fontsize=7.5, color=colour)
    # Group strips must not expand the heatmap data limits relative to the statistics.
    heat.set_ylim(len(y_stats) - .5, -.5)
    heat.spines[:].set_visible(False)
    table = fig.add_axes([.705, .085, .265, .218])
    table.set(xlim=(0, 2), ylim=(len(y_stats) - .5, -.5))
    table.axis('off')
    for x, label in [(.45, 'Nominal P'), (1.5, 'FDR')]:
        table.text(x, -1.25, label, ha='center', va='bottom', fontsize=7.5)
    for i, row in enumerate(y_stats.itertuples()):
        for x, value in [(.45, row.PValue), (1.5, row.FDR)]:
            label = f'{value:.2e}' if value < .001 else f'{value:.3f}'
            table.text(x, i, label, ha='center', va='center', fontsize=7,
                       fontweight='bold' if row.gene_name == 'KDM5D' else 'normal')
    fig.text(.03, .320, 'f', fontsize=11, weight='semibold')
    colour_axis = fig.add_axes([.235, .038, .44, .010])
    colourbar = fig.colorbar(im, cax=colour_axis, orientation='horizontal')
    colourbar.set_label('TMM-normalised log₂ CPM', fontsize=7, labelpad=1)
    colourbar.ax.tick_params(labelsize=7, pad=1, length=2)
    pd.DataFrame(qc_display).to_csv(DATA / 'panels_abc_display.csv', index=False)
    pd.DataFrame(saved_labels).to_csv(DATA / 'panels_de_label_values.csv', index=False)
    exports = {}
    for ext in ['pdf', 'svg', 'png']:
        path = OUT / f'{STEM}.{ext}'
        fig.savefig(path, dpi=300, facecolor='white')
        exports[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
    plt.close(fig)
    provenance = {
        'scope': 'S01 b-d become a-c; S02 a-b become d-e; S04b becomes f. Original S01 a/e and S04a omitted. No refitting.',
        'font': font, 'size_inches': [7.2, total_height], 'samples': 8, 'gene_family': 22843,
        'heatmap': {'genes': y_stats.gene_name.tolist(), 'samples': samples,
                    'colour_limits': [vmin, vmax], 'row_scaled': False, 'RINe_adjusted_display': False},
        'inputs': {name: {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
                   for name, path in inputs.items()},
        'volcano_limits': {'x': xlim, 'y': ylim}, 'saved_BH_nominal_boundaries': cutoffs,
        'exports': exports,
    }
    (DATA / 'provenance.json').write_text(json.dumps(provenance, indent=2), encoding='utf-8')
    print('Wrote six-panel S01: QC, two model comparisons and the 13-gene Y-expression heatmap.')


if __name__ == '__main__':
    render()
