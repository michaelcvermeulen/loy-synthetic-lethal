"""Compact PAR barplot for F01; frozen CPM, matching the A549/DU145 display."""
import hashlib
import json
import numpy as np
import pandas as pd
from matplotlib.patches import Patch


def add_par_panel(fig, figure_root, blue, red):
    source = figure_root / 'source_data/shared_sex_chromosome/HT1080'
    inputs = [source / 'gene_period_values.csv', source / 'region_sample_sums.csv']
    hashes = {str(p.relative_to(figure_root)): hashlib.sha256(p.read_bytes()).hexdigest()
              for p in inputs}
    genes = pd.read_csv(inputs[0])
    genes = genes[genes.region.isin(['PAR1', 'PAR2']) &
                  ~genes.gene_name.str.startswith('ENSG')].copy()
    assert not genes.duplicated(['gene_id', 'batch']).any()
    # Author-requested display exclusion; all source features remain in regional sums.
    names = genes.loc[genes.CTRL.ge(5) & ~genes.gene_name.eq('CSF2RA')].sort_values(['region', 'start']).gene_name.drop_duplicates().tolist()
    plotted = genes[genes.gene_name.isin(names)].copy()
    plotted['displayed'] = plotted.CTRL.ge(5)
    plotted['kind'] = 'gene'
    sums = pd.read_csv(inputs[1])
    sums = sums[sums.region.isin(['PAR1', 'PAR2'])]
    assert sums.groupby(['batch', 'region', 'treatment']).size().eq(3).all()
    summary = sums.groupby(['batch', 'region', 'treatment']).CPM.mean().unstack('treatment').reset_index()
    summary['percent'] = 100 * summary.LOY / summary.CTRL
    pooled = sums.groupby(['batch', 'sample', 'treatment']).CPM.sum().reset_index()
    combined = pooled.groupby(['batch', 'treatment']).CPM.mean().unstack('treatment').reset_index()
    combined['percent'] = 100 * combined.LOY / combined.CTRL
    combined['gene_name'] = 'Combined PAR'
    combined['region'] = 'Combined PAR'
    combined['kind'] = 'combined sum'
    combined['displayed'] = True
    data = pd.concat([plotted, combined], ignore_index=True)
    combined_x = len(names) + .45
    data['x'] = data.gene_name.map(dict(zip(names + ['Combined PAR'],
                                          list(np.arange(len(names))) + [combined_x])))
    colours = {'PAR1': '#4DBBD5', 'PAR2': '#8491B4', 'Combined PAR': '#555555'}
    ax = fig.add_axes([.725, .35035, .25, .2769])
    ax.axhline(100, color='#999999', ls='--', lw=.6, zorder=1)
    ax.axvline(len(names) - .25, color='#CCCCCC', lw=.6)
    for batch, offset, hatch in [('AUG_2025', -.19, ''), ('MAY_2026', .19, '///')]:
        q = data[data.batch.eq(batch) & data.displayed]
        ax.bar(q.x + offset, q.percent, width=.36,
               color=[colours[r] for r in q.region], hatch=hatch,
               edgecolor='white', lw=.25, zorder=2)
    for r in data[~data.displayed].itertuples():
        ax.text(r.x + .19, -.025, '–', transform=ax.get_xaxis_transform(),
                ha='center', va='top', color='#777777', fontsize=7)
    ax.set_xticks(list(np.arange(len(names))) + [combined_x], names + ['Combined PAR'],
                  rotation=90, fontsize=5.6)
    for label in ax.get_xticklabels()[:-1]:
        label.set_fontstyle('italic')
    ax.set(xlim=(-.65, combined_x + .8), ylim=(0, 180), ylabel='% of sgCTRL')
    ax.set_yticks([0, 50, 100, 150])
    ax.tick_params(axis='y', labelsize=7)
    ax.yaxis.label.set_size(8)
    fig.text(.69, .645, 'f', fontsize=12, weight='semibold')
    fig.text(.725, .645, 'Pseudoautosomal expression', fontsize=7.5, weight='semibold')
    handles = [Patch(facecolor=colours[r], edgecolor='none', label=r) for r in ['PAR1', 'PAR2']]
    handles += [Patch(facecolor='white', edgecolor='#777777', linewidth=.4, hatch=h, label=l)
                for h, l in [('', 'Aug'), ('///', 'May')]]
    ax.legend(handles=handles, ncol=4, fontsize=5.5, frameon=False,
              loc='upper center', bbox_to_anchor=(.5, 1.01),
              handlelength=.9, handletextpad=.3, columnspacing=.7)
    for i, region in enumerate(['PAR1', 'PAR2']):
        values = summary[summary.region.eq(region)].set_index('batch').percent
        ax.text(.98, .86 - i * .058,
                f'{region}: {values.AUG_2025:.1f}% / {values.MAY_2026:.1f}%',
                transform=ax.transAxes, ha='right', fontsize=5.5, color='#555555')
    values = combined.set_index('batch').percent
    ax.text(.98, .744, f'Combined: {values.AUG_2025:.1f}% / {values.MAY_2026:.1f}%',
            transform=ax.transAxes, ha='right', fontsize=5.5, color='#555555')
    output = figure_root / 'source_data/ht1080_rna_expression_changes'
    (output / 'panel_f_PAR_provenance.json').write_text(json.dumps({
        'date': '2026-10-09', 'inputs_sha256': hashes,
        'comparison': 'Ratio of mean TMM-normalised CPM, LOY / CTRL, within sampling period',
        'regional_rule': 'All saved annotated PAR features included in regional and combined sums',
        'gene_rule': 'Named chrX PAR genes with mean sgCTRL CPM >=5 in the displayed period, excluding CSF2RA at author request',
        'omission': 'CSF2RA omitted from individual-gene bars in both periods; retained in regional and combined sums',
        'layout': 'Third column of middle row, beside d and e; former bottom PAR row removed',
        'style': 'PAR1 cyan, PAR2 blue-grey and combined PAR grey; solid August, hatched May',
        'sampling': 'Three passage samples per condition per period; one established population per condition',
        'statistics': 'Descriptive display; no inferential tests or primary DE refitting',
        'PAR_origin': 'Combined X- and Y-derived transcripts; not Y-specific',
        'figure_size_inches': list(fig.get_size_inches()),
    }, indent=2) + chr(10), encoding='utf-8')
    assert all(hashlib.sha256(p.read_bytes()).hexdigest() == hashes[str(p.relative_to(figure_root))]
               for p in inputs)
    return data, summary
