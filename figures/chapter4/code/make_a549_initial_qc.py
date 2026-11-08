"""Initial A549 screen QC review from frozen count-derived tables."""
from pathlib import Path
import hashlib
import json
import os
import sys

THESIS = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(THESIS.parent / '.figure_python_deps'))
OUT = THESIS / 'figures/chapter4'
DATA = OUT / 'source_data/a549_initial_qc'
os.environ.setdefault('MPLCONFIGDIR', str(DATA / 'mpl_cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pymupdf
from matplotlib.lines import Line2D

sys.path.insert(0, str(THESIS / 'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style, condition_colours
apply_style(7.2)
COL = condition_colours()
COLOURS = {'CTRL': COL['Y-retaining'], 'LOY': COL['Y-depleted']}
STEM = 'ch04_supp_a549_screen_initial_qc_draft'


def main():
    tables = [p for p in DATA.glob('*.csv')]
    hashes = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in tables}
    qc = pd.read_csv(DATA / 'sample_qc.csv')
    corr = pd.read_csv(DATA / 'spearman_all_samples.csv', index_col=0)
    pca = pd.read_csv(DATA / 'pca_all_samples.csv')
    genes = pd.read_csv(DATA / 'baseline_gene_qc.csv')
    branches = pd.read_csv(DATA / 'saved_drugz_branch_comparison.csv', index_col=0).dropna()
    summary = json.loads((DATA / 'summary.json').read_text(encoding='utf-8'))
    labels = [s.replace('_', ' ').replace(' R', ' ') for s in qc['sample']]
    colors = qc.genotype.map(COLOURS).tolist()
    fig, axes = plt.subplots(3, 2, figsize=(183 / 25.4, 218 / 25.4))
    fig.subplots_adjust(left=.10, right=.985, top=.95, bottom=.065, hspace=.64, wspace=.40)
    titles = ['Guide mapping', 'Read concentration', 'Sample concordance',
              'Count-profile PCA', 'D0-to-P4 depletion', 'Saved drugZ branches']
    for letter, ax, title in zip('abcdef', axes.flat, titles):
        ax.text(-.17, 1.06, letter, transform=ax.transAxes, fontweight='semibold', fontsize=10)
        ax.set_title(title, fontsize=8, fontweight='semibold', pad=7)
        ax.spines[['top', 'right']].set_visible(False)
    x = np.arange(len(qc))
    for ax, values, ylabel in [
        (axes[0, 0], qc.mapping_percent, 'Mapped reads (%)'),
        (axes[0, 1], qc.top100_read_fraction * 100, 'Reads in top 100 guides (%)'),
    ]:
        ax.bar(x, values, color=colors, width=.72)
        ax.set_xticks(x, labels, rotation=70, ha='right', fontsize=6.2)
        ax.set(ylabel=ylabel, ylim=(0, 100 if ax is axes[0, 0] else 1.2))
        for split in (3.5,):
            ax.axvline(split, color='#DDDDDD', lw=.6)
        ax.tick_params(axis='x', length=0)
    ax = axes[1, 0]
    im = ax.imshow(corr.to_numpy(), vmin=.7, vmax=1, cmap='viridis', interpolation='nearest')
    ax.set_xticks(x, labels, rotation=70, ha='right', fontsize=6.5)
    ax.set_yticks(x, labels, fontsize=6.5)
    ax.tick_params(length=0)
    for split in (3.5,):
        ax.axhline(split, color='white', lw=.65)
        ax.axvline(split, color='white', lw=.65)
    for tick, sample in zip(ax.get_yticklabels(), qc['sample']):
        tick.set_color(COLOURS[sample.split('_')[0]])
    cb = fig.colorbar(im, ax=ax, fraction=.045, pad=.03, ticks=[.7, .8, .9, 1])
    cb.set_label('Spearman ρ', fontsize=6.5)
    cb.ax.tick_params(labelsize=6)
    ax = axes[1, 1]
    for _, row in pca.iterrows():
        ax.scatter(row.PC1, row.PC2, s=28, marker={'D0': 's', 'P4': 'o'}[row.stage],
                   c=COLOURS[row.genotype], edgecolors='white', linewidths=.5, zorder=3)
    annotations = [
        ('CTRL', 'D0', 'CTRL D0 1–2', (0, -15)),
        ('LOY', 'D0', 'LOY D0 1–2', (0, -15)),
        ('CTRL', 'P4', 'CTRL P4 2', (3, 8)),
        ('LOY', 'P4', 'LOY P4 1–2', (0, 10)),
    ]
    for genotype, stage, label, offset in annotations:
        values = pca.loc[pca.genotype.eq(genotype) & pca.stage.eq(stage)]
        ax.annotate(label, (values.PC1.mean(), values.PC2.mean()), xytext=offset,
                    textcoords='offset points', fontsize=6,
                    ha='left', color=COLOURS[genotype])
    variance = summary['pca_explained_variance']
    ax.set(xlabel=f'PC1 ({variance[0]:.1%})', ylabel=f'PC2 ({variance[1]:.1%})')
    ax.set_xlim(pca.PC1.min() - .15 * np.ptp(pca.PC1), pca.PC1.max() + .45 * np.ptp(pca.PC1))
    ax.set_ylim(pca.PC2.min() - .25 * np.ptp(pca.PC2), pca.PC2.max() + .25 * np.ptp(pca.PC2))
    ax = axes[2, 0]
    groups = [('CTRL', 'Essential'), ('CTRL', 'Non-essential'), ('LOY', 'Essential'), ('LOY', 'Non-essential')]
    values = [genes.loc[genes.genotype.eq(g) & genes.category.eq(c), 'lfc'] for g, c in groups]
    boxes = ax.boxplot(values, positions=[1, 2, 3.6, 4.6], widths=.55, patch_artist=True, showfliers=True,
        medianprops={'color': '#222222', 'linewidth': .9},
        boxprops={'linewidth': .65}, whiskerprops={'linewidth': .65}, capprops={'linewidth': .65},
        flierprops={'marker': 'o', 'markerfacecolor': '#777777', 'markersize': 1.2,
                    'markeredgewidth': 0, 'alpha': .35})
    for patch, (genotype, category) in zip(boxes['boxes'], groups):
        patch.set_facecolor(COLOURS[genotype]); patch.set_alpha(.75 if category == 'Essential' else .25)
    ax.axhline(0, lw=.6, ls='--', color='#999999', zorder=0)
    ax.set_xticks([1, 2, 3.6, 4.6], ['CTRL\nE', 'CTRL\nNE', 'LOY\nE', 'LOY\nNE'])
    ax.set_ylabel('Gene log₂(P4 / D0)')
    quality = summary['biological_qc']
    ax.text(.02, .03, 'Essential/non-essential ROC AUC\n' +
        ' · '.join(f"{v['genotype']} {v['roc_auc']:.3f}" for v in quality),
        transform=ax.transAxes, fontsize=6, va='bottom')
    ax = axes[2, 1]
    ax.scatter(branches.normZ_original, branches.normZ_rerun, s=1.6,
               c='#777E86', alpha=.35, edgecolors='none', rasterized=True)
    limits = (-5.9, 3.9)
    ax.plot(limits, limits, ls='--', color='#AAB0B7', lw=.7, zorder=0)
    for gene, offset in [('EIF1AX', (13, -10)), ('DDX3X', (-15, 23))]:
        row = branches.loc[gene]
        ax.scatter(row.normZ_original, row.normZ_rerun, c=COLOURS['LOY'], s=20, zorder=4)
        ax.annotate(gene, (row.normZ_original, row.normZ_rerun), xytext=offset,
            textcoords='offset points', fontsize=6.5,
            arrowprops={'arrowstyle': '-', 'lw': .5, 'color': '#777777'})
    ax.text(.03, .97, f"ρ = {summary['saved_drugz_normz_spearman']:.5f}\n19,112 matched genes",
            transform=ax.transAxes, va='top', fontsize=6.3)
    ax.set(xlim=limits, ylim=limits, xlabel='Original drugZ normZ', ylabel='Rerun drugZ normZ')
    for ext in ('pdf', 'svg', 'png'):
        fig.savefig(OUT / f'{STEM}.{ext}', dpi=600)
    preview = OUT / 'review' / f'{STEM}_preview.png'
    fig.savefig(preview, dpi=190)
    plt.close(fig)
    pdf = pymupdf.open(OUT / f'{STEM}.pdf')
    page = pdf[0]
    assert 'P3' not in page.get_text()
    outside = [w for w in page.get_text('words') if
        w[0] < -.5 or w[1] < -.5 or w[2] > page.rect.width + .5 or w[3] > page.rect.height + .5]
    assert not outside, outside
    assert all('OpenSans' in font[3].replace(' ', '') for font in page.get_fonts(full=True))
    for path, digest in hashes.items():
        assert hashlib.sha256(Path(path).read_bytes()).hexdigest() == digest
    validation = {'figure_mm': [183, 218], 'png_dpi': 600, 'pdf_pages': len(pdf),
        'pdf_text_outside_page': 0, 'all_fonts_open_sans': True,
        'svg_editable_text': '<text' in (OUT / f'{STEM}.svg').read_text(encoding='utf-8'),
        'source_tables_unchanged': True, 'preview': str(preview),
        'retained_libraries': len(qc), 'p3_present_in_figure': False}
    (DATA / 'figure_validation.json').write_text(json.dumps(validation, indent=2), encoding='utf-8')
    print(json.dumps(validation, indent=2))


if __name__ == '__main__':
    main()
