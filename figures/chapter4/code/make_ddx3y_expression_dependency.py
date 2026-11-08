"""Plot the saved expression association and descriptive partner correlations."""
from pathlib import Path
import json
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import pymupdf as fitz

THESIS = Path(__file__).resolve().parents[3]
OUT = THESIS / 'figures/chapter4'
DATA = OUT / 'source_data/ddx3y_expression_dependency'
STEM = 'ch04_supp_ddx3y_expression_dependency_draft'
sys.path.insert(0, str(THESIS / 'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style


def main():
    apply_style(8.5)
    manifest = json.loads((DATA / 'manifest.json').read_text())
    models = pd.read_csv(DATA / 'analysis_model_values.csv', index_col='ModelID')
    xy = models.loc[models.group.eq('XY')]
    statistics = pd.read_csv(DATA / 'ddx3x_expression_associations.csv').set_index('analysis')
    residuals = pd.read_csv(DATA / 'expression_dependency_residuals.csv')
    genes = manifest['partner_order']
    fig = plt.figure(figsize=(8.4, 10.1))
    for col, analysis in enumerate(['raw', 'primary']):
        ax = fig.add_axes([.11 + .47 * col, .727, .355, .235])
        row = statistics.loc[analysis]
        if analysis == 'raw':
            x, y = xy.DDX3Y_expression.to_numpy(), xy.DDX3X.to_numpy()
            intercept = y.mean() - row.effect * x.mean()
            ax.set_xlabel('DDX3Y expression, log₂(TPM + 1)')
            ax.set_ylabel('DDX3X gene effect')
        else:
            vals = residuals.loc[residuals.analysis.eq(analysis)].set_index('ModelID').loc[xy.index]
            x, y = vals.expression_residual.to_numpy(), vals.dependency_residual.to_numpy()
            intercept = 0.
            ax.set_xlabel('DDX3Y expression residual')
            ax.set_ylabel('DDX3X gene-effect residual')
        ax.scatter(x, y, s=11, color='#3C5488', alpha=.40, edgecolors='none')
        grid = np.linspace(x.min(), x.max(), 100)
        ax.plot(grid, intercept + row.effect * grid, color='#222222', lw=1.2)
        ax.margins(x=.06, y=.17)
        p_label = f'{row.p_value:.1e}'.replace('e-0', 'e−').replace('e-', 'e−')
        text = (f"{'Pearson' if analysis == 'raw' else 'Partial'} r = {row.partial_r:.2f}"
                f" [{row.r_ci_low:.2f}, {row.r_ci_high:.2f}]\n"
                f"P = {p_label}; n = {len(xy)}")
        ax.text(.035, .965, text, transform=ax.transAxes, va='top', fontsize=8.3,
                bbox=dict(facecolor='white', edgecolor='none', alpha=.88, pad=2))
        fig.text(.037 + .47 * col, .971, chr(ord('a') + col),
                 fontsize=12, fontweight='semibold')
    cols = ['XY_low_DDX3Y', 'XY_high_DDX3Y', 'LOY', 'Female']
    values = pd.read_csv(DATA / 'stratum_correlations.csv')
    matrix = values.pivot(index='gene', columns='stratum', values='pearson_r').loc[genes, cols]
    interactions = pd.read_csv(DATA / 'expression_partner_interactions.csv')
    q = interactions.loc[interactions.analysis.eq('primary')].set_index('gene').loc[
        genes, 'fdr_bh_31_within_analysis']
    ax = fig.add_axes([.18, .093, .55, .516])
    cmap = LinearSegmentedColormap.from_list('signed_correlation',
                                             ['#3C5488', '#F7F7F7', '#C96C2D'])
    im = ax.imshow(matrix, cmap=cmap, vmin=-.5, vmax=.5, aspect='auto')
    ax.set_yticks(range(len(genes)), genes)
    labels = ['XY: lower\nDDX3Y\nn = 153', 'XY: higher\nDDX3Y\nn = 152',
              'LOY\n\nn = 216', 'Female\n\nn = 448']
    ax.set_xticks(range(4), labels)
    ax.xaxis.tick_top()
    ax.tick_params(axis='both', length=0)
    ax.tick_params(axis='x', pad=7, labelsize=8.4)
    ax.tick_params(axis='y', pad=7, labelsize=8.2)
    for tick, gene in zip(ax.get_yticklabels(), genes):
        if gene in ['GBF1', 'IGBP1', 'IK']:
            tick.set_fontweight('semibold')
    for i, gene in enumerate(genes):
        for j in range(4):
            value = matrix.iloc[i, j]
            ax.text(j, i, f'{value:.2f}', ha='center', va='center',
                    fontsize=7.7, color='white' if abs(value) > .38 else '#222222')
        ax.text(4.35, i, f'{q.loc[gene]:.2f}', ha='center', va='center',
                fontsize=7.7, clip_on=False, color='#555555')
    ax.set_xticks(np.arange(-.5, 4, 1), minor=True)
    ax.set_yticks(np.arange(-.5, len(genes), 1), minor=True)
    ax.grid(which='minor', color='white', lw=.8)
    ax.tick_params(which='minor', bottom=False, left=False, top=False)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.text(4.35, 1.014, 'DDX3Y ×\npartner\nFDR',
            transform=ax.get_xaxis_transform(), ha='center', va='bottom', fontsize=8.3)
    fig.text(.037, .641, 'c', fontsize=12, fontweight='semibold')
    color_ax = fig.add_axes([.30, .046, .31, .012])
    cb = fig.colorbar(im, cax=color_ax, orientation='horizontal', ticks=[-.5, 0, .5])
    cb.outline.set_visible(False)
    cb.ax.tick_params(labelsize=8, length=2)
    cb.set_label('Adjusted Pearson r with DDX3X', fontsize=8.4, labelpad=3)
    for ext in ['pdf', 'svg', 'png']:
        fig.savefig(OUT / f'{STEM}.{ext}', dpi=600 if ext == 'png' else 300)
    review = OUT / 'review'
    review.mkdir(exist_ok=True)
    fig.savefig(review / f'{STEM}_preview.png', dpi=160)
    plt.close(fig)
    with fitz.open(OUT / f'{STEM}.pdf') as doc:
        page = doc[0]
        assert all('OpenSans' in f[3] for f in page.get_fonts())
        spans = [s for b in page.get_text('dict')['blocks'] if 'lines' in b
                 for line in b['lines'] for s in line['spans'] if s['text'].strip()]
        outside = [s['text'] for s in spans if not page.rect.contains(fitz.Rect(s['bbox']))]
        overlaps = []
        for i, a in enumerate(spans):
            for b in spans[i + 1:]:
                ra, rb = fitz.Rect(a['bbox']), fitz.Rect(b['bbox'])
                area = (ra & rb).get_area()
                if area > 1 and area / min(ra.get_area(), rb.get_area()) > .1:
                    overlaps.append([a['text'], b['text']])
        assert not outside, outside
        assert not overlaps, overlaps
        page.get_pixmap(matrix=fitz.Matrix(1.5, 1.5), colorspace=fitz.csGRAY).save(
            review / f'{STEM}_greyscale.png')
    record = dict(stem=STEM, dimensions_inches=[8.4, 10.1], font='Open Sans',
                  scatter_models=len(xy), scatter_donors=int(xy.PatientID.nunique()),
                  heatmap_cells=int(matrix.size), gene_order=genes, columns=cols,
                  scale=[-.5, .5], no_color_clipping=bool(matrix.abs().max().max() <= .5),
                  intervals='1000 donor-bootstrap percentile intervals for Pearson r.',
                  p_values='Donor-clustered two-sided regression tests of expression term.',
                  fdr='BH over 31 continuous expression-by-partner slope-interaction tests.',
                  median_xy_expression=manifest['summary']['xy_expression_median'],
                  validation=dict(text_outside_page=outside, text_overlaps=overlaps))
    (DATA / 'render_manifest.json').write_text(json.dumps(record, indent=2), encoding='utf-8')
    print(json.dumps(record['validation']))


if __name__ == '__main__':
    main()
