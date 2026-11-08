"""Forest plots of all 31 partners, with DDX3X on a separate reference scale."""
from pathlib import Path
import json
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pymupdf as fitz

THESIS = Path(__file__).resolve().parents[3]
OUT = THESIS / 'figures/chapter4'
DATA = OUT / 'source_data/ddx3x_partner_differential_dependency'
STEM = 'ch04_supp_ddx3x_partner_differential_dependency_draft'
sys.path.insert(0, str(THESIS / 'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style


def points(ax, rows, is_reference=False):
    for y, row in enumerate(rows.itertuples()):
        passed = row.fdr_bh_genomewide_joint < .05
        ax.plot([row.ci_low, row.ci_high], [y, y], color='#7A3E85' if passed else '#777777',
                lw=1.1, solid_capstyle='round', zorder=2)
        ax.plot(row.effect, y, marker='o', ms=4.2, linestyle='none',
                markerfacecolor='#7A3E85' if passed else 'white',
                markeredgecolor='#7A3E85' if passed else '#555555',
                markeredgewidth=.8, zorder=3)
        q = row.fdr_bh_genomewide_joint
        label = f'{q:.1g}' if q < .001 else f'{q:.2f}'
        ax.text(1.035, y, label, transform=ax.get_yaxis_transform(),
                va='center', ha='left', fontsize=7.4)
        if not is_reference and y % 2 == 0:
            ax.axhspan(y - .5, y + .5, color='#F4F4F4', zorder=0)
    ax.axvline(0, color='#777777', lw=.7, ls='--', zorder=1)
    ax.spines['left'].set_visible(False)
    ax.tick_params(axis='y', length=0, pad=7)
    ax.tick_params(axis='x', labelsize=8, pad=3)


def main():
    apply_style(8.5)
    data = pd.read_csv(DATA / 'panel_primary_results.csv')
    manifest = json.loads((DATA / 'manifest.json').read_text())
    genes = manifest['panel_genes'][1:]
    assert len(genes) == 31 and len(data) == 96
    specs = [('LOY_minus_XY', 'LOY − XY'),
             ('Female_minus_XY', 'Female − XY'),
             ('LOY_minus_Female', 'LOY − Female')]
    width, height = 10.1, 10.2
    fig = plt.figure(figsize=(width, height))
    lefts = [.105, .405, .705]
    axes_width = .235
    lim = .15
    assert data.loc[data.is_partner, ['ci_low', 'ci_high']].abs().max().max() < lim
    plotted = []
    for col, ((contrast, title), left) in enumerate(zip(specs, lefts)):
        ax = fig.add_axes([left, .285, axes_width, .625])
        rows = data.loc[data.contrast.eq(contrast) & data.is_partner].set_index(
            'gene').loc[genes].reset_index()
        assert rows.gene.tolist() == genes
        points(ax, rows)
        ax.set(xlim=(-lim, lim), ylim=(len(genes) - .5, -.5),
               xticks=[-.1, 0, .1], yticks=np.arange(len(genes)),
               yticklabels=genes if col == 0 else [''] * len(genes))
        ax.set_xlabel('Adjusted difference', labelpad=5)
        ax.text(.5, 1.035, title, transform=ax.transAxes, ha='center',
                va='bottom', fontsize=10, fontweight='semibold')
        ax.text(1.035, 1.015, 'FDR', transform=ax.transAxes, ha='left',
                va='bottom', fontsize=7.4)
        fig.text(left - .045, .947, chr(ord('a') + col), fontsize=12,
                 fontweight='semibold')
        if col == 0:
            for tick, gene in zip(ax.get_yticklabels(), genes):
                if gene in ['IGBP1', 'GBF1', 'IK']:
                    tick.set_fontweight('semibold')
        ref_ax = fig.add_axes([left, .125, axes_width, .055])
        reference = data.loc[data.contrast.eq(contrast) & data.is_reference_ddx3x]
        assert len(reference) == 1
        points(ref_ax, reference, is_reference=True)
        ref_ax.set(xlim=(-1, .2), ylim=(.6, -.6), xticks=[-1, -.5, 0],
                   yticks=[0], yticklabels=['DDX3X'] if col == 0 else [''])
        plotted.extend(rows.gene.tolist() + ['DDX3X'])
    fig.text(.105, .207, 'DDX3X reference (separate x scale)',
             fontsize=8.5, fontweight='semibold')
    fig.text(.52, .075, 'Difference in Chronos gene effect: negative = stronger dependency in the first group',
             ha='center', fontsize=8.2)
    handles = [
        Line2D([], [], marker='o', color='#777777', markerfacecolor='white',
               lw=1, ms=4, label='Genome-wide FDR ≥0.05'),
        Line2D([], [], marker='o', color='#7A3E85', markerfacecolor='#7A3E85',
               lw=1, ms=4, label='Genome-wide FDR <0.05'),
    ]
    fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(.52, .025),
               ncol=2, fontsize=8.2, columnspacing=2)
    for ext in ['pdf', 'svg', 'png']:
        fig.savefig(OUT / f'{STEM}.{ext}', dpi=600 if ext == 'png' else 300)
    fig.savefig(OUT / 'review' / f'{STEM}_preview.png', dpi=170)
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
            OUT / 'review' / f'{STEM}_greyscale.png')
    record = dict(
        stem=STEM, dimensions_inches=[width, height], font='Open Sans',
        partner_order=genes, contrasts=[c for c, _ in specs], plotted_estimates=len(plotted),
        partner_axis_limits=[-lim, lim], reference_axis_limits=[-1, .2],
        uncertainty='Nominal 95% donor-clustered t intervals; not multiplicity-adjusted.',
        significance='BH across 17,787 genes x 3 contrasts = 53,361 tests.',
        annotation='Actual genome-wide FDR shown to two decimals, or one significant '
                   'figure below 0.001. Original three partner labels are bold.',
        validation=dict(all_estimates_plotted=len(plotted) == 96,
                        no_clipped_intervals=True, text_outside_page=outside,
                        text_overlaps=overlaps))
    (DATA / 'render_manifest.json').write_text(json.dumps(record, indent=2),
                                               encoding='utf-8')
    print(json.dumps(record['validation'], indent=2))


if __name__ == '__main__':
    main()
