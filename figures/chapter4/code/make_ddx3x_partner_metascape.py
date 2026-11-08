"""Thesis-style display of the original full-31-gene Metascape enrichment."""
from pathlib import Path
import json
import sys
import textwrap
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
import pymupdf as fitz

THESIS = Path(__file__).resolve().parents[3]
OUT = THESIS / 'figures/chapter4'
DATA = OUT / 'source_data/ddx3x_partner_metascape'
STEM = 'ch04_supp_ddx3x_partner_metascape_draft'
sys.path.insert(0, str(THESIS / 'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style


def main():
    apply_style(9)
    data = pd.read_csv(DATA / 'cluster_representatives.csv')
    assert len(data) == 14 and data.GROUP_ID.is_unique
    assert data['#GeneInHitList'].eq(31).all()
    width, height = 12.4, 7.5
    fig = plt.figure(figsize=(width, height))
    grid = fig.add_gridspec(1, 4, left=.018, right=.985, top=.88, bottom=.12,
                           width_ratios=[3.45, 1.85, .90, 4.45], wspace=.10)
    axes = [fig.add_subplot(grid[0, i]) for i in range(4)]
    for ax in axes:
        ax.set_ylim(len(data) - .5, -.5)
        ax.set_yticks([])
        for row in range(len(data)):
            if row % 2 == 0:
                ax.axhspan(row - .5, row + .5, color='#f4f4f4', zorder=0)
    for i in [0, 2, 3]:
        axes[i].set_xlim(0, 1)
        axes[i].axis('off')
    ax = axes[1]
    scores = -data['Log(q-value)'].to_numpy()
    strong = data.reported_q_clearly_lt_005.to_numpy(bool)
    colours = np.where(strong, '#7A3E85', '#A9A9A9')
    bars = ax.barh(np.arange(len(data)), scores, height=.43,
                   color=colours, edgecolor='none', zorder=3)
    ax.axvline(-np.log10(.05), color='#666666', lw=.8, ls='--', zorder=4)
    ax.set_xlim(0, max(scores) * 1.06)
    ax.xaxis.set_major_locator(MaxNLocator(4))
    ax.set_xlabel('−log₁₀(q)', labelpad=6)
    ax.spines['left'].set_visible(False)
    ax.tick_params(axis='y', left=False)
    for row, (_, term) in enumerate(data.iterrows()):
        label = term.Description.strip()
        axes[0].text(.99, row, textwrap.fill(label, 43), ha='right', va='center',
                     fontsize=8.5, linespacing=1.15)
        axes[2].text(.5, row, f'{term.approx_q_value:.2g}', ha='center', va='center',
                     fontsize=8.5)
        genes = term.input_gene_symbols.split('|')
        assert len(genes) == term['#GeneInGOAndHitList']
        text = f'{len(genes)}/31: ' + ', '.join(genes)
        axes[3].text(.02, row, textwrap.fill(text, 64), ha='left', va='center',
                     fontsize=8.1, linespacing=1.15)
        assert abs(bars[row].get_width() - scores[row]) < 1e-12
    for i, title in [(0, 'Term'), (1, 'Enrichment significance'),
                     (2, 'q ≈'), (3, 'Contributing genes')]:
        axes[i].text(.99 if i == 0 else (.02 if i == 3 else .5), 1.025,
                     title, transform=axes[i].transAxes,
                     ha='right' if i == 0 else ('left' if i == 3 else 'center'),
                     va='bottom', fontsize=9, fontweight='semibold')
    handles = [
        Line2D([], [], marker='s', linestyle='none', color='#7A3E85',
               label='Reported q <0.05', markersize=6),
        Line2D([], [], marker='s', linestyle='none', color='#A9A9A9',
               label='Reported q ≥0.05', markersize=6),
    ]
    fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(.52, .98),
               ncol=2, fontsize=9, columnspacing=2)
    fig.text(.5, .035, 'Original 31-gene enrichment; q values retain the export’s rounding.',
             ha='center', fontsize=8, color='#555555')
    for ext in ['pdf', 'svg', 'png']:
        fig.savefig(OUT / f'{STEM}.{ext}', dpi=600 if ext == 'png' else 300)
    (OUT / 'review').mkdir(exist_ok=True)
    fig.savefig(OUT / 'review' / f'{STEM}_preview.png', dpi=160)
    plt.close(fig)
    with fitz.open(OUT / f'{STEM}.pdf') as doc:
        page = doc[0]
        assert all('OpenSans' in f[3] for f in page.get_fonts())
        spans = [s for b in page.get_text('dict')['blocks'] if 'lines' in b
                 for line in b['lines'] for s in line['spans'] if s['text'].strip()]
        outside = [s['text'] for s in spans
                   if not page.rect.contains(fitz.Rect(s['bbox']))]
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
        stem=STEM, font='Open Sans', dimensions_inches=[width, height],
        statistical_source='Original GO_AllLists.csv, all 31 genes',
        displayed_terms=data.GO.tolist(),
        source_row_indices=data.source_row.tolist(),
        colour_mapping={'reported_q_lt_005': '#7A3E85', 'reported_q_ge_005': '#A9A9A9'},
        validation={'bar_values_match': True, 'overlap_counts_match': True,
                    'text_outside_page': outside, 'text_overlaps': overlaps})
    (DATA / 'render_manifest.json').write_text(json.dumps(record, indent=2),
                                               encoding='utf-8')
    print(json.dumps(record, indent=2))


if __name__ == '__main__':
    main()
