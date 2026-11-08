"""Rough visual comparison of the currently displayed D0 and endpoint guide LFCs."""
from pathlib import Path
import sys, os, json, hashlib

THESIS = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(THESIS.parent / '.figure_python_deps'))
OUT = THESIS / 'figures/chapter4'
SOURCE = OUT / 'source_data/ht1080_y_guides'
DATA = OUT / 'source_data/ht1080_y_guide_slopes_draft'
DATA.mkdir(parents=True, exist_ok=True)
os.environ.setdefault('MPLCONFIGDIR', str(DATA / 'mpl_cache'))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.backends.backend_pdf import PdfPages
import pymupdf
sys.path.insert(0, str(THESIS / 'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
apply_style(9)

STEM = 'ch04_supp_ht1080_y_guides_d0_endpoint_slopes_draft'
# Guide identity, not genotype. Shapes provide a redundant identity cue.
GUIDE_COLOURS = {1: '#009E73', 2: '#CC79A7', 3: '#E69F00', 4: '#56B4E9'}
GUIDE_MARKERS = {1: 'o', 2: '^', 3: 's', 4: 'D'}
KEY_GENES = ['DDX3Y', 'EIF1AY', 'RPS4Y1', 'RPS4Y2', 'ZFY', 'KDM5D', 'UTY', 'USP9Y']


def prepare():
    source = pd.read_csv(SOURCE / 'all_y_guides.csv')
    order = json.loads((SOURCE / 'summary.json').read_text())['gene_order']
    keep = ['gene', 'guide_number', 'guide_id', 'sequence', 'D0_CTRL', 'D0_LOY',
        'D0_CTRL_CPM', 'D0_LOY_CPM', 'D0_log2FC_CPM_plus_1',
        'endpoint_mean_sgCTRL_CPM', 'endpoint_mean_LOY_CPM', 'mageck_LFC',
        'D0_zero_both', 'endpoint_zero', 'endpoint_samples_detected']
    plotted = source[keep].copy()
    plotted['plotted_D0_LFC'] = plotted.D0_log2FC_CPM_plus_1
    plotted['plotted_endpoint_LFC'] = plotted.mageck_LFC.mask(plotted.endpoint_zero)
    plotted['connected'] = plotted[['plotted_D0_LFC', 'plotted_endpoint_LFC']].notna().all(axis=1)
    plotted.to_csv(DATA / 'guide_slope_values.csv', index=False)
    assert len(plotted) == 188 and plotted.groupby('gene').size().eq(4).all()
    assert plotted.connected.sum() == 186
    assert plotted.plotted_D0_LFC.notna().sum() == 187
    assert plotted.plotted_endpoint_LFC.notna().sum() == 186
    assert plotted[['plotted_D0_LFC', 'plotted_endpoint_LFC']].abs().max().max() < 3.3
    return plotted, order


def draw_page(data, genes, page_number):
    fig, axes = plt.subplots(6, 4, figsize=(10.8, 11.8), sharey=True)
    fig.subplots_adjust(left=.07, right=.985, top=.918, bottom=.12, wspace=.17, hspace=.67)
    line_count, dot_count = 0, 0
    for i, ax in enumerate(axes.flat):
        if i >= len(genes):
            ax.set_visible(False)
            continue
        gene = genes[i]
        rows = data.loc[data.gene.eq(gene)].sort_values('guide_number')
        ax.axhline(0, color='#BBBBBB', lw=.65, ls='--', zorder=0)
        for row in rows.itertuples():
            n = row.guide_number
            vals = np.array([row.plotted_D0_LFC, row.plotted_endpoint_LFC], dtype=float)
            finite = np.isfinite(vals)
            if finite.all():
                ax.plot([0, 1], vals, color=GUIDE_COLOURS[n], lw=1.25, alpha=.85, zorder=2)
                line_count += 1
            ax.scatter(np.array([0, 1])[finite], vals[finite], color=GUIDE_COLOURS[n],
                marker=GUIDE_MARKERS[n], s=20, edgecolor='white', lw=.35, zorder=3)
            dot_count += int(finite.sum())
        if gene == 'PRY':
            ax.text(.98, .97, 'G2*: D0 only', transform=ax.transAxes, ha='right', va='top', fontsize=7, color='#666666')
        if gene == 'XKRY':
            ax.text(.98, .97, 'G3*: no reads', transform=ax.transAxes, ha='right', va='top', fontsize=7, color='#666666')
        ax.set_xlim(-.17, 1.17)
        ax.set_ylim(-3.3, 3.3)
        ax.set_xticks([0, 1], ['D0', 'Endpoint'])
        ax.set_yticks([-3, 0, 3])
        ax.tick_params(labelsize=8, length=2.5, pad=2)
        if i % 4: ax.tick_params(labelleft=False)
        ax.set_title(gene, fontsize=9, pad=4, weight='semibold' if gene in KEY_GENES else 'normal')
    handles = [Line2D([], [], color=GUIDE_COLOURS[n], marker=GUIDE_MARKERS[n], lw=1.2,
        markersize=4, label=f'Guide {n}') for n in range(1, 5)]
    fig.legend(handles=handles, loc='upper left', bbox_to_anchor=(.065, .994), ncol=4,
        handlelength=1.8, columnspacing=1.5)
    fig.text(.985, .976, f'Rough draft · {page_number}/2', ha='right', fontsize=9, color='#666666')
    fig.text(.07, .948, 'Each line connects the same guide; guide numbers match the heatmaps.', fontsize=9)
    fig.supylabel('Guide log₂ fold change (LOY / sgCTRL)', x=.012, fontsize=10)
    fig.text(.07, .070, 'D0: log₂[(LOY CPM + 1)/(sgCTRL CPM + 1)]     Endpoint: saved MAGeCK LFC', fontsize=8.5)
    fig.text(.07, .050, 'Normalisation differs between D0 and endpoint; slopes are a visual comparison, not baseline-corrected effects.', fontsize=8.5)
    fig.text(.07, .030, '* PRY guide 2 has D0 only; XKRY guide 3 has no observed LFC at either stage. Missing points are omitted.', fontsize=8)
    fig.text(.07, .013, 'D0: one sample per genotype. Endpoint: two outgrowth cultures from one population per genotype.', fontsize=8)
    return fig, {'genes': genes, 'connected_guides': line_count, 'points': dot_count}


def main():
    data, order = prepare()
    review = OUT / 'review'
    review.mkdir(exist_ok=True)
    pages = []
    with PdfPages(OUT / f'{STEM}.pdf') as report:
        for part in range(2):
            fig, record = draw_page(data, order[part*24:(part+1)*24], part+1)
            report.savefig(fig)
            stem = f'{STEM}_{part+1}'
            fig.savefig(OUT / f'{stem}.pdf')
            fig.savefig(OUT / f'{stem}.svg')
            fig.savefig(OUT / f'{stem}.png', dpi=600)
            fig.savefig(review / f'{stem}_preview.png', dpi=160)
            plt.close(fig)
            pages.append(record)
    assert sum(p['connected_guides'] for p in pages) == 186
    assert sum(p['points'] for p in pages) == 373
    doc = pymupdf.open(OUT / f'{STEM}.pdf')
    text = '\n'.join(p.get_text() for p in doc)
    assert len(doc) == 2 and all(g in text for g in order)
    fonts = sorted({font[3] for page in doc for font in page.get_fonts()})
    assert all('OpenSans' in font for font in fonts)
    clipped = []
    for p, page in enumerate(doc):
        for block in page.get_text('dict')['blocks']:
            for line in block.get('lines', []):
                for span in line['spans']:
                    x0, y0, x1, y1 = span['bbox']
                    if x0 < 0 or y0 < 0 or x1 > page.rect.width or y1 > page.rect.height:
                        clipped.append([p+1, span['text']])
    assert not clipped, clipped
    sources = {name: {'path': str(SOURCE/name), 'sha256': hashlib.sha256((SOURCE/name).read_bytes()).hexdigest()}
        for name in ['all_y_guides.csv', 'summary.json']}
    manifest = {'status': 'rough visual draft', 'source_files': sources,
        'D0': 'Existing D0_log2FC_CPM_plus_1; one sample per genotype; both-zero unavailable.',
        'endpoint': 'Existing saved mageck_LFC; two outgrowth cultures per genotype; all-zero endpoint guides omitted.',
        'interpretation': 'Two existing, differently normalised estimators. Connecting segments are not a matched-normalisation time-course change or baseline correction. No new significance tests.',
        'guide_colours': GUIDE_COLOURS, 'guide_markers': GUIDE_MARKERS,
        'axes_limits': [-3.3, 3.3], 'font': fonts, 'figure_inches': [10.8, 11.8],
        'pages': pages, 'validation': {'guides_retained': 188, 'genes': 47, 'connected_guides': 186,
            'D0_points': 187, 'endpoint_points': 186, 'total_points': 373, 'clipped_text': False}}
    (DATA / 'source_manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    print(json.dumps({'report': str(OUT/f'{STEM}.pdf'), 'validation': manifest['validation']}, indent=2))


if __name__ == '__main__':
    main()
