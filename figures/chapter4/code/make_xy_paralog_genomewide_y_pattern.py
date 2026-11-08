"""Compare matching-partner counts across all 16 X-linked gametolog anchors."""
from pathlib import Path
import json
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.patches import Rectangle
import pymupdf as fitz

THESIS = Path(__file__).resolve().parents[3]
OUT = THESIS / 'figures/chapter4'
DATA = OUT / 'source_data/xy_paralog_genomewide_y_pattern'
STEM = 'ch04_supp_xy_paralog_y_pattern_counts_draft'
sys.path.insert(0, str(THESIS / 'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style


def main():
    apply_style(9)
    manifest = json.loads((DATA / 'manifest.json').read_text())
    summary = pd.read_csv(DATA / 'anchor_summary.csv').set_index('anchor')
    sensitivity = pd.read_csv(DATA / 'threshold_sensitivity_counts.csv')
    fixed = sensitivity.loc[sensitivity.max_abs_xy_r.eq(.15)
                            & sensitivity.min_delta.eq(.20)]
    counts = fixed.pivot(index='anchor', columns='min_positive_r',
                         values='selected_partners').loc[manifest['anchors'], [.20, .25, .30]]
    assert counts[.25].tolist() == summary.loc[counts.index, 'selected_partners'].tolist()
    counts.rename(columns={.20: 'min_r_020', .25: 'min_r_025',
                           .30: 'min_r_030'}).to_csv(DATA / 'displayed_counts.csv')
    width, height = 5.1, 5.9
    fig = plt.figure(figsize=(width, height))
    ax = fig.add_axes([.235, .105, .53, .79])
    mat = counts.to_numpy()
    cmap = plt.get_cmap('magma_r')
    norm = Normalize(vmin=0, vmax=max(1, mat.max()))
    im = ax.imshow(mat, cmap=cmap, norm=norm, aspect='auto', interpolation='none')
    np.testing.assert_array_equal(im.get_array(), mat)
    ax.set(yticks=np.arange(len(counts)), yticklabels=counts.index,
           xticks=np.arange(3), xticklabels=['≥0.20', '≥0.25', '≥0.30'])
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')
    ax.set_xlabel('Minimum adjusted r in both X0 and XX', labelpad=9)
    ax.tick_params(axis='both', which='both', length=0, pad=6)
    ax.set_xticks(np.arange(-.5, 3, 1), minor=True)
    ax.set_yticks(np.arange(-.5, len(counts), 1), minor=True)
    ax.grid(which='minor', color='white', linewidth=1)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.add_patch(Rectangle((.5, -.5), 1, len(counts), fill=False,
                           edgecolor='#444444', linewidth=1.2, clip_on=False))
    ax.get_xticklabels()[1].set_fontweight('semibold')
    for label, anchor in zip(ax.get_yticklabels(), counts.index):
        if anchor == 'DDX3X':
            label.set_fontweight('semibold')
    for i in range(len(counts)):
        for j in range(3):
            rgba = cmap(norm(mat[i, j]))
            lum = np.dot(rgba[:3], [.2126, .7152, .0722])
            ax.text(j, i, str(mat[i, j]), ha='center', va='center',
                    color='white' if lum < .5 else '#222222', fontsize=9)
    cax = fig.add_axes([.82, .38, .027, .28])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label('Matching partners', labelpad=7)
    cb.ax.tick_params(labelsize=8, width=.6, length=2)
    cb.outline.set_linewidth(.5)
    fig.text(.5, .047, '|rXY| ≤0.15; both increases over XY ≥0.20',
             ha='center', va='center', fontsize=8)
    for ext in ['pdf', 'svg', 'png']:
        fig.savefig(OUT / f'{STEM}.{ext}', dpi=600 if ext == 'png' else 300)
    (OUT / 'review').mkdir(exist_ok=True)
    fig.savefig(OUT / 'review' / f'{STEM}_preview.png', dpi=180)
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
        page.get_pixmap(matrix=fitz.Matrix(2, 2), colorspace=fitz.csGRAY).save(
            OUT / 'review' / f'{STEM}_greyscale.png')
    record = dict(
        stem=STEM, dimensions_inches=[width, height], font='Open Sans',
        matrix=counts.reset_index().to_dict('records'),
        row_order=manifest['anchors'],
        colour_scale=dict(palette='magma_r', minimum=0, maximum=int(mat.max()),
                          clipped_values=0),
        primary_column='Minimum adjusted r >=0.25; outlined and bold header.',
        inference='Counts from descriptive effect-size rules, not significant discoveries.',
        validation=dict(all_plotted_counts_match=True, text_outside_page=outside,
                        text_overlaps=overlaps))
    (DATA / 'render_manifest.json').write_text(json.dumps(record, indent=2),
                                              encoding='utf-8')
    print(counts.to_string())


if __name__ == '__main__':
    main()
