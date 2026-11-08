"""Overlay the existing Y-guide slope draft on one transparent axis."""
from pathlib import Path
import json, hashlib
from make_ht1080_y_guide_slopes_draft import OUT
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pymupdf

STEM = 'ch04_supp_ht1080_y_guides_d0_endpoint_overlay_draft'
SOURCE = OUT / 'source_data/ht1080_y_guide_slopes_draft/guide_slope_values.csv'
DATA = OUT / 'source_data/ht1080_y_guide_overlay_draft'
GUIDE_COLOUR = '#555555'


def main():
    data = pd.read_csv(SOURCE)
    assert len(data) == 188 and data.gene.nunique() == 47
    DATA.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.2, 6.4))
    fig.subplots_adjust(left=.13, right=.965, top=.94, bottom=.27)
    ax.axhline(0, color='#AAAAAA', lw=.7, ls='--', zorder=0)
    connected = 0
    # Preserve the previous draft's fixed drawing order.
    for row in data.sample(frac=1, random_state=17).itertuples():
        y = np.array([row.plotted_D0_LFC, row.plotted_endpoint_LFC])
        if np.isfinite(y).all():
            ax.plot([0, 1], y, color=GUIDE_COLOUR, lw=.85, alpha=.15, zorder=1)
            connected += 1
    point_count = 0
    for x, field in [(0, 'plotted_D0_LFC'), (1, 'plotted_endpoint_LFC')]:
        vals = data[field].dropna()
        ax.scatter(np.repeat(x, len(vals)), vals, color=GUIDE_COLOUR, marker='o',
            s=13, alpha=.32, linewidths=0, zorder=2)
        point_count += len(vals)
    ax.set_xlim(-.10, 1.10)
    ax.set_ylim(-3.3, 3.3)
    ax.set_xticks([0, 1], ['D0\nCPM + 1 contrast', 'Endpoint\nMAGeCK'])
    ax.set_yticks(range(-3, 4))
    ax.set_ylabel('Guide log₂ fold change (LOY / sgCTRL)', labelpad=8)
    ax.tick_params(axis='x', length=0, pad=7)
    fig.text(.13, .968, '47 genes · 186 complete guide trajectories', fontsize=9)
    fig.text(.965, .968, 'Rough draft', ha='right', color='#666666', fontsize=8)
    fig.text(.13, .155, 'D0 LFC = log₂[(LOY CPM + 1)/(sgCTRL CPM + 1)]', fontsize=8.5)
    fig.text(.13, .117, 'Different normalisation at D0 and endpoint; descriptive comparison.', fontsize=8.5)
    fig.text(.13, .079, 'Lines: 15% opacity. Each line connects the same guide.', fontsize=8)
    fig.text(.13, .041, 'PRY guide 2: D0 only. XKRY guide 3: no observed LFC at either stage.', fontsize=8)
    assert connected == 186 and point_count == 373
    fig.savefig(OUT / f'{STEM}.pdf')
    fig.savefig(OUT / f'{STEM}.svg')
    fig.savefig(OUT / f'{STEM}.png', dpi=600)
    fig.savefig(OUT / f'review/{STEM}_preview.png', dpi=160)
    plt.close(fig)
    doc = pymupdf.open(OUT / f'{STEM}.pdf')
    fonts = sorted({f[3] for f in doc[0].get_fonts()})
    assert all('OpenSans' in f for f in fonts)
    for block in doc[0].get_text('dict')['blocks']:
        for line in block.get('lines', []):
            for span in line['spans']:
                x0, y0, x1, y1 = span['bbox']
                assert x0 >= 0 and y0 >= 0 and x1 <= doc[0].rect.width and y1 <= doc[0].rect.height, span['text']
    manifest = {
        'status': 'rough single-axis visual draft',
        'input': {'path': str(SOURCE), 'sha256': hashlib.sha256(SOURCE.read_bytes()).hexdigest()},
        'estimators': 'Existing D0 CPM+1 log2 ratio and saved endpoint MAGeCK LFC; different normalisations; not a baseline correction.',
        'colour_semantics': 'One neutral colour for all guides; guide identity is not encoded',
        'line_colour': GUIDE_COLOUR, 'point_colour': GUIDE_COLOUR, 'point_marker': 'o',
        'line_alpha': .15, 'point_alpha': .32, 'drawing_order_seed': 17,
        'figure_inches': [7.2, 6.4], 'y_limits': [-3.3, 3.3], 'font': fonts,
        'validation': {'source_guides': 188, 'genes': 47, 'lines': connected,
            'D0_points': 187, 'endpoint_points': 186, 'total_points': point_count,
            'no_axis_clipping': bool(data[['plotted_D0_LFC', 'plotted_endpoint_LFC']].abs().max().max() < 3.3),
            'no_text_clipping': True}}
    (DATA/'source_manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    print(json.dumps({'figure': str(OUT/f'{STEM}.pdf'), 'validation': manifest['validation']}, indent=2))


if __name__ == '__main__':
    main()
