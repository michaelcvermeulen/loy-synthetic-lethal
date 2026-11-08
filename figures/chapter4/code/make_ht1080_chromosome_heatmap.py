"""Restyle saved chromosome summaries as an annotated heatmap; no refitting."""
from pathlib import Path
import json
import hashlib
import shutil
from datetime import datetime

import make_ht1080_chromosome_guide_bias as source
from matplotlib.colors import LinearSegmentedColormap, Normalize

pd, np, plt = source.pd, source.np, source.plt
OUT = source.OUT
DATA = OUT / 'source_data/ht1080_chromosome_heatmap'
STEM = 'ch04_supp_ht1080_chromosome_guide_bias_heatmap'
INPUT = source.DATA / 'chromosome_summary.csv'
COLUMNS = [
    'median_D0_LFC', 'median_endpoint_matched_LFC',
    'median_delta_relative_NTC', 'median_delta_D0_ge50_relative_NTC',
]


def main():
    DATA.mkdir(parents=True, exist_ok=True)
    summary = pd.read_csv(INPUT, dtype={'chromosome': str})
    assert summary.chromosome.tolist() == source.GROUPS
    values = summary[COLUMNS].to_numpy()
    assert values.shape == (25, 4) and np.isfinite(values).all()
    limit = .4
    assert np.abs(values).max() <= limit
    existing = [OUT / f'{STEM}.{ext}' for ext in ['pdf', 'svg', 'png']]
    existing += [OUT / f'review/{STEM}_preview.png']
    if any(p.exists() for p in existing):
        archive = OUT / 'archive' / ('chromosome_heatmap_' + datetime.now().strftime('%Y%m%d_%H%M%S'))
        archive.mkdir(parents=True)
        for p in existing:
            if p.exists():
                shutil.copy2(p, archive / p.name)

    colours = ['#3B6FB6', '#F7F7F7', '#D88036']
    cmap = LinearSegmentedColormap.from_list('chromosome_signed_effect', colours)
    norm = Normalize(vmin=-limit, vmax=limit)
    fig = plt.figure(figsize=(6.5, 8.3))
    a = fig.add_axes([.12, .185, .25, .73])
    b = fig.add_axes([.415, .185, .25, .73])
    count_ax = fig.add_axes([.715, .185, .26, .73])
    for ax, block, labels in [
        (a, values[:, :2], ['D0', 'Endpoint']),
        (b, values[:, 2:], ['All guides', 'D0 ≥50']),
    ]:
        mesh = ax.pcolormesh(np.arange(3), np.arange(26), block,
                            cmap=cmap, norm=norm, edgecolors='white', linewidth=.85)
        ax.set_xlim(0, 2)
        ax.set_ylim(25, 0)
        ax.set_xticks([.5, 1.5], labels)
        ax.xaxis.tick_top()
        ax.tick_params(axis='both', length=0, pad=5, labelsize=9)
        ax.set_yticks(np.arange(25)+.5)
        ax.axhline(24, color='#777777', lw=.65)
        for spine in ax.spines.values():
            spine.set_visible(False)
        for row in range(25):
            for col in range(2):
                value = block[row, col]
                label = f'{value:+.3f}' if abs(value) >= .0005 else '0.000'
                rgba = cmap(norm(value))
                # Numeric labels retain magnitude and sign when viewed in greyscale.
                luminance = np.dot(rgba[:3], [.2126, .7152, .0722])
                ax.text(col+.5, row+.5, label, ha='center', va='center',
                        fontsize=8.4, color='white' if luminance < .48 else '#222222')
    a.set_yticklabels([f'Chr {c}' if c != 'NTC' else 'NTC' for c in summary.chromosome])
    b.set_yticklabels([])
    count_ax.set_xlim(0, 1)
    count_ax.set_ylim(25, 0)
    count_ax.axis('off')
    for row, item in summary.iterrows():
        count_ax.text(.5, row+.5, f'{item.n_paired_guides:,} / {item.n_D0_ge50_both:,}',
                      ha='center', va='center', fontsize=8.3, color='#555555')
    count_ax.axhline(24, color='#BBBBBB', lw=.65)
    fig.text(.12, .966, 'a', fontsize=13, weight='semibold')
    fig.text(.16, .966, 'Median guide LFC', fontsize=9, weight='semibold')
    fig.text(.415, .966, 'b', fontsize=13, weight='semibold')
    fig.text(.455, .966, 'ΔLFC relative to NTC', fontsize=9, weight='semibold')
    fig.text(.845, .938, 'n guides\nall / D0 ≥50', ha='center', va='center', fontsize=8.5, color='#555555')
    cb_ax = fig.add_axes([.18, .122, .425, .015])
    cb = fig.colorbar(mesh, cax=cb_ax, orientation='horizontal', ticks=[-.4, -.2, 0, .2, .4])
    cb.outline.set_visible(False)
    cb.ax.tick_params(labelsize=8, length=2, pad=2)
    cb.set_label('log₂ units (shared scale)', fontsize=8.5, labelpad=2)
    fig.text(.12, .045, 'ΔLFC = endpoint − D0; centred on NTCs within each guide set.', fontsize=8)
    fig.text(.12, .023, 'D0 ≥50: at least 50 reads in both D0 samples. Descriptive medians.', fontsize=8)

    for ext in ['pdf', 'svg', 'png']:
        fig.savefig(OUT / f'{STEM}.{ext}', dpi=600)
    fig.savefig(OUT / f'review/{STEM}_preview.png', dpi=160)
    plt.close(fig)
    plotted = summary[['chromosome'] + COLUMNS + ['n_paired_guides', 'n_D0_ge50_both']]
    plotted.to_csv(DATA / 'plotted_values.csv', index=False)
    with source.pymupdf.open(OUT / f'{STEM}.pdf') as doc:
        assert len(doc) == 1
        fonts = sorted({f[3] for page in doc for f in page.get_fonts()})
        assert all('OpenSans' in f for f in fonts)
        for page in doc:
            for block in page.get_text('dict')['blocks']:
                for line in block.get('lines', []):
                    for span in line['spans']:
                        x0, y0, x1, y1 = span['bbox']
                        assert x0 >= 0 and y0 >= 0 and x1 <= page.rect.width and y1 <= page.rect.height, span['text']
    manifest = {
        'input': str(INPUT),
        'sha256': hashlib.sha256(INPUT.read_bytes()).hexdigest(),
        'source_columns_in_display_order': COLUMNS,
        'chromosome_order': summary.chromosome.tolist(),
        'colour_scale': {'colours': colours, 'limits': [-limit, limit], 'centre': 0, 'clipped_cells': 0},
        'cell_label_decimal_places': 3,
        'fonts': fonts,
        'original_analytical_values_retained': True,
        'new_statistical_tests': False,
        'panel_a': 'Stage-wise median guide LFCs using log2[(LOY CPM+1)/(sgCTRL CPM+1)].',
        'panel_b': 'Median within-guide endpoint-minus-D0 LFC minus the eligible NTC median delta, separately for each guide set. This is not the difference of the two panel-a medians.',
        'guide_counts': 'Observed at both stages; second count also requires at least 50 raw reads in each D0 sample.',
        'validation': {'all_100_cells_finite': True, 'shared_symmetric_limits_no_clipping': True, 'no_text_outside_page': True},
    }
    (DATA / 'source_manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    print(json.dumps({'pdf': str(OUT / f'{STEM}.pdf'), 'cells': 100, 'range': [float(values.min()), float(values.max())]}, indent=2))


if __name__ == '__main__':
    main()
