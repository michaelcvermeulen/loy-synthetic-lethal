"""Standalone supplement from panel b of the chromosome summary heatmap."""
import json
import hashlib
import shutil
from datetime import datetime

import make_ht1080_chromosome_heatmap as original
from matplotlib.colors import LinearSegmentedColormap, Normalize

pd, np, plt = original.pd, original.np, original.plt
OUT = original.OUT
DATA = OUT / 'source_data/ht1080_chromosome_shift_heatmap'
STEM = 'ch04_supp_ht1080_chromosome_shift_heatmap'
COLUMNS = original.COLUMNS[2:]


def main():
    summary = pd.read_csv(original.INPUT, dtype={'chromosome': str})
    assert summary.chromosome.tolist() == original.source.GROUPS
    values = summary[COLUMNS].to_numpy()
    assert values.shape == (25, 2) and np.isfinite(values).all()
    assert np.abs(values).max() <= .4
    DATA.mkdir(parents=True, exist_ok=True)
    existing = [OUT / f'{STEM}.{ext}' for ext in ['pdf', 'svg', 'png']]
    existing += [OUT / f'review/{STEM}_preview.png']
    if any(p.exists() for p in existing):
        archive = OUT / 'archive' / ('chromosome_shift_heatmap_' + datetime.now().strftime('%Y%m%d_%H%M%S'))
        archive.mkdir(parents=True)
        for p in existing:
            if p.exists():
                shutil.copy2(p, archive / p.name)

    colours = ['#3B6FB6', '#F7F7F7', '#D88036']
    cmap = LinearSegmentedColormap.from_list('chromosome_signed_effect', colours)
    # 22 Oct 2026: narrower canvas with the count column beside the heatmap (was 5.3 in wide).
    fig = plt.figure(figsize=(4.2, 7.8))
    ax = fig.add_axes([.183, .27, .448, .66])
    count_ax = fig.add_axes([.655, .27, .214, .66])
    mesh = ax.pcolormesh(np.arange(3), np.arange(26), values,
                        cmap=cmap, norm=Normalize(-.4, .4),
                        edgecolors='white', linewidth=.85)
    ax.set_xlim(0, 2); ax.set_ylim(25, 0)
    ax.set_xticks([.5, 1.5], ['All guides', 'D0 ≥50'])
    ax.xaxis.tick_top()
    ax.set_yticks(np.arange(25)+.5,
                 [f'Chr {c}' if c != 'NTC' else 'NTC' for c in summary.chromosome])
    ax.tick_params(axis='both', length=0, pad=5, labelsize=9)
    ax.axhline(24, color='#777777', lw=.65)
    for spine in ax.spines.values():
        spine.set_visible(False)
    for row in range(25):
        for col in range(2):
            value = values[row, col]
            label = f'{value:+.3f}'.replace('-', '−') if abs(value) >= .0005 else '0.000'
            ax.text(col+.5, row+.5, label, ha='center', va='center',
                    fontsize=8.5, color='#222222')
    count_ax.set_xlim(0, 1); count_ax.set_ylim(25, 0); count_ax.axis('off')
    for row, item in summary.iterrows():
        count_ax.text(.5, row+.5, f'{item.n_paired_guides:,} / {item.n_D0_ge50_both:,}',
                      ha='center', va='center', fontsize=8.3, color='#555555')
    count_ax.axhline(24, color='#BBBBBB', lw=.65)
    fig.text(.762, .952, 'n guides\nall / D0 ≥50', ha='center', va='center',
             fontsize=8.5, color='#555555')

    fig.text(.526, .222, 'Median ΔLFC relative to non-targeting controls',
             ha='center', fontsize=8.8)
    cb_ax = fig.add_axes([.183, .175, .686, .018])
    cb = fig.colorbar(mesh, cax=cb_ax, orientation='horizontal',
                      ticks=[-.4, -.2, 0, .2, .4])
    cb.outline.set_visible(False)
    cb.ax.tick_params(labelsize=8.5, length=2, pad=2)
    fig.text(.183, .115, 'Relative guide depletion\nin LOY vs CTRL',
             ha='left', va='top', fontsize=8.4, color='#305A92')
    fig.text(.869, .115, 'Relative guide enrichment\nin LOY vs CTRL',
             ha='right', va='top', fontsize=8.4, color='#97521F')
    fig.text(.183, .047, 'ΔLFC = endpoint − D0 (log₂ units); 0 = NTC median shift.', fontsize=8)
    fig.text(.183, .024, 'D0 ≥50: at least 50 reads in each D0 sample.', fontsize=8)

    for ext in ['pdf', 'svg', 'png']:
        fig.savefig(OUT / f'{STEM}.{ext}', dpi=600)
    fig.savefig(OUT / f'review/{STEM}_preview.png', dpi=160)
    plt.close(fig)
    summary[['chromosome']+COLUMNS+['n_paired_guides', 'n_D0_ge50_both']].to_csv(DATA/'plotted_values.csv', index=False)
    with original.source.pymupdf.open(OUT/f'{STEM}.pdf') as doc:
        fonts = sorted({font[3] for page in doc for font in page.get_fonts()})
        assert all('OpenSans' in font for font in fonts)
        for page in doc:
            for block in page.get_text('dict')['blocks']:
                for line in block.get('lines', []):
                    for span in line['spans']:
                        x0, y0, x1, y1 = span['bbox']
                        assert x0>=0 and y0>=0 and x1<=page.rect.width and y1<=page.rect.height, span['text']
    manifest = {
        'input': str(original.INPUT),
        'sha256': hashlib.sha256(original.INPUT.read_bytes()).hexdigest(),
        'source_columns_in_display_order': COLUMNS,
        'chromosome_order': summary.chromosome.tolist(),
        'colour_scale': {'colours': colours, 'limits': [-.4, .4], 'centre': 0, 'clipped_cells': 0},
        'direction_labels': {'negative': 'Relative guide depletion in LOY vs CTRL',
                             'positive': 'Relative guide enrichment in LOY vs CTRL'},
        'statistic': 'Median within-guide endpoint-minus-D0 LFC minus the eligible NTC median delta, separately for each guide set.',
        'fonts': fonts, 'unchanged_values_from_original_panel_b': True,
        'new_statistical_tests': False,
        'validation': {'all_50_cells_finite': True, 'no_text_outside_page': True},
    }
    (DATA/'source_manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    print(OUT/f'{STEM}.pdf')


if __name__ == '__main__':
    main()
