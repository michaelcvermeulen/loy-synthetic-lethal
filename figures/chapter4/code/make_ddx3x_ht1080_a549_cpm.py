"""Compare all four DDX3X guides in HT1080 and retained A549 libraries."""
from pathlib import Path
import hashlib
import json
import os
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT.parent / '.figure_python_deps'))
OUT = ROOT / 'figures/chapter4'
DATA = OUT / 'source_data/ddx3x_ht1080_a549_cpm'
DATA.mkdir(parents=True, exist_ok=True)
os.environ.setdefault('MPLCONFIGDIR', str(DATA / 'mpl_cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter
import numpy as np
import pandas as pd
import pymupdf

sys.path.insert(0, str(ROOT / 'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style, condition_colours
apply_style(8)
PALETTE = condition_colours()
COLOURS = {'CTRL': PALETTE['Y-retaining'], 'LOY': PALETTE['Y-depleted']}
STEM = 'ch04_supp_ddx3x_ht1080_a549_guide_cpm_draft'
SOURCES = {
    'HT1080': ROOT.parent / 'crispr_screen/ht1080_loy_screen_2/count/ht1080_screen2_R1.count.txt',
    'A549': ROOT.parent / 'crispr_screen/a549_loy_screen1/CRISPR_screen/Screen1/Mageck/A549_full.count.txt',
}
SAMPLES = {
    'HT1080': ['D0_CTRL', 'D0_LOY', 'POST_CTRL_REP1', 'POST_CTRL_REP2', 'POST_LOY_REP1', 'POST_LOY_REP2'],
    'A549': ['CTRL_D0_1', 'CTRL_D0_2', 'LOY_D0_1', 'LOY_D0_2', 'CTRL_P4_R2', 'LOY_P4_R1', 'LOY_P4_R2'],
}
LABELS = {
    'HT1080': ['CTRL\nD0', 'LOY\nD0', 'CTRL\nPOST-R1', 'CTRL\nPOST-R2', 'LOY\nPOST-R1', 'LOY\nPOST-R2'],
    'A549': ['CTRL\nD0-1', 'CTRL\nD0-2', 'LOY\nD0-1', 'LOY\nD0-2', 'CTRL\nP4-R2', 'LOY\nP4-R1', 'LOY\nP4-R2'],
}
GUIDES = [
    'DDX3X_4637_CGTGGACGGAGTGATTACGA',
    'DDX3X_4638_GCACCACCATAAACCACGCA',
    'DDX3X_4639_TAGCACTAAACATCATAGTG',
    'DDX3X_4640_TCTCGGTTCCTTAAATGAGG',
]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prepare():
    hashes = {line: digest(path) for line, path in SOURCES.items()}
    records, totals_records = [], []
    for line, source in SOURCES.items():
        raw = pd.read_csv(source, sep='\t').set_index('sgRNA')
        assert len(raw) == 77441 and raw.index.is_unique
        assert sorted(raw.index[raw.Gene.eq('DDX3X')]) == GUIDES
        samples = SAMPLES[line]
        assert not any('_P3_' in sample for sample in samples)
        totals = raw[samples].sum()
        assert (totals > 0).all()
        for sample in samples:
            totals_records.append(dict(cell_line=line, sample=sample,
                                       total_guide_reads=int(totals[sample]), denominator_guides=len(raw)))
        for guide in GUIDES:
            for sample in samples:
                count = int(raw.loc[guide, sample])
                records.append(dict(cell_line=line, sgrna=guide, guide=guide.split('_')[1],
                                    sequence=guide.rsplit('_', 1)[1], sample=sample,
                                    genotype='CTRL' if 'CTRL' in sample else 'LOY',
                                    stage='D0' if 'D0' in sample else 'Endpoint',
                                    raw_count=count, cpm=count / totals[sample] * 1e6))
    values = pd.DataFrame(records)
    assert len(values) == 52 and (values.cpm >= 0).all()
    old = pd.read_csv(OUT / 'source_data/a549_ddx3x_cpm/plotted_cpm.csv').set_index('sgRNA')
    for row in values.loc[values.cell_line.eq('A549')].itertuples():
        np.testing.assert_allclose(row.cpm, old.loc[row.sgrna, row.sample], rtol=1e-12)
    values.to_csv(DATA / 'plotted_cpm_and_raw_counts.csv', index=False)
    pd.DataFrame(totals_records).to_csv(DATA / 'library_totals.csv', index=False)
    pd.DataFrame({'sgrna': GUIDES, 'sequence': [g.rsplit('_', 1)[1] for g in GUIDES]}).to_csv(DATA / 'guide_key.csv', index=False)
    return values, hashes


def point_label(value):
    return f'{value:,.0f}' if value >= 1000 else f'{value:.3g}'


def endpoint_zoom(ax, sub):
    inset = ax.inset_axes([.47, .46, .51, .41])
    ends = sub.loc[sub.stage.eq('Endpoint')]
    for i, row in enumerate(ends.itertuples()):
        colour = COLOURS[row.genotype]
        inset.scatter(i, row.cpm, s=16, color=colour, edgecolor='white', lw=.35, clip_on=False, zorder=3)
        inset.annotate(f'{row.cpm:.3g}', (i, row.cpm), xytext=(0, 3), textcoords='offset points',
                       ha='center', fontsize=5.5, color=colour)
    inset.set(xlim=(-.5, 3.5), ylim=(0, 55))
    inset.set_yticks([0, 25, 50])
    inset.set_xticks(range(4), ['CTRL\nR1', 'CTRL\nR2', 'LOY\nR1', 'LOY\nR2'])
    inset.tick_params(axis='both', labelsize=5.5, length=2, pad=1)
    for label, row in zip(inset.get_xticklabels(), ends.itertuples()):
        label.set_color(COLOURS[row.genotype])
    inset.set_title('Endpoint zoom (linear CPM)', fontsize=6.1, pad=3)
    inset.patch.set_alpha(1)
    return inset


def main():
    values, hashes = prepare()
    fig, axes = plt.subplots(4, 2, figsize=(183 / 25.4, 210 / 25.4))
    fig.subplots_adjust(left=.105, right=.985, bottom=.095, top=.932, wspace=.25, hspace=.60)
    for col, line in enumerate(SOURCES):
        box = axes[0, col].get_position()
        fig.text((box.x0 + box.x1) / 2, .982, line, ha='center', va='top', fontsize=11, fontweight='semibold')
    for row_index, guide in enumerate(GUIDES):
        for col, line in enumerate(SOURCES):
            ax = axes[row_index, col]
            sub = values.loc[values.cell_line.eq(line) & values.sgrna.eq(guide)]
            assert sub['sample'].tolist() == SAMPLES[line]
            extreme = line == 'HT1080' and guide == GUIDES[1]
            ax.axvline(1.5 if line == 'HT1080' else 3.5, color='#B7BDC5', lw=.7, ls='--', zorder=0)
            for i, row in enumerate(sub.itertuples()):
                colour = COLOURS[row.genotype]
                ax.scatter(i, row.cpm, s=27, marker='D' if row.stage == 'D0' else 'o',
                           color=colour, edgecolor='white', lw=.5, zorder=3, clip_on=False)
                ax.annotate(point_label(row.cpm), (i, row.cpm), xytext=(0, 5),
                            textcoords='offset points', ha='center', va='bottom', fontsize=6.2, color=colour)
            ax.set_yscale('linear')
            ax.set(xlim=(-.55, len(sub) - .45), ylim=(0, 26000 if extreme else 50))
            ax.set_yticks([0, 10000, 20000] if extreme else [0, 10, 20, 30, 40, 50])
            ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f'{value:,.0f}'))
            ax.tick_params(axis='y', labelsize=6.7, length=3, pad=3)
            ax.set_xticks(range(len(sub)), LABELS[line], fontsize=6.1)
            ax.tick_params(axis='x', length=0, pad=4)
            for label, r in zip(ax.get_xticklabels(), sub.itertuples()):
                label.set_color(COLOURS[r.genotype])
            ax.set_title('DDX3X guide ' + guide.split('_')[1], fontsize=8.5, pad=7)
            ax.text(-.14, 1.07, chr(97 + 2 * row_index + col), transform=ax.transAxes,
                    fontsize=10.5, fontweight='semibold')
            ax.grid(axis='y', color='#E9EBEE', lw=.4, zorder=0)
            if sub.raw_count.sum() == 0:
                ax.text(.5, .48, 'Zero reads in all samples', transform=ax.transAxes,
                        ha='center', fontsize=7.3, color='#666666')
            if extreme:
                endpoint_zoom(ax, sub)
    fig.supylabel('Guide abundance (CPM)', x=.008, fontsize=9)
    fig.legend(handles=[
        Line2D([], [], marker='D', ls='', color='#666666', label='D0 libraries'),
        Line2D([], [], marker='o', ls='', color='#666666', label='Endpoint cultures')],
        loc='lower center', bbox_to_anchor=(.53, .026), ncol=2, fontsize=7.2, handletextpad=.4)
    fig.text(.53, .015, 'Linear scales throughout; panel c has an expanded range and endpoint zoom.',
             ha='center', fontsize=6.8, color='#666666')
    for ext in ['pdf', 'svg', 'png']:
        fig.savefig(OUT / f'{STEM}.{ext}', dpi=600)
    preview = OUT / 'review' / f'{STEM}_preview.png'
    fig.savefig(preview, dpi=180)
    hashed_preview = preview.with_name(f'{STEM}_{digest(preview)[:12]}.png')
    hashed_preview.write_bytes(preview.read_bytes())
    plt.close(fig)
    with pymupdf.open(OUT / f'{STEM}.pdf') as doc:
        page = doc[0]
        assert all('OpenSans' in f[3].replace(' ', '') for f in page.get_fonts(full=True))
        assert not [w for w in page.get_text('words') if
                    w[0] < -.5 or w[1] < -.5 or w[2] > page.rect.width + .5 or w[3] > page.rect.height + .5]
        assert 'P3' not in page.get_text()
    assert hashes == {line: digest(path) for line, path in SOURCES.items()}
    manifest = {
        'sources': {line: {'path': str(path), 'sha256': hashes[line]} for line, path in SOURCES.items()},
        'normalisation': 'Raw guide count / sum of all 77,441 guide counts in that library * 1e6; no pseudocount.',
        'n_guides': 4, 'n_unique_observations': 52, 'samples': SAMPLES,
        'a549_values_match_previous_plot': True, 'source_unchanged': True,
        'scale': 'Linear. All main panels 0–50 CPM except HT1080 guide 4638 (0–26,000); its endpoint inset is 0–55.',
        'sampling': 'HT1080: one pooled D0 library per genotype and two endpoint cultures per genotype. A549: two technical D0 libraries per genotype; one CTRL and two LOY P4 cultures.',
        'exclusions': 'All A549 P3 libraries excluded. No DDX3X guide excluded.',
        'zero_read_guide_both_screens': GUIDES[3],
        'fonts_open_sans': True, 'text_within_page': True, 'figure_mm': [183, 210], 'png_dpi': 600,
        'preview': str(hashed_preview),
    }
    (DATA / 'manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    print(json.dumps(manifest, indent=2))


if __name__ == '__main__':
    main()
