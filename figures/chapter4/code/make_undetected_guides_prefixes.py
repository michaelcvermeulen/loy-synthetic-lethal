"""Compact exploratory supplement: shared undetected guides and 5-prime prefixes."""
from pathlib import Path
import hashlib
import json
import os
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT.parent / '.figure_python_deps'))
OUT = ROOT / 'figures/chapter4'
DATA = OUT / 'source_data/undetected_guides_prefixes'
DATA.mkdir(parents=True, exist_ok=True)
os.environ.setdefault('MPLCONFIGDIR', str(DATA / 'mpl_cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Patch
import numpy as np
import pandas as pd
import pymupdf

sys.path.insert(0, str(ROOT / 'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
apply_style(8)
AUDIT = ROOT / 'research/shared_undetected_guides_20261020'
STEM = 'ch04_supp_undetected_guides_prefixes_draft'
COLOURS = {'TCTC': '#00A087', 'AGACG': '#8491B4', 'Other': '#E9EBEE'}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prepare():
    audit = json.loads((AUDIT / 'summary.json').read_text(encoding='utf-8'))
    raw = {}
    for line, source in audit['sources'].items():
        path = Path(source['path'])
        assert digest(path) == source['sha256']
        raw[line] = pd.read_csv(path, sep='\t').set_index('sgRNA').sort_index()
    ht, a549 = raw['HT1080'], raw['A549']
    pd.testing.assert_index_equal(ht.index, a549.index)
    pd.testing.assert_series_equal(ht.Gene, a549.Gene)
    assert len(ht) == 77441 and ht.index.is_unique
    flags = {}
    for line, frame in raw.items():
        columns = audit['samples'][line]
        assert not any('_P3_' in c for c in columns)
        flags[line] = frame[columns].eq(0).all(axis=1)
    sequences = pd.Series(ht.index.str.rsplit('_', n=1).str[-1], index=ht.index)
    masks = {
        'HT1080 only': flags['HT1080'] & ~flags['A549'],
        'Both screens': flags['HT1080'] & flags['A549'],
        'A549 only': ~flags['HT1080'] & flags['A549'],
        'Detected in each': ~flags['HT1080'] & ~flags['A549'],
    }
    assert np.array_equal(sum(mask.astype(int) for mask in masks.values()), np.ones(len(ht)))
    expected = pd.read_csv(AUDIT / 'guides_undetected_in_either_screen.csv')
    assert set(ht.index[flags['HT1080'] | flags['A549']]) == set(expected.sgRNA)
    prefix_a, prefix_b = sequences.str.startswith('TCTC'), sequences.str.startswith('AGACG')
    assert not (prefix_a & prefix_b).any()
    assert sequences.loc['DDX3X_4640_TCTCGGTTCCTTAAATGAGG'].startswith('TCTC')
    assert masks['Both screens'].loc['DDX3X_4640_TCTCGGTTCCTTAAATGAGG']
    records = []
    for name, mask in masks.items():
        n = int(mask.sum())
        first, second = int((mask & prefix_a).sum()), int((mask & prefix_b).sum())
        records.append(dict(group=name, n_guides=n, TCTC=first, AGACG=second,
                            Other=n-first-second, either_prefix=first+second,
                            TCTC_percent=100*first/n, AGACG_percent=100*second/n,
                            Other_percent=100*(n-first-second)/n,
                            either_prefix_percent=100*(first+second)/n))
    table = pd.DataFrame(records).set_index('group')
    assert table.n_guides.tolist() == [13, 223, 84, 77121]
    assert table.loc['Both screens', ['TCTC', 'AGACG']].tolist() == [163, 33]
    assert table.loc['Detected in each', ['TCTC', 'AGACG']].tolist() == [180, 41]
    assert np.allclose(table[['TCTC_percent', 'AGACG_percent', 'Other_percent']].sum(axis=1), 100)
    table.to_csv(DATA / 'prefix_composition.csv')
    shared = table.loc['Both screens']
    background = table.loc['Detected in each']
    effects = []
    for prefix in ['TCTC', 'AGACG', 'either_prefix']:
        freq_shared = shared[prefix] / shared.n_guides
        freq_background = background[prefix] / background.n_guides
        effects.append(dict(prefix=prefix, shared_count=int(shared[prefix]), shared_total=int(shared.n_guides),
                            comparison_count=int(background[prefix]), comparison_total=int(background.n_guides),
                            shared_frequency=float(freq_shared), comparison_frequency=float(freq_background),
                            frequency_ratio=float(freq_shared/freq_background)))
    effects = pd.DataFrame(effects).set_index('prefix')
    effects.to_csv(DATA / 'descriptive_enrichment.csv')
    membership = ht[['Gene']].copy()
    membership['sequence'] = sequences
    for name, mask in masks.items():
        membership.loc[mask, 'group'] = name
    membership['prefix_class'] = np.select([prefix_a, prefix_b], ['TCTC', 'AGACG'], default='Other')
    membership.loc[flags['HT1080'] | flags['A549']].to_csv(DATA / 'undetected_guide_membership.csv')
    return table, effects, audit


def main():
    table, effects, audit = prepare()
    fig = plt.figure(figsize=(183/25.4, 76/25.4))
    fig.text(.018, .955, 'a', fontweight='semibold', fontsize=11, va='top')
    fig.text(.413, .955, 'b', fontweight='semibold', fontsize=11, va='top')
    fig.text(.205, .945, 'Undetected-guide overlap', ha='center', va='top', fontsize=9)
    fig.text(.713, .945, '5′ sequence prefixes', ha='center', va='top', fontsize=9)

    # A schematic Venn makes the very small exclusive subset legible; counts,
    # rather than circle areas, encode set size.
    venn = fig.add_axes([.023, .185, .36, .66])
    venn.set(xlim=(0, 3.15), ylim=(0, 2.43), aspect='equal')
    venn.axis('off')
    for centre in [(1.05, 1.14), (2.08, 1.14)]:
        venn.add_patch(Circle(centre, .91, fc='#E6E8EC', ec='#69717B', alpha=.65, lw=.8))
    venn.text(.83, 2.24, 'HT1080\n236', ha='center', va='center', fontsize=8, linespacing=1.25)
    venn.text(2.28, 2.24, 'A549\n307', ha='center', va='center', fontsize=8, linespacing=1.25)
    venn.text(.56, 1.19, '13', ha='center', va='center', fontsize=13)
    venn.text(1.565, 1.19, '223', ha='center', va='center', fontsize=13, fontweight='semibold')
    venn.text(2.57, 1.19, '84', ha='center', va='center', fontsize=13)
    venn.annotate('Includes DDX3X 4640', xy=(1.565, .65), xytext=(1.565, .11),
                  ha='center', va='center', fontsize=7,
                  arrowprops=dict(arrowstyle='-', color='#777777', lw=.6))
    fig.text(.203, .105, 'Zero reads at D0 and endpoint', ha='center', fontsize=7, color='#555555')

    ax = fig.add_axes([.575, .23, .28, .52])
    groups = table.index.tolist()
    y = np.arange(len(groups))[::-1]
    left = np.zeros(len(groups))
    for prefix in ['TCTC', 'AGACG', 'Other']:
        width = table[f'{prefix}_percent'].to_numpy()
        ax.barh(y, width, left=left, height=.56, color=COLOURS[prefix], edgecolor='none', zorder=2)
        left += width
    ax.set(xlim=(0, 100), ylim=(-.6, 3.6))
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xlabel('Guides (%)', fontsize=7.5, labelpad=3)
    ax.set_yticks(y, groups, fontsize=7.2)
    ax.tick_params(axis='y', length=0, pad=6)
    ax.tick_params(axis='x', length=3, labelsize=7)
    ax.spines['left'].set_visible(False)
    ax.grid(axis='x', color='#E9EBEE', lw=.4, zorder=0)
    ax.axhline(.5, color='#B8BDC4', lw=.6, ls='--', xmin=-.04, xmax=1.03, clip_on=False)
    for yi, (_, row) in zip(y, table.iterrows()):
        percent = f'{row.either_prefix_percent:.3f}%' if row.either_prefix_percent < 1 else f'{row.either_prefix_percent:.1f}%'
        ax.text(104, yi, f'{int(row.either_prefix)}/{int(row.n_guides):,}\n{percent}',
                va='center', fontsize=6.7, linespacing=1.2, clip_on=False)
    ax.text(104, 3.83, 'Either prefix', fontsize=6.5, clip_on=False)
    fig.legend(handles=[Patch(fc=COLOURS[p], label=p if p != 'Other' else 'Other') for p in COLOURS],
               loc='upper center', bbox_to_anchor=(.715, .866), ncol=3, frameon=False,
               fontsize=7.2, handlelength=1.2, handletextpad=.4, columnspacing=1.1)
    ratio = effects.loc['either_prefix', 'frequency_ratio']
    fig.text(.715, .058, f'Shared set: {ratio:.0f}× prefix enrichment', ha='center', fontsize=7.5, fontweight='semibold')
    for ext in ['pdf', 'svg', 'png']:
        fig.savefig(OUT / f'{STEM}.{ext}', dpi=600)
    preview = OUT / 'review' / f'{STEM}_preview.png'
    fig.savefig(preview, dpi=200)
    hashed = preview.with_name(f'{STEM}_{digest(preview)[:12]}.png')
    hashed.write_bytes(preview.read_bytes())
    plt.close(fig)
    with pymupdf.open(OUT / f'{STEM}.pdf') as doc:
        page = doc[0]
        assert all('OpenSans' in f[3].replace(' ', '') for f in page.get_fonts(full=True))
        assert not [w for w in page.get_text('words') if
                    w[0] < -.5 or w[1] < -.5 or w[2] > page.rect.width + .5 or w[3] > page.rect.height + .5]
    for source in audit['sources'].values():
        assert digest(Path(source['path'])) == source['sha256']
    manifest = {
        'sources': audit['sources'], 'samples': audit['samples'],
        'audit_summary': str(AUDIT / 'summary.json'), 'audit_summary_sha256': digest(AUDIT / 'summary.json'),
        'definition': 'Undetected means zero raw counts throughout every retained D0 and endpoint library in a screen.',
        'overlap': {'HT1080_only': 13, 'shared': 223, 'A549_only': 84, 'neither': 77121},
        'venn_areas': 'Schematic, not proportional to counts.',
        'comparison_group': '77,121 guides with at least one read somewhere in each screen; excludes all 320 guides undetected throughout either screen.',
        'prefix_definition': 'Exact TCTC or AGACG match at the 5-prime start of the 20-nt reference sequence; mutually exclusive.',
        'frequency_ratio_either_prefix_shared_vs_comparison': float(ratio),
        'statistics': 'Descriptive frequency ratio; post hoc prefix selection, no P values or causal claims.',
        'exclusions': 'A549 P3 libraries remain excluded. No guide-level count or candidate analysis modified.',
        'figure_mm': [183, 76], 'png_dpi': 600, 'fonts_open_sans': True, 'text_within_page': True,
        'preview': str(hashed), 'source_hashes_unchanged': True,
    }
    (DATA / 'manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    print(table.to_string())
    print(effects.to_string())
    print('PREVIEW:', hashed)


if __name__ == '__main__':
    main()
