"""Combined Methods supplement (22 Oct 2026): shared guide non-detection (a,b) and
shared jackpotting guides with the motif scan (c-e).

Merges ch04_supp_undetected_guides_prefixes_draft and
ch04_supp_shared_jackpot_guides_sequence_motifs without changing any value: panels
a,b reuse prepare() from make_undetected_guides_prefixes.py, and c-e read the saved
outputs of research/jackpot_sequence_motifs_20261021/analyse.py.
"""
from pathlib import Path
import hashlib
import json
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import make_undetected_guides_prefixes as U  # sets style, paths and matplotlib backend

plt, np, pd, pymupdf = U.plt, U.np, U.pd, U.pymupdf
from matplotlib.patches import Circle, Ellipse, Patch

OUT = U.OUT
JDATA = OUT / 'source_data/jackpot_sequence_motifs_20261021'
DATA = OUT / 'source_data/guide_artefacts_combined_20261022'
DATA.mkdir(parents=True, exist_ok=True)
STEM = 'ch04_supp_guide_artefacts_combined'
W, H_TOP, H_BOTTOM = 9.6, 3.0, 5.45
LETTER = dict(fontsize=13, fontweight='semibold', va='top')


def draw_undetected(sf, table, effects):
    """Panels a,b on a 9.6 x 3.0 in subfigure (layout of the former standalone figure)."""
    sf.text(.018, .955, 'a', **LETTER)
    sf.text(.413, .955, 'b', **LETTER)
    sf.text(.205, .945, 'Undetected-guide overlap', ha='center', va='top', fontsize=9)
    sf.text(.713, .945, '5′ sequence prefixes', ha='center', va='top', fontsize=9)
    venn = sf.add_axes([.023, .185, .36, .66])
    venn.set(xlim=(0, 3.15), ylim=(0, 2.43), aspect='equal')
    venn.axis('off')
    for centre in [(1.05, 1.14), (2.08, 1.14)]:
        venn.add_patch(Circle(centre, .91, fc='#E6E8EC', ec='#69717B', alpha=.65, lw=.8))
    venn.text(.83, 2.24, 'HT1080\n236', ha='center', va='center', fontsize=8, linespacing=1.25)
    venn.text(2.28, 2.24, 'A549\n307', ha='center', va='center', fontsize=8, linespacing=1.25)
    venn.text(.56, 1.19, '13', ha='center', va='center', fontsize=13)
    venn.text(1.565, 1.19, '223', ha='center', va='center', fontsize=13, fontweight='semibold')
    venn.text(2.57, 1.19, '84', ha='center', va='center', fontsize=13)
    venn.annotate('Includes DDX3X 4640', xy=(1.565, .65), xytext=(1.565, .11), ha='center', va='center',
                  fontsize=7, arrowprops=dict(arrowstyle='-', color='#777777', lw=.6))
    sf.text(.203, .105, 'Zero reads at D0 and endpoint', ha='center', fontsize=7, color='#555555')
    ax = sf.add_axes([.575, .23, .28, .52])
    groups = table.index.tolist()
    y = np.arange(len(groups))[::-1]
    left = np.zeros(len(groups))
    for prefix in ['TCTC', 'AGACG', 'Other']:
        width = table[f'{prefix}_percent'].to_numpy()
        ax.barh(y, width, left=left, height=.56, color=U.COLOURS[prefix], edgecolor='none', zorder=2)
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
        pct = f'{row.either_prefix_percent:.3f}%' if row.either_prefix_percent < 1 else f'{row.either_prefix_percent:.1f}%'
        ax.text(104, yi, f'{int(row.either_prefix)}/{int(row.n_guides):,}\n{pct}', va='center', fontsize=6.7,
                linespacing=1.2, clip_on=False)
    ax.text(104, 3.83, 'Either prefix', fontsize=6.5, clip_on=False)
    sf.legend(handles=[Patch(fc=U.COLOURS[p], label=p) for p in U.COLOURS], loc='upper center',
              bbox_to_anchor=(.715, .866), ncol=3, frameon=False, fontsize=7.2, handlelength=1.2,
              handletextpad=.4, columnspacing=1.1)
    ratio = effects.loc['either_prefix', 'frequency_ratio']
    sf.text(.715, .058, f'Shared set: {ratio:.0f}× prefix enrichment', ha='center', fontsize=7.5, fontweight='semibold')
    return ratio


def draw_jackpot(sf):
    """Panels c-e on a 9.6 x 5.45 in subfigure (layout of the former standalone figure)."""
    audit = pd.read_csv(JDATA / 'guide_audit.csv').set_index('sgRNA')
    m = json.loads((JDATA / 'manifest.json').read_text())
    summary = pd.read_csv(JDATA / 'motif_scan_summary.csv')
    shared = audit[audit.HT1080_mask173 & audit.A549_mask342]
    assert len(shared) == 172
    a = sf.add_axes([.055, .49, .37, .44]); a.set(xlim=(0, 1), ylim=(0, 1)); a.axis('off')
    b = sf.add_axes([.57, .54, .39, .40])
    c = sf.add_axes([.075, .13, .885, .25]); c.axis('off')
    sf.text(.018, .975, 'c', **LETTER)
    sf.text(.505, .975, 'd', **LETTER)
    sf.text(.018, .405, 'e', **LETTER)
    for centre in [.35, .65]:
        a.add_patch(Ellipse((centre, .52), .60, .58, fc='#E4E0EA', ec='#6A6176', alpha=.75, lw=1))
    a.text(.20, .96, 'HT1080 D0\n173 targeting guides', ha='center', va='top', fontsize=8.5)
    a.text(.80, .96, 'A549 P3\n272 targeting guides', ha='center', va='top', fontsize=8.5)
    for x, label in [(.20, '1'), (.50, '172'), (.80, '100')]:
        a.text(x, .52, label, ha='center', va='center', fontsize=17, fontweight='semibold',
               color='#5E3C73' if x == .5 else '#333333')
    a.text(.50, .12, '172 of 173 HT1080 outliers shared (99.4%)', ha='center', fontsize=8)
    a.text(.50, .025, 'A549 exclusion also contains 70 NTCs, analysed separately.', ha='center', fontsize=7, color='#555555')
    x = np.log10(shared.HT1080_D0_mean_CPM + 1); y = np.log10(shared.A549_P3_mean_CPM + 1)
    b.scatter(x, y, s=15, color='#6D4C7D', alpha=.72, edgecolor='white', lw=.3, rasterized=True)
    b.plot([0, 4.7], [0, 4.7], ls=':', color='#AAAAAA', lw=.8, zorder=0)
    b.set(xlim=(0, 4.7), ylim=(0, 4.7), xlabel='HT1080 D0 log₁₀(mean CPM + 1)', ylabel='A549 P3 log₁₀(mean CPM + 1)')
    b.text(.035, .965, f"172 shared guides\nSpearman ρ = {m['shared_abundance_spearman']:.3f}\n70 genes; 30 with all four guides",
           transform=b.transAxes, va='top', fontsize=8, linespacing=1.4)
    rows = []
    for line in ['HT1080_D0', 'A549_P3']:
        for feature in ['prefix', 'suffix', 'anywhere']:
            r = summary[(summary['set'] == line) & (summary.feature == feature)].iloc[0]
            rows.append([r['set'].replace('_', ' '), {'prefix': '5′ prefix', 'suffix': '3′ suffix', 'anywhere': 'Anywhere'}[r.feature],
                         r.best_motif, f'{int(r.best_foreground_with)}/{int(r.best_foreground_n)}',
                         f'{int(r.best_background_with)}/{int(r.best_background_n):,}', f'{r.best_fold:.1f}×', f'{r.best_BH:.3f}'])
    table = c.table(cellText=rows, colLabels=['Guide set', 'Feature', 'Top nominal motif', 'Flagged', 'Remaining', 'Ratio', 'BH FDR'],
                    colWidths=[.16, .12, .21, .12, .16, .105, .125], bbox=[0, 0, 1, 1], cellLoc='center')
    table.auto_set_font_size(False); table.set_fontsize(7.5)
    for (i, j), cell in table.get_celld().items():
        cell.set_edgecolor('white'); cell.set_facecolor('#EAEAEA' if i == 0 else '#F5F5F5' if i % 2 else 'white')
        if j == 0:
            cell.set_text_props(ha='left')
    sf.text(.075, .095, 'No tested motif passed FDR <0.05. The table shows the lowest nominal P value in each feature class.', fontsize=7.6)
    sf.text(.075, .055, 'Targeting guides only; one-sided exact enrichment tests against remaining targeting guides; BH across all 32,880 primary tests.', fontsize=7, color='#555555')
    sf.text(.075, .020, 'Guide identity and abundance recur across affected libraries; this analysis does not establish the source of the distortion.', fontsize=7, color='#555555')
    return rows


def main():
    table, effects, audit = U.prepare()
    fig = plt.figure(figsize=(W, H_TOP + H_BOTTOM))
    top, bottom = fig.subfigures(2, 1, height_ratios=[H_TOP, H_BOTTOM])
    ratio = draw_undetected(top, table, effects)
    rows = draw_jackpot(bottom)
    for ext in ['pdf', 'svg', 'png']:
        fig.savefig(OUT / f'{STEM}.{ext}', dpi=600)
    preview = OUT / 'review' / f'{STEM}_preview.png'
    fig.savefig(preview, dpi=170)
    plt.close(fig)
    hashed = preview.with_name(f'{STEM}_{hashlib.sha256(preview.read_bytes()).hexdigest()[:12]}.png')
    preview.replace(hashed)
    with pymupdf.open(OUT / f'{STEM}.pdf') as doc:
        page = doc[0]
        assert all('OpenSans' in f[3].replace(' ', '') for f in page.get_fonts(full=True))
        assert not [w for w in page.get_text('words') if
                    w[0] < -.5 or w[1] < -.5 or w[2] > page.rect.width + .5 or w[3] > page.rect.height + .5]
    (DATA / 'manifest.json').write_text(json.dumps({
        'stem': STEM, 'preview': str(hashed), 'size_inches': [W, H_TOP + H_BOTTOM],
        'panels': {'a,b': 'former ch04_supp_undetected_guides_prefixes_draft a,b',
                   'c-e': 'former ch04_supp_shared_jackpot_guides_sequence_motifs a-c'},
        'either_prefix_frequency_ratio': float(ratio), 'motif_table_rows': rows,
        'values_changed': False, 'fonts_open_sans': True, 'text_within_page': True}, indent=2), encoding='utf-8')
    print('PREVIEW:', hashed)


if __name__ == '__main__':
    main()
