"""Render ch04_supp_depmap_model_positions from the frozen source data."""
from pathlib import Path
import hashlib
import json
import sys

THESIS = Path(__file__).resolve().parents[3]
DEPS = THESIS.parent / '.figure_python_deps'
if DEPS.is_dir():
    sys.path.insert(0, str(DEPS))
sys.path.insert(0, str(THESIS / 'skills/thesis-figure-style/scripts'))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import pymupdf
from thesis_style import apply_style
apply_style(7)

FIGS = THESIS / 'figures/chapter4'
DATA = FIGS / 'source_data/depmap_model_positions'
REVIEW = FIGS / 'review'
STEM = 'ch04_supp_depmap_model_positions'
BLUE, RED, GREY, DARK = '#3C5488', '#E64B35', '#C9CDD2', '#222222'
XLIM = (-4.2, 1.7)
EDGE = XLIM[1] - .04
COLS = (1.05, 1.22)  # percentile columns, axes fraction
LABELS = {'ICK': 'ICK/CILK1'}


def strip_panel(ax, pos, long, letter):
    genes = pos.gene.tolist()
    rng = np.random.default_rng(7)
    for i, g in enumerate(genes):
        v = long.loc[long.gene == g, 'chronos'].to_numpy()
        # Values beyond the shared axis (TP53 only) are pinned at the edge and counted.
        over = v > EDGE
        ax.scatter(np.minimum(v, EDGE), i + rng.uniform(-.2, .2, len(v)), s=.6, color=GREY, lw=0,
                   rasterized=True, zorder=1)
        if over.any():
            ax.text(XLIM[1], i - .32, f'{over.sum()} >1.7', ha='right', va='bottom', fontsize=5.3, color='.35')
        p = pos.iloc[i]
        for q1, med, q3, col, off in [(p.XY_q1, p.XY_median, p.XY_q3, BLUE, -.2),
                                      (p.LOY_q1, p.LOY_median, p.LOY_q3, RED, .2)]:
            ax.plot([q1, q3], [i + off] * 2, color=col, lw=1.3, solid_capstyle='butt', zorder=2)
            ax.plot([med], [i + off], marker='|', ms=4.5, mew=1.2, color=col, zorder=3)
        ax.scatter(p.HT1080_chronos, i, marker='D', s=13, color=DARK, lw=0, zorder=4)
        ax.scatter(p.A549_chronos, i, marker='o', s=15, facecolor='white', edgecolor=DARK, lw=.8, zorder=4)
        # Percent of Y-retaining males more dependent than each parental line.
        for x, val in [(COLS[0], p.HT1080_pct_XY), (COLS[1], p.A549_pct_XY)]:
            ax.text(x, i, f'{val:.0f}', transform=ax.get_yaxis_transform(), ha='center', va='center', fontsize=6.3)
    ax.axvline(0, color='.75', lw=.5, zorder=0)
    ax.axvline(-1, color='.75', lw=.5, ls=(0, (3, 2)), zorder=0)
    ax.set(xlim=XLIM, ylim=(len(genes) - .4, -.6), xticks=[-4, -3, -2, -1, 0, 1])
    ax.set_yticks(range(len(genes)), [LABELS.get(g, g) for g in genes])
    for lab, s in zip(ax.get_yticklabels(), pos.saga_reference):
        lab.set_fontweight('semibold' if s else 'normal')
    ax.tick_params(axis='y', length=0, pad=2, labelsize=6.5)
    ax.tick_params(axis='x', labelsize=6.5)
    ax.spines['left'].set_visible(False)
    for x, lab in [(COLS[0], 'HT1080'), (COLS[1], 'A549')]:
        ax.text(x, -.9, lab, transform=ax.get_yaxis_transform(), ha='center', va='bottom', fontsize=6.3,
                fontweight='semibold')
    if letter in 'ac':
        ax.text(sum(COLS) / 2, -1.95, '% more dependent', transform=ax.get_yaxis_transform(), ha='center',
                va='bottom', fontsize=5.8, color='.3')
    ax.text(-.30, 1.0, letter, transform=ax.transAxes, fontsize=10, fontweight='semibold', va='bottom')


def p_text(p):
    if p >= .001:
        return f'{p:.3f}'
    m, e = f'{p:.1e}'.split('e')
    return f'{m} × 10$^{{{int(e)}}}$'


def arm_panel(ax, d, ntc, stats):
    for _, row in ntc.iterrows():
        col = BLUE if row.arm == 'CTRL' else RED
        ax.axhspan(row.q05, row.q95, color=col, alpha=.10, lw=0, zorder=0)
        ax.axhline(row['median'], color=col, lw=.6, ls=(0, (3, 2)), alpha=.8, zorder=0)
    for _, r in d.iterrows():
        ax.plot([r.HT1080_chronos] * 2, [r.CTRL_lfc, r.LOY_lfc], color='.65', lw=.5, zorder=1)
        m = 's' if r.saga_reference else 'o'
        ax.scatter(r.HT1080_chronos, r.CTRL_lfc, marker=m, s=14, color=BLUE, lw=0, zorder=3)
        ax.scatter(r.HT1080_chronos, r.LOY_lfc, marker=m, s=14, color=RED, lw=0, zorder=3)
    for g, dx, dy, ha in [('TP53', -.06, .05, 'right'), ('DDX3X', .06, .02, 'left'), ('ITGB5', 0, .16, 'center'),
                          ('SREBF1', .06, -.02, 'left'), ('MED12', .06, -.02, 'left')]:
        r = d.set_index('gene').loc[g]
        y = r.LOY_lfc if g != 'TP53' else r.CTRL_lfc
        ax.annotate(g, (r.HT1080_chronos, y), xytext=(r.HT1080_chronos + dx, y + dy), fontsize=6, ha=ha, va='center')
    ax.axvline(0, color='.75', lw=.5, zorder=0)
    ax.set(xlim=(-1.45, 1.65), ylim=(-2.45, 1.25), xlabel='HT1080 parental Chronos (DepMap)',
           ylabel='Change from D0 to endpoint (log$_2$ FC)')
    ax.tick_params(labelsize=6.5)
    s1, s2 = stats['CTRL'], stats['LOY']
    ax.text(.02, .97, f"Spearman ρ, n = {s1['n_genes']} genes\n"
            f"CTRL {s1['spearman_rho']:.2f}, P = {p_text(s1['p_two_sided'])}\n"
            f"LOY {s2['spearman_rho']:.2f}, P = {p_text(s2['p_two_sided'])}",
            transform=ax.transAxes, fontsize=6.2, va='top', linespacing=1.3)
    ax.legend(handles=[Line2D([], [], marker='o', ls='', color=BLUE, ms=3.5, label='CTRL arm'),
                       Line2D([], [], marker='o', ls='', color=RED, ms=3.5, label='LOY arm'),
                       Line2D([], [], marker='s', ls='', color='.45', ms=3.5, label='SAGA member'),
                       Patch(color='.5', alpha=.2, lw=0, label='NTC 5th–95th pct.')],
              loc='lower right', frameon=False, fontsize=6.2, handletextpad=.3, labelspacing=.3, borderaxespad=.2)
    ax.text(-.30, 1.0, 'd', transform=ax.transAxes, fontsize=10, fontweight='semibold', va='bottom')


def main():
    pos = pd.read_csv(DATA / 'gene_positions.csv')
    long = pd.read_csv(DATA / 'plotted_model_values.csv.gz')
    d = pd.read_csv(DATA / 'ht1080_hit_arm_changes.csv')
    d = d[d.plotted]
    ntc = pd.read_csv(DATA / 'ht1080_ntc_arm_reference.csv')
    stats = json.loads((DATA / 'source_manifest.json').read_text(encoding='utf-8'))['spearman']
    assert long.chronos.min() > XLIM[0] and pos[['HT1080_chronos', 'A549_chronos']].lt(EDGE).all().all()
    clipped = long[long.chronos > EDGE].groupby('gene').size().to_dict()
    assert clipped == {'TP53': 47}, clipped
    a, b, c = (pos[pos.panel == k] for k in 'abc')

    W, H = 7.2, 7.15
    fig = plt.figure(figsize=(W, H))
    left, width, top = .105, .30, .92
    right_x = .60
    rows_l = len(a) + len(b)
    gap = .045
    unit = (top - .075 - gap) / rows_l
    ax_a = fig.add_axes([left, top - unit * len(a), width, unit * len(a)])
    ax_b = fig.add_axes([left, .075, width, unit * len(b)])
    ax_c = fig.add_axes([right_x, top - unit * len(c), width, unit * len(c)])
    ax_d = fig.add_axes([right_x, .075, .365, top - unit * len(c) - .075 - .105])
    strip_panel(ax_a, a, long, 'a')
    strip_panel(ax_b, b, long, 'b')
    strip_panel(ax_c, c, long, 'c')
    ax_a.tick_params(axis='x', labelbottom=False)
    for ax in (ax_b, ax_c):
        ax.set_xlabel('Chronos gene effect (DepMap 25Q3)', fontsize=7)
    for ax, lab in [(ax_a, 'X-linked gametologs'), (ax_b, 'HT1080 depleted hits'), (ax_c, 'A549 depleted hits')]:
        ax.text(0, 1.0, lab, transform=ax.transAxes, fontsize=7, color='.3', va='bottom', ha='left')
    arm_panel(ax_d, d, ntc, stats)
    fig.legend(handles=[Line2D([], [], marker='D', ls='', color=DARK, ms=3.5, label='HT1080 (parental)'),
                        Line2D([], [], marker='o', ls='', mfc='white', mec=DARK, ms=4, label='A549 (parental)'),
                        Line2D([], [], color=BLUE, lw=1.3, label='Y-retaining males, median and IQR (n = 324)'),
                        Line2D([], [], color=RED, lw=1.3, label='Y-depleted males (n = 229)'),
                        Line2D([], [], marker='o', ls='', color=GREY, ms=2.5, label='All models (n = 1,186)')],
               loc='upper left', bbox_to_anchor=(left - .09, .998), ncol=5, frameon=False, fontsize=6.2,
               handletextpad=.3, columnspacing=1.1)

    for ext in ['pdf', 'svg', 'png']:
        fig.savefig(FIGS / f'{STEM}.{ext}', dpi=600)
    preview = REVIEW / f'{STEM}_preview.png'
    fig.savefig(preview, dpi=170)
    plt.close(fig)
    digest = hashlib.sha256(preview.read_bytes()).hexdigest()[:12]
    hashed = REVIEW / f'{STEM}_{digest}.png'
    hashed.write_bytes(preview.read_bytes())
    with pymupdf.open(FIGS / f'{STEM}.pdf') as pdf:
        page = pdf[0]
        fonts = [f[3] for f in page.get_fonts()]
        assert all('OpenSans' in f.replace(' ', '') for f in fonts), fonts
        outside = [w[4] for w in page.get_text('words')
                   if w[0] < 0 or w[1] < 0 or w[2] > page.rect.width or w[3] > page.rect.height]
        assert not outside, outside
    val = dict(stem=STEM, size_inches=[W, H], rows={'a': len(a), 'b': len(b), 'c': len(c)}, panel_d_genes=len(d),
               fonts=fonts, pdf_text_within_page=True, preview=str(hashed),
               clipped_at_axis_edge=clipped, x_limits=list(XLIM))
    (DATA / 'render_validation.json').write_text(json.dumps(val, indent=2), encoding='utf-8')
    print(json.dumps(val, indent=2))


if __name__ == '__main__':
    main()
