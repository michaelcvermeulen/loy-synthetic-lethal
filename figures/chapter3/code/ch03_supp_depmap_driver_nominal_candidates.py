"""Two dependency panels with four driver-by-Y groups; no models refitted."""
from pathlib import Path
import hashlib
import json
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

ROOT = Path(__file__).resolve().parents[3]
FIG = ROOT / 'figures/chapter3'
DATA = FIG / 'source_data/depmap_driver_nominal_candidates'
STEM = 'ch03_supp_depmap_driver_nominal_candidates'
sys.path.insert(0, str(ROOT / 'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
font = apply_style(8.5)
points = pd.read_csv(DATA / 'plotted_points.csv')
saved = pd.read_csv(DATA / 'saved_tests.csv').set_index('symbol')
assert len(points) == 106 and set(saved.index) == {'INSM1', 'ATRX'}
COLOURS = {'Y-positive': '#3C5488', 'LOY': '#E64B35'}
fig = plt.figure(figsize=(9.8, 3.75))
rng = np.random.default_rng(20261008)


def p_text(value):
    if value >= .001:
        return f'{value:.3f}'
    exponent = int(np.floor(np.log10(value)))
    return f'{value / 10 ** exponent:.2f} × 10$^{{{exponent}}}$'


def boxes(ax, part, groups):
    labels = []
    for x, context, status in groups:
        colour = COLOURS[status]
        vals = part.loc[part.ctx.eq(context) & part.y_status.eq(status), 'chronos'].to_numpy()
        labels.append(f'{status}\nn = {len(vals)}')
        if len(vals) == 0:
            ax.text(x, .39, 'No models', transform=ax.get_xaxis_transform(),
                    ha='center', va='center', fontsize=8, color='#777777')
            continue
        ax.boxplot([vals], positions=[x], widths=.48, manage_ticks=False,
                   showfliers=False, patch_artist=True,
                   boxprops={'facecolor': colour, 'edgecolor': colour, 'alpha': .20, 'lw': .8},
                   medianprops={'color': colour, 'lw': 1.5},
                   whiskerprops={'color': colour, 'lw': .8}, capprops={'color': colour, 'lw': .8})
        ax.scatter(x + rng.uniform(-.13, .13, len(vals)), vals, s=21, c=colour,
                   alpha=.85, edgecolors='white', linewidths=.4, zorder=3)
    low, high = min(part.chronos.min(), 0), max(part.chronos.max(), 0)
    spread = high - low
    ax.set_ylim(low - .08 * spread, high + .49 * spread)
    ax.axhline(0, color='#999999', lw=.65, ls=(0, (3, 3)), zorder=0)
    ax.set_xticks([group[0] for group in groups], labels)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=4))
    ax.tick_params(length=3, pad=3, labelsize=8)


def comparison(ax, left, right, effect, p, fdr=None):
    centre = (left + right) / 2
    trans = ax.get_xaxis_transform()
    ax.plot([left, left, right, right], [.745, .765, .765, .745], transform=trans,
            color='#222222', lw=.8, solid_capstyle='butt')
    ax.text(centre, .985, f'Δ = {effect:+.3f}', transform=trans, ha='center', va='top', fontsize=8.2)
    ax.text(centre, .905, f'P = {p_text(p)}', transform=trans, ha='center', va='top', fontsize=8.2)
    if fdr is not None:
        ax.text(centre, .825, f'FDR = {fdr:.3f}', transform=trans, ha='center', va='top', fontsize=8.2)


for i, gene in enumerate(['INSM1', 'ATRX']):
    ax = fig.add_axes([.075, .255, .397, .665] if i == 0 else [.587, .255, .397, .665])
    part = points[points.symbol.eq(gene)]
    hit = saved.loc[gene]
    altered = part[part.ctx.eq(1)]
    assert sum(altered.y_status.eq('Y-positive')) == hit.n_pos_in
    assert sum(altered.y_status.eq('LOY')) == hit.n_neg_in
    groups = [(0, 1, 'Y-positive'), (1, 1, 'LOY'), (2.6, 0, 'Y-positive'), (3.6, 0, 'LOY')]
    assert sum(part.ctx.eq(0) & part.y_status.eq('Y-positive')) == hit.n_pos_out
    assert sum(part.ctx.eq(0) & part.y_status.eq('LOY')) == hit.n_neg_out
    boxes(ax, part, groups)
    comparison(ax, 0, 1, hit.inside, hit.inside_p, hit.inside_fdr)
    ax.set_ylabel(f'{gene} gene effect (Chronos)')
    ax.set_xlim(-.58, 4.18)
    driver = 'MYC' if gene == 'INSM1' else 'CDKN2A'
    for centre, label in [(.5, f'{driver}-altered'), (3.1, f'{driver}-unaltered')]:
        ax.text(centre, -.205, label, transform=ax.get_xaxis_transform(), ha='center', va='top', fontsize=8.5)
    ax.text(.5, -.325, 'Lung cancer' if gene == 'INSM1' else 'Melanoma',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.3)
    if gene == 'INSM1':
        comparison(ax, 2.6, 3.6, hit.outside, hit.outside_p)
    else:
        assert hit.n_pos_out == 9 and hit.n_neg_out == 0
        assert pd.isna(hit.outside_p) and pd.isna(hit.interaction_p)
    fig.text(.012 if i == 0 else .524, .99, 'ab'[i], ha='left', va='top', fontsize=11, weight='semibold')
for extension in ['pdf', 'svg', 'png']:
    fig.savefig(FIG / f'{STEM}.{extension}', dpi=300, facecolor='white', bbox_inches='tight', pad_inches=.035)
plt.close(fig)

derived = ROOT / 'research/depmap_loy_reanalysis_v1/loy_reanalysis_v1_export_2026-10-07/derived/loy_reanalysis_v1'
sources = [derived / name for name in ['final/driver_type_screen.csv.gz',
           'final/driver_type_split_lines.csv', 'cohort.csv', 'matrices/crispr_effect.rds']]
provenance = {
    'working_id': 'CH3-S3.3.1-S07', 'font': font,
    'selection': 'The two non-core tests among 14 distinct gene-by-altered-cell-set tests with primary P < 1e-5.',
    'panels': {'a': 'INSM1 in lung: MYC-altered, 7 Y-positive and 13 LOY; MYC-unaltered, 21 Y-positive and 31 LOY',
               'b': 'ATRX in melanoma: CDKN2A-altered, 20 Y-positive and 5 LOY; CDKN2A-unaltered, 9 Y-positive and 0 LOY'},
    'display': 'Untransformed Chronos scores: four driver-by-Y groups for both genes; the empty CDKN2A-unaltered LOY group is labelled n = 0. Only within-group comparisons are annotated; interaction annotations removed at author request on 8 October 2026.',
    'effect': 'Saved adjusted raw-scale LOY-minus-Y-positive effects within altered or unaltered models, as labelled.',
    'inference': 'Saved limma moderated tests on globally rank-normalised scores, fitted jointly with available unaltered models and library covariates.',
    'fdr_family': 'Displayed within-altered FDRs: 414052 tests across 68 comparisons. The MYC-unaltered P value is nominal, without a new FDR calculation. Saved interaction results remain in source data but are not displayed.',
    'statistics_refitted': False, 'jitter_seed': 20261008,
    'input_sha256': {p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
    'plotted_points_sha256': hashlib.sha256((DATA / 'plotted_points.csv').read_bytes()).hexdigest(),
    'saved_tests_sha256': hashlib.sha256((DATA / 'saved_tests.csv').read_bytes()).hexdigest(),
    'exports': {ext: hashlib.sha256((FIG / f'{STEM}.{ext}').read_bytes()).hexdigest() for ext in ['pdf','svg','png']},
}
(DATA / 'provenance.json').write_text(json.dumps(provenance, indent=2) + '\n', encoding='utf-8')
(DATA / 'README.md').write_text(
    '# INSM1 and ATRX nominal associations in cancer-driver subsets\n\n'
    'Selected as CH3-S3.3.1-S07. These are the two non-core tests with P < 1e-5 among the '
    '14 distinct tests cited in Results. Both saved effects are positive, indicating weaker '
    'dependency in LOY; neither passed the full-search FDR threshold.\n\n'
    'The 17 original rows reduce to 14 after deduplicating dependency gene and altered-cell-set ID. '
    'Twelve involve DDX3X, EIF1AX, RPS4X or ZFX. No statistical models were refitted. '
    'Group membership and counts are checked against the saved screen. Panel a now shows all '
    'four MYC-by-Y groups in lung: altered 7 Y-positive/13 LOY and unaltered 21/31. Panel b '
    'shows all four CDKN2A-by-Y groups in melanoma: altered 20 Y-positive/5 LOY and '
    'unaltered 9 Y-positive/0 LOY. The empty group is labelled n = 0; neither a within-unaltered '
    'comparison nor a driver interaction is estimable. plotted_points.csv and '
    'all_context_points.csv contain the same 106 observations.\n\n'
    'The saved MYC-unaltered effect is +0.024830 Chronos units with nominal P = 0.788013. '
    'Within-altered FDRs retain their original 414,052-test family. Interaction annotations '
    'were removed from the figure and caption at the author request on 8 October 2026; '
    'saved interaction statistics remain in the source data. No new correction or statistical fit was performed.\n\n'
    '[Central caption](../../../../Chapter3_figure_captions.md#supp--nominal-dependency-associations-in-driver-altered-cancer-subsets)\n\n'
    '[Extraction script](../../code/ch03_supp_depmap_driver_nominal_extract.R) (Rscript --vanilla), then '
    '[plotting script](../../code/ch03_supp_depmap_driver_nominal_candidates.py).\n', encoding='utf-8')
print(FIG / f'{STEM}.png')
