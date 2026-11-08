"""Selected six-panel supplementary figure from frozen tables; no models are refitted."""
from pathlib import Path
import hashlib
import json
import shutil
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator

ROOT = Path(__file__).resolve().parents[3]
FIG = ROOT / 'figures/chapter3'
DATA = FIG / 'source_data/depmap_combined_checks'
EXPORT = ROOT / 'research/depmap_loy_reanalysis_v1/loy_reanalysis_v1_export_2026-10-07/results/loy_reanalysis_v1/figure_source_data'
STEM = 'ch03_supp_depmap_combined_checks'
DATA.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(ROOT / 'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
font = apply_style(8)
plt.rcParams.update({'xtick.major.pad': 2, 'ytick.major.pad': 3})
BLUE, RED, GREY, BLACK, TEAL, PURPLE = '#3C5488', '#E64B35', '#B8B8B8', '#222222', '#00A087', '#8172B3'

inputs = {
    'fig1c_call_tiers.csv': FIG / 'source_data/depmap_y_status_dependencies/fig1c_call_tiers.csv',
    'fig3a_hit_robustness.csv': EXPORT / 'fig3a_hit_robustness.csv',
}
for name in ('fig8a_nested_accuracy.csv', 'fig8c_y_status_forest.csv'):
    inputs[name] = FIG / 'source_data/depmap_random_forests' / name
for name in ('fig5d_arhgef2_cin_points.csv',
             'fig5d_arhgef2_discovery.csv', 'fig5e_cib1_points.csv', 'fig5e_cib1_discovery.csv'):
    inputs[name] = FIG / 'source_data/depmap_context_screens' / name
hashes = {}
for name, source in inputs.items():
    shutil.copy2(source, DATA / name)
    hashes[name] = {'source': source.relative_to(ROOT).as_posix(),
                    'sha256': hashlib.sha256(source.read_bytes()).hexdigest()}
tables = {name: pd.read_csv(DATA / name) for name in inputs}

WIDTH, HEIGHT = 8.2, 7.65
fig = plt.figure(figsize=(WIDTH, HEIGHT))


def axes(x, y, w, h):
    ax = fig.add_axes([x / WIDTH, y / HEIGHT, w / WIDTH, h / HEIGHT])
    ax.tick_params(length=2.6, labelsize=7.4)
    return ax


def label(x, y, text):
    fig.text(x / WIDTH, y / HEIGHT, text, va='top', ha='left',
             fontsize=11, weight='semibold')


def text(x, y, value, **kwargs):
    return fig.text(x / WIDTH, y / HEIGHT, value, **kwargs)


# a: All male models, retaining the exact eight evidence categories.
tiers = tables['fig1c_call_tiers.csv']
assert tiers.N.sum() == 1143
assert tiers.groupby('status').N.sum().to_dict() == {'Y-negative': 383, 'Y-positive': 588, 'set aside': 172}
tier_order = ['RNA and DNA agree', 'DNA only', 'RNA only', 'RNA, DNA partial',
              'RNA intermediate', 'DNA partial, no RNA', 'no data', 'RNA and DNA disagree']
tier_labels = ['RNA + DNA agree', 'DNA only', 'RNA only', 'RNA; partial DNA',
               'Intermediate RNA', 'Partial DNA; no RNA', 'No data', 'RNA + DNA disagree']
a = axes(1.37, 5.37, 1.14, 1.96)
left = np.zeros(len(tier_order))
ys = np.arange(len(tier_order))
for status, colour in [('Y-positive', BLUE), ('Y-negative', RED), ('set aside', GREY)]:
    values = tiers[tiers.status.eq(status)].groupby('tier').N.sum().reindex(tier_order, fill_value=0).to_numpy()
    a.barh(ys, values, left=left, height=.72, color=colour, edgecolor='white', linewidth=.3)
    for y, start, value in zip(ys, left, values):
        if value >= 40:
            # Place the small red DNA-only count outside its narrow segment.
            outside = status == 'Y-negative' and value < 70
            a.text(start + value + 15 if outside else start + value / 2, y, str(value),
                   ha='left' if outside else 'center', va='center',
                   color=colour if outside else BLACK if status == 'set aside' else 'white', fontsize=6.4)
    left += values
a.set(yticks=ys, yticklabels=tier_labels, ylim=(7.6, -.6), xlim=(0, 800),
      xticks=[0, 400, 800], xlabel='Male models')
a.tick_params(axis='y', length=0, labelsize=7)
a.spines['left'].set_visible(False)
fig.legend(handles=[Patch(facecolor=c, label=s) for s, c in
                    [('Y-positive', BLUE), ('LOY', RED), ('Set aside', GREY)]],
           loc='lower left', bbox_to_anchor=(.10 / WIDTH, 7.39 / HEIGHT), ncol=3,
           fontsize=7.1, handlelength=.9, handletextpad=.4, columnspacing=.75, borderaxespad=0)
label(.04, 7.63, 'a')

# b: Former robustness panel a only. The redundant within-lineage panel is omitted.
robust = tables['fig3a_hit_robustness.csv']
assert len(robust) == 68 and set(robust['check']) == {'Y-status definition', 'one model per donor', 'covariates'}
variant_order = robust.variant.drop_duplicates().tolist()
short_variants = {
    'primary call': 'Primary call', 'RNA and DNA agree only': 'RNA + DNA agree only',
    'RNA only': 'RNA only', 'DNA only': 'DNA only', 'one-copy controls only': 'One-copy controls',
    'relative Y gain removed': 'Relative Y gain removed', 'lenient RNA rule': 'Lenient RNA rule',
    'XIST-high lines removed': 'XIST-high removed', 'threshold-stable calls only': 'Threshold-stable calls',
    'first model per donor': 'One model per donor', 'lineage only': 'Lineage only',
    'lineage + library (primary)': 'Lineage + library', '+ WGD, TP53, aneuploidy': '+ WGD, TP53, aneuploidy',
    '+ LOH fraction, ploidy': '+ LOH fraction, ploidy', '+ screen quality': '+ Screen quality',
    '+ Cas9 activity, doubling time': '+ Cas9 activity, doubling time',
    'primary disease instead of lineage': 'Disease instead of lineage',
}
check_colours = {'Y-status definition': BLACK, 'one model per donor': PURPLE, 'covariates': TEAL}
for i, gene in enumerate(['DDX3X', 'EIF1AX', 'RPS4X', 'ZFX']):
    b = axes(4.47 + i * .91, 5.37, .76, 1.96)
    q = robust[robust.symbol.eq(gene)].set_index('variant').loc[variant_order]
    low, high = q.estimate - 1.96 * q.se, q.estimate + 1.96 * q.se
    for j, row in enumerate(q.itertuples()):
        b.errorbar(row.estimate, j, xerr=1.96 * row.se, fmt='o',
                   color=check_colours[row.check], markersize=2.8, markeredgewidth=0,
                   elinewidth=.72, capsize=0)
    # Preserve free x scales; a break between panels avoids implying equal scales.
    pad = (high.max() - low.min()) * .13
    b.set_xlim(low.min() - pad, high.max() + pad)
    b.set_ylim(len(q) - .4, -.6)
    b.set_yticks(range(len(q)), [short_variants[v] for v in variant_order] if i == 0 else [])
    b.tick_params(axis='y', length=0, labelsize=6.9, pad=4)
    b.tick_params(axis='x', labelsize=7.2)
    b.xaxis.set_major_locator(MaxNLocator(nbins=2, min_n_ticks=2))
    b.spines['left'].set_visible(False)
    b.text(.5, 1.028, gene, transform=b.transAxes, ha='center', va='bottom',
           fontsize=8.1, weight='semibold', fontstyle='italic')
label(2.67, 7.63, 'b')
text(6.22, 5.04, 'Chronos effect: LOY − Y-positive (95% CI)', ha='center', fontsize=7.6)
fig.legend(handles=[Line2D([], [], color=c, marker='o', markersize=3, lw=.75, label=s)
                    for s, c in [('Y-status definition', BLACK), ('One model/donor', PURPLE), ('Covariates', TEAL)]],
           loc='center', bbox_to_anchor=(6.22 / WIDTH, 4.85 / HEIGHT), ncol=3,
           fontsize=7, handlelength=1.25, handletextpad=.4, columnspacing=1.2)

# c: Nested prediction, using the original fitted values and colour mapping.
c = axes(.47, 2.99, 3.43, 1.60)
accuracy = tables['fig8a_nested_accuracy.csv']
colour_by_gene = {'DDX3X': BLACK, 'EIF1AX': RED, 'RPS4X': TEAL, 'ZFX': BLUE, 'COMPOSITE': PURPLE}
named = list(colour_by_gene)
for gene, frame in sorted(accuracy.groupby('dependency'), key=lambda item: item[0] in named):
    frame = frame.set_index('model_short').loc[['m0', 'm1', 'm2', 'm3']]
    c.plot(range(4), frame.r2_mean, color=colour_by_gene.get(gene, GREY), lw=.85, marker='o', ms=2.7)
offsets = {'DDX3X': -.022, 'EIF1AX': .011, 'ZFX': 0, 'RPS4X': 0, 'COMPOSITE': .012}
for gene in named:
    y = accuracy.loc[accuracy.dependency.eq(gene) & accuracy.model_short.eq('m3'), 'r2_mean'].item()
    c.annotate(gene, (3, y), xytext=(3.20, y + offsets[gene]), ha='left', va='center',
               color=colour_by_gene[gene], fontsize=7,
               arrowprops={'arrowstyle': '-', 'lw': .45, 'color': colour_by_gene[gene]})
c.axhline(0, color=BLACK, lw=.55, ls=(0, (1, 2)))
c.set(xlim=(-.18, 4.47), ylim=(-.13, .70), yticks=[0, .2, .4, .6], ylabel='Cross-validated R²')
c.set_xticks(range(4), ['Lineage\n+ library', '+ Y-status\ncategory', '+ 7 Y-gene\nexpression', '+ broader\nexpression'])
c.tick_params(axis='x', labelsize=7.1)
label(.04, 4.67, 'c')

# d: Stable feature selection, preserving both observed and shuffled-label counts.
d = axes(5.34, 2.99, 2.62, 1.60)
selected = tables['fig8c_y_status_forest.csv']
assert selected.loc[selected.dataset.eq('crispr_effect'), 'stable'].item() == 4
assert selected.loc[~selected.dataset.eq('crispr_effect'), 'stable'].sum() == 0
dataset_order = ['CRISPR dependencies', 'REP.PRIMARY', 'REP.1M', 'REP.300', 'PRISM secondary', 'GDSC1', 'GDSC2', 'CTRP v2']
frame = selected.set_index('dataset_label').loc[dataset_order]
ys = np.arange(len(frame))
for column, shift, colour, legend in [('stable', .175, BLACK, 'Observed Y status'),
                                    ('mean_null_confirmed', -.175, GREY, 'Shuffled Y (mean)')]:
    d.barh(ys + shift, frame[column], height=.325, color=colour, label=legend)
    for y, value in zip(ys + shift, frame[column]):
        d.text(value + .05, y, f'{value:.1f}'.removesuffix('.0'), va='center', fontsize=6.5)
d.set(yticks=ys, yticklabels=dataset_order, ylim=(7.6, -.6), xlim=(0, 4.8), xticks=range(5), xlabel='Features confirmed')
d.tick_params(axis='y', labelsize=7, length=0)
d.legend(loc='upper right', bbox_to_anchor=(1.03, -.23), ncol=2, fontsize=6.9,
         handlelength=1, columnspacing=.8, handletextpad=.45, borderaxespad=0)
label(4.03, 4.67, 'd')

# e: CIN thirds are a display stratification of the saved continuous interaction.
e = axes(.59, .50, 4.88, 1.75)
points = tables['fig5d_arhgef2_cin_points.csv']
arh_hit = tables['fig5d_arhgef2_discovery.csv'].iloc[0]
assert len(points) == 473 and points.groupby('y_status').size().to_dict() == {'Y-negative': 183, 'Y-positive': 290}
cuts = points.cin.quantile([1 / 3, 2 / 3]).to_numpy()
rng = np.random.default_rng(20261007)
for i, cin_group in enumerate(['Low CIN', 'Middle CIN', 'High CIN']):
    for j, (status, colour) in enumerate([('Y-positive', BLUE), ('Y-negative', RED)]):
        vals = points.loc[points.cin_group.eq(cin_group) & points.y_status.eq(status), 'adjusted_dependency'].to_numpy()
        x = i + (-.19 if j == 0 else .19)
        e.scatter(x + rng.uniform(-.1131, .1131, len(vals)), vals, s=7.5, color=colour, alpha=.36, edgecolors='none', zorder=2)
        e.boxplot([vals], positions=[x], widths=.29, showfliers=False, patch_artist=True, manage_ticks=False,
                  boxprops={'facecolor': colour, 'edgecolor': colour, 'alpha': .22, 'lw': .8},
                  medianprops={'color': colour, 'lw': 1.3}, whiskerprops={'color': colour, 'lw': .7},
                  capprops={'color': colour, 'lw': .7}, zorder=3)
        e.text(x, 3.07, f'{len(vals)}', ha='center', fontsize=6.7, color=colour)
e.set(xlim=(-.55, 2.55), ylim=(-3.2, 3.5), xticks=[0, 1, 2], yticks=[-3, -1, 1, 3],
      xticklabels=[f'Low CIN\n≤{cuts[0]:.3f}', f'Middle CIN\n{cuts[0]:.3f}–{cuts[1]:.3f}', f'High CIN\n>{cuts[1]:.3f}'],
      ylabel='Adjusted ARHGEF2 dependency\n(rank-normalised)')
e.tick_params(axis='x', labelsize=7.1)
e.text(0, 1.035, f'LOY × CIN: FDR = {arh_hit.interaction_fdr_class:.3f}', transform=e.transAxes, fontsize=7.5, va='bottom')
e.legend(handles=[Patch(facecolor=c, edgecolor=c, label=s, alpha=.65)
                  for s, c in [('Y-positive', BLUE), ('LOY', RED)]],
         loc='lower right', bbox_to_anchor=(1.02, 1.015), ncol=2, fontsize=7,
         handlelength=.85, handletextpad=.4, columnspacing=.7, borderaxespad=0)
label(.04, 2.51, 'e')

# f: Untransformed CIB1 scores and the saved within-subset FDR.
f = axes(6.12, .50, 1.84, 1.75)
cib = tables['fig5e_cib1_points.csv']
cib_hit = tables['fig5e_cib1_discovery.csv'].iloc[0]
assert cib.groupby('y_status').size().to_dict() == {'LOY': 8, 'Y-positive': 10}
assert np.isclose(cib.loc[cib.y_status.eq('LOY'), 'chronos'].mean() - cib.loc[cib.y_status.eq('Y-positive'), 'chronos'].mean(), cib_hit.inside, atol=1e-9)
rng = np.random.default_rng(20261007)
for i, (status, colour) in enumerate([('Y-positive', BLUE), ('LOY', RED)]):
    vals = cib.loc[cib.y_status.eq(status), 'chronos'].to_numpy()
    f.boxplot([vals], positions=[i], widths=.49, showfliers=False, patch_artist=True, manage_ticks=False,
              boxprops={'facecolor': colour, 'edgecolor': colour, 'alpha': .20, 'lw': .8},
              medianprops={'color': colour, 'lw': 1.4}, whiskerprops={'color': colour, 'lw': .8},
              capprops={'color': colour, 'lw': .8})
    f.scatter(i + rng.uniform(-.14, .14, len(vals)), vals, color=colour, s=19,
              alpha=.8, edgecolors='white', linewidth=.4, zorder=3)
f.axhline(0, color='#999999', lw=.6, ls=(0, (3, 3)), zorder=0)
f.set(xlim=(-.55, 1.55), ylim=(-.5, 1.52), xticks=[0, 1],
      xticklabels=['Y-positive\nn = 10', 'LOY\nn = 8'], yticks=[-.4, 0, .4, .8, 1.2],
      ylabel='CIB1 gene effect (Chronos)', xlabel='Bladder/urinary tract')
f.tick_params(axis='x', labelsize=7.1)
f.set_xlabel(f.get_xlabel(), fontsize=7.4)
f.plot([0, 0, 1, 1], [1.28, 1.33, 1.33, 1.28], color=BLACK, lw=.65)
f.text(.5, 1.36, f'FDR = {cib_hit.own_fdr:.3f}', ha='center', va='bottom', fontsize=7.6)
label(5.62, 2.51, 'f')

for extension in ('pdf', 'svg', 'png'):
    fig.savefig(FIG / f'{STEM}.{extension}', dpi=300, facecolor='white', bbox_inches='tight', pad_inches=.035)
plt.close(fig)
provenance = {
    'status': 'Selected as CH3-S3.4.1-S01 in Results_edit.md on 8 October 2026; superseded source figures are preserved.',
    'working_id': 'CH3-S3.4.1-S01', 'selected_pdf': 'figures/chapter3/selected/CH3_S3_4_1_S01_depmap_combined_checks/ch03_supp_depmap_combined_checks.pdf',
    'font': font, 'figure_inches_before_tight_crop': [WIDTH, HEIGHT], 'source_files': hashes,
    'panel_mapping': {
        'a': 'CH3-S3.4.1-S01: complete Y-call evidence panel',
        'b': 'CH3-S3.4.2-S01: panel a only; within-lineage panel b omitted',
        'c': 'CH3-S3.3.3-S04: panel a (nested random forests)',
        'd': 'CH3-S3.3.3-S04: panel c (Boruta feature selection); formerly combined panel e',
        'e': 'CH3-S3.3.2-S02: panel b (ARHGEF2 / CIN); formerly combined panel g',
        'f': 'CH3-S3.3.2-S02: panel c (CIB1 bladder subset); formerly combined panel h',
    },
    'statistics_refitted': False, 'all_source_values_preserved': True,
    'excluded': ['within-lineage robustness panel',
                 'former combined panel d: confounder-matched forest check (author request, 8 October 2026)',
                 'former combined panel f: approximate detectable interaction sizes (author request, 8 October 2026)',
                 'previously removed genetic cancer-type interactions',
                 'previously removed drug interactions'],
    'display_changes': ['Six panels in a compact three-row layout', 'Shorter labels without changing category definitions',
                        'Free x scales for the four robustness gene facets, with 95% intervals unchanged',
                        'Counts above ARHGEF2 boxes denote individual models'],
    'robustness_source_note': '5 October and 7 October fig3a source files have identical SHA-256 hashes.',
    'exports': {extension: hashlib.sha256((FIG / f'{STEM}.{extension}').read_bytes()).hexdigest()
                for extension in ('pdf', 'svg', 'png')},
}
(DATA / 'provenance.json').write_text(json.dumps(provenance, indent=2), encoding='utf-8')
(DATA / 'README.md').write_text("# Combined DepMap supplementary figure\n\nSelected as **CH3-S3.4.1-S01** in Results_edit.md section 3.3.1 on 8 October 2026. The six-panel figure replaces the four separate Y-call evidence, robustness, random-forest and genetic-context figure blocks.\n\nPanels a-f show Y-call evidence, robustness of the four core dependencies, cross-validated prediction, Boruta feature counts, ARHGEF2/CIN and CIB1 in bladder/urinary tract models. Former combined panels d (confounder-matched forest check) and f (approximate detectable interaction sizes) were removed at the author's request. The redundant within-lineage robustness panel remains omitted; those data appear in main F01g.\n\nAll retained panel values, intervals, counts, FDRs and eligibility definitions come from the eight copied source tables listed in provenance.json. No models were refitted. Superseded source tables and figure exports are retained for historical documents.\n\n[Caption](../../../../Chapter3_figure_captions.md#supp--combined-depmap-evidence-and-robustness)\n\n[Selected PDF](../../selected/CH3_S3_4_1_S01_depmap_combined_checks/ch03_supp_depmap_combined_checks.pdf)\n\n[Reproduction script](../../code/ch03_supp_depmap_combined_checks.py)\n", encoding='utf-8')
print(FIG / f'{STEM}.png')
