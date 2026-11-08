"""Descriptive RCOR3 follow-up using existing measurements and saved inference."""
from pathlib import Path
import hashlib
import json
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

ROOT = Path(__file__).resolve().parents[3]
FIG = ROOT / 'figures/chapter3'
DATA = FIG / 'source_data/rcor3_exploratory'
STEM = 'ch03_supp_rcor3_exploratory'
sys.path.insert(0, str(ROOT / 'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
font = apply_style(8.4)
BLUE, RED, BLACK = '#3C5488', '#E64B35', '#222222'
COLOURS = {'Y-positive': BLUE, 'LOY': RED}
crispr = pd.read_csv(DATA / 'crispr_points.csv')
rnai = pd.read_csv(DATA / 'rnai_points.csv')
pan = pd.read_csv(DATA / 'saved_gene_pan_results.csv').iloc[0]
rnai_result = pd.read_csv(DATA / 'saved_gene_pan_rnai.csv').iloc[0]
continuous = pd.read_csv(DATA / 'saved_gene_pan_continuous.csv').iloc[0]
forest = pd.read_csv(DATA / 'saved_forest_confounders.csv').iloc[0]
assert len(crispr) == pan.n == 568 and len(rnai) == rnai_result.n == 306
assert crispr.dependency_probability.ge(.5).sum() == 1
assert crispr.loc[crispr.dependency_probability.ge(.5), 'cell_line'].item() == 'MAVER1'
fig = plt.figure(figsize=(9.0, 6.2))
a = fig.add_axes([.075, .575, .225, .34])
b = fig.add_axes([.075, .10, .225, .34])
c = fig.add_axes([.57, .10, .405, .815])
for ax in (a, b, c):
    ax.tick_params(length=3, pad=3, labelsize=7.7)
for x, y, letter in [(.018, .974, 'a'), (.018, .493, 'b'), (.347, .974, 'c')]:
    fig.text(x, y, letter, va='top', ha='left', fontsize=11, weight='semibold')
fig.legend(handles=[Patch(facecolor=colour, edgecolor=colour, label=status, alpha=.65)
                    for status, colour in COLOURS.items()],
           loc='upper center', bbox_to_anchor=(.52, .99), ncol=2, frameon=False,
           handlelength=1.1, columnspacing=1.2, fontsize=8.2)


def pooled(ax, df, score_label, saved):
    rng = np.random.default_rng(20261008)
    for i, (status, colour) in enumerate(COLOURS.items()):
        vals = df.loc[df.y_status.eq(status), 'score'].to_numpy()
        ax.scatter(i + rng.uniform(-.17, .17, len(vals)), vals, s=9,
                   color=colour, alpha=.36, edgecolors='none', zorder=2)
        ax.boxplot([vals], positions=[i], widths=.48, showfliers=False,
                   patch_artist=True, manage_ticks=False,
                   boxprops={'facecolor': colour, 'edgecolor': colour, 'alpha': .20, 'lw': .8},
                   medianprops={'color': colour, 'lw': 1.6}, whiskerprops={'color': colour, 'lw': .8},
                   capprops={'color': colour, 'lw': .8}, zorder=3)
    counts = df.y_status.value_counts()
    ax.axhline(0, color='#999999', ls=(0, (3, 3)), lw=.65, zorder=0)
    ax.set(xlim=(-.62, 1.62), xticks=[0, 1],
           xticklabels=[f'Y-positive\nn = {counts["Y-positive"]}', f'LOY\nn = {counts["LOY"]}'],
           ylabel=score_label)
    ax.text(.5, 1.09, f'Adjusted effect = {saved.estimate:+.3f}', transform=ax.transAxes,
            ha='center', va='bottom', fontsize=8)
    ax.text(.5, 1.018, f'P = {saved.p:.3f}; FDR = {saved.fdr:.3f}', transform=ax.transAxes,
            ha='center', va='bottom', fontsize=8)


pooled(a, crispr, 'RCOR3 gene effect (Chronos)', pan)
a.set_ylim(-.66, .64)
a.set_yticks([-.5, 0, .5])
maver = crispr.loc[crispr.cell_line.eq('MAVER1')].iloc[0]
a.annotate('MAVER1', xy=(0, maver.score), xytext=(.65, -.55),
           fontsize=7.3, ha='left', va='center',
           arrowprops={'arrowstyle': '-', 'lw': .7, 'color': BLUE}, color=BLUE)
pooled(b, rnai, 'RCOR3 dependency (RNAi)', rnai_result)
b.set_ylim(-.44, .50)
b.set_yticks([-.4, -.2, 0, .2, .4])

counts = crispr.groupby(['lineage', 'y_status']).size().unstack(fill_value=0)
eligible = counts.loc[counts.min(axis=1).ge(5)].copy()
eligible['total'] = eligible.sum(axis=1)
eligible = eligible.sort_values(['total', 'lineage'], ascending=[False, True])
assert len(eligible) == 14
display_names = {'Bladder/Urinary Tract': 'Bladder/urinary tract', 'Esophagus/Stomach': 'Oesophagus/stomach',
                 'Head and Neck': 'Head and neck', 'Soft Tissue': 'Soft tissue', 'Biliary Tract': 'Biliary tract'}
rng = np.random.default_rng(20261008)
for i, lineage in enumerate(eligible.index):
    for j, (status, colour) in enumerate(COLOURS.items()):
        vals = crispr.loc[crispr.lineage.eq(lineage) & crispr.y_status.eq(status), 'score'].to_numpy()
        y = i + (-.18 if j == 0 else .18)
        c.scatter(vals, y + rng.uniform(-.105, .105, len(vals)), s=8, color=colour,
                   alpha=.38, edgecolors='none', zorder=2)
        c.boxplot([vals], positions=[y], widths=.25, orientation='horizontal',
                   showfliers=False, patch_artist=True, manage_ticks=False,
                   boxprops={'facecolor': colour, 'edgecolor': colour, 'alpha': .23, 'lw': .65},
                   medianprops={'color': colour, 'lw': 1.3}, whiskerprops={'color': colour, 'lw': .65},
                   capprops={'color': colour, 'lw': .65}, zorder=3)
c.axvline(0, color='#999999', ls=(0, (3, 3)), lw=.65, zorder=0)
c.set(ylim=(13.6, -.6), xlim=(-.57, .60), xticks=[-.4, -.2, 0, .2, .4, .6],
      yticks=range(14), yticklabels=[f'{display_names.get(lineage, lineage)}  ({row["Y-positive"]}/{row["LOY"]})'
                                  for lineage, row in eligible.iterrows()],
      xlabel='RCOR3 gene effect (Chronos)')
c.tick_params(axis='y', length=0, pad=5, labelsize=7.6)
c.text(0, 1.019, 'Lineages with ≥5 models per Y group; counts: Y-positive / LOY',
       transform=c.transAxes, fontsize=7.2, ha='left', va='bottom')

for ext in ('png', 'pdf', 'svg'):
    fig.savefig(FIG / f'{STEM}.{ext}', dpi=300, bbox_inches='tight', pad_inches=.04)
plt.close(fig)

summaries = []
for assay, df in [('CRISPR', crispr), ('RNAi', rnai)]:
    summary = df.groupby('y_status').score.agg(n='size', mean='mean', median='median', sd='std', minimum='min', maximum='max').reset_index()
    summary.insert(0, 'assay', assay)
    summaries.append(summary)
pd.concat(summaries).to_csv(DATA / 'group_summaries.csv', index=False)
lineage_summary = crispr.groupby(['lineage', 'y_status']).score.agg(n='size', mean='mean', median='median').reset_index()
lineage_summary['displayed'] = lineage_summary.lineage.isin(eligible.index)
lineage_summary.to_csv(DATA / 'lineage_summaries.csv', index=False)
crispr.sort_values('score').head(15).to_csv(DATA / 'most_negative_crispr_scores.csv', index=False)
summary = {
    'status': 'Exploratory follow-up of the largest average shadow-win fraction among random control genes; not a new discovery.',
    'font': font, 'statistical_tests_refitted': False, 'within_lineage_inference': 'Descriptive distributions only; no new P values or FDRs.',
    'primary_crispr': {key: float(pan[key]) for key in ['estimate', 'ci_low', 'ci_high', 'p', 'fdr', 'n', 'mean_negative', 'mean_positive']},
    'rna_i': {key: float(rnai_result[key]) for key in ['estimate', 'ci_low', 'ci_high', 'p', 'fdr', 'n']},
    'continuous_y_score': {key: float(continuous[key]) for key in ['estimate', 'p', 'fdr', 'n']},
    'forest': {key: float(forest[key]) for key in ['success_seed1', 'success_seed2', 'y_importance', 'y_rank']},
    'dep_probability_at_least_0_5': {'Y-positive': 1, 'LOY': 0},
    'only_dependent_model': {'name': 'MAVER1', 'Y_status': 'Y-positive', 'chronos': float(maver.score),
                            'dependency_probability': float(maver.dependency_probability)},
    'context_discovery_eligibility_note': 'Only one line reaches dependency probability 0.5; the context discovery selection requires at least five dependent lines.',
    'input_files': {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in DATA.glob('*.csv')},
    'plot': 'Unadjusted individual scores. P values, confidence intervals and effect estimates come from the saved adjusted pan-cancer models.',
    'selected_lineages': eligible.index.tolist(),
}
(DATA / 'provenance.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
(DATA / 'README.md').write_text(
    '# Exploratory RCOR3 follow-up\n\n'
    'Requested after RCOR3 had the largest average shadow-win fraction among the 200 random non-hit genes. '
    'It won 23/30 and 26/30 donor subsamples, so did not reach 80% under both seeds.\n\n'
    'Plots show unadjusted individual measurements; all reported inferential statistics are copied from '
    'the saved 7 October 2026 HPC export. The primary CRISPR model adjusts for lineage and library with donor blocking; '
    'RNAi adjusts for lineage with donor blocking. No model was refitted. Within-lineage plots are descriptive '
    'and require five measured models per Y group. CRISPR and RNAi scales differ.\n\n'
    'RCOR3 has only one primary-cohort model with dependency probability ≥0.5, below the five-dependent-model '
    'selection criterion for the context discovery family. The selected model is Y-positive MAVER1. '
    'This does not prove why the random forest assigned Y status predictive importance.\n\n'
    '[Caption](../../../../Chapter3_figure_captions.md#supp--exploratory-rcor3-dependency-comparisons)\n\n'
    '[R extraction](../../code/ch03_supp_rcor3_exploratory_extract.R) | '
    '[Python plot](../../code/ch03_supp_rcor3_exploratory.py)\n', encoding='utf-8')
print(FIG / f'{STEM}.png')
