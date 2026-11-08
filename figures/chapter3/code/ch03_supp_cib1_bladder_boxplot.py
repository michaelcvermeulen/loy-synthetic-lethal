"""Measured CIB1 Chronos scores in bladder/urinary tract models; exploratory preview."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager

root = Path(__file__).resolve().parents[3]
out = root / 'figures/chapter3'
data = out / 'source_data/cib1_bladder_boxplot'
source = root / 'research/depmap_loy_reanalysis_v1/loy_reanalysis_v1_export_2026-10-07/derived/loy_reanalysis_v1'
p = pd.read_csv(data / 'plotted_points.csv')
hit = pd.read_csv(data / 'saved_within_cancer_result.csv').iloc[0]
assert len(p) == 18 and p.groupby('y_status').size().to_dict() == {'LOY': 8, 'Y-positive': 10}
assert np.isclose(p.loc[p.y_status.eq('LOY'), 'chronos'].mean() - p.loc[p.y_status.eq('Y-positive'), 'chronos'].mean(), hit.inside, atol=1e-9)
for font in ['Regular', 'SemiBold']:
    font_manager.fontManager.addfont(str(root / f'skills/thesis-figure-style/assets/fonts/OpenSans-{font}.ttf'))
plt.rcParams.update({'font.family': 'Open Sans', 'font.size': 9,
    'axes.spines.top': False, 'axes.spines.right': False, 'axes.linewidth': .7,
    'pdf.fonttype': 42, 'svg.fonttype': 'none', 'text.color': '#222222', 'axes.labelcolor': '#222222'})
fig = plt.figure(figsize=(4.3, 3.8))
ax = fig.add_axes([.185, .21, .77, .66])
colours = {'Y-positive': '#3C5488', 'LOY': '#E64B35'}
rng = np.random.default_rng(20261007)
summaries = []
for i, (group, colour) in enumerate(colours.items()):
    vals = p.loc[p.y_status.eq(group), 'chronos'].to_numpy()
    ax.boxplot([vals], positions=[i], widths=.49, showfliers=False, manage_ticks=False, patch_artist=True,
        boxprops={'facecolor': colour, 'edgecolor': colour, 'alpha': .20, 'linewidth': 1},
        medianprops={'color': colour, 'linewidth': 1.6}, whiskerprops={'color': colour, 'linewidth': .9},
        capprops={'color': colour, 'linewidth': .9})
    ax.scatter(i + rng.uniform(-.14, .14, len(vals)), vals, s=28, c=colour, alpha=.8, edgecolors='white', linewidth=.45, zorder=3)
    summaries.append({'y_status': group, 'n': len(vals), 'mean': vals.mean(), 'median': np.median(vals),
                      'q1': np.quantile(vals, .25), 'q3': np.quantile(vals, .75), 'min': vals.min(), 'max': vals.max()})
ax.axhline(0, color='#999999', linestyle=(0, (3, 3)), linewidth=.65, zorder=0)
ax.set(xlim=(-.58, 1.58), ylim=(-.50, 1.52), xticks=[0, 1],
       xticklabels=['Y-positive\nn = 10', 'LOY\nn = 8'],
       yticks=[-.4, 0, .4, .8, 1.2], ylabel='CIB1 gene effect (Chronos)')
ax.plot([0, 0, 1, 1], [1.28, 1.33, 1.33, 1.28], color='#333333', linewidth=.75)
ax.text(.5, 1.36, f'Within-cancer FDR = {hit.own_fdr:.3f}', ha='center', va='bottom', fontsize=8.5)
fig.text(.185, .08, 'Bladder/urinary tract models', fontsize=8.5)
fig.text(.185, .03, 'Lower scores indicate stronger dependency.', fontsize=8, color='#555555')
stem = 'ch03_supp_cib1_bladder_boxplot'
for ext in ['png', 'pdf', 'svg']:
    fig.savefig(out / f'{stem}.{ext}', dpi=300, facecolor='white')
plt.close(fig)
pd.DataFrame(summaries).to_csv(data / 'box_summaries.csv', index=False)
(data / 'provenance.json').write_text(json.dumps({
    'source_files': {str(source / name): hashlib.sha256((source / name).read_bytes()).hexdigest()
                     for name in ['cohort.csv', 'matrices/crispr_effect.rds', 'power_alone_hits.csv']},
    'extractor': 'figures/chapter3/code/ch03_supp_cib1_bladder_boxplot_extract.R',
    'plotter': str(Path(__file__).resolve()),
    'outcome': 'Untransformed measured CIB1 (10519) Chronos gene-effect scores',
    'cohort': 'Primary male Y-classified CRISPR models in Oncotree lineage Bladder/Urinary Tract; 10 Y-positive and 8 Y-negative',
    'saved_raw_effect': float(hit.inside), 'saved_rank_based_p': float(hit.inside_p),
    'saved_within_cancer_fdr': float(hit.own_fdr),
    'inference': 'Existing globally rank-normalised limma screen, with BH correction across 6089 genes within this cancer subset. Library covariates are constant in these 18 models. No new statistical tests.',
    'display': 'Raw scores, boxes median/IQR, whiskers to furthest observations within 1.5 IQR, every model plotted; dashed line at zero.',
    'status': 'Standalone source preview, incorporated into selected Supplement S02 as panel c.',
    'font': 'Bundled Open Sans'
}, indent=2) + '\n', encoding='utf-8')
(data / 'README.md').write_text('''# CIB1 bladder/urinary tract boxplot preview

Individual CIB1 Chronos measurements are extracted from the 7 October 2026 HPC export using the original primary cohort and cancer-lineage labels. Group sizes and the raw mean difference are checked against the saved within-cancer result. The plot shows untransformed scores; its FDR annotation is the saved rank-based screen result, corrected across the 6,089 genes tested within this subset. The library covariates are constant in these models. No interaction test or new hypothesis test is used.

The standalone preview is incorporated into selected Supplement S02 as panel c. Extraction and plotting scripts are in ../../code/. The caption is held in [central captions](../../../../Chapter3_figure_captions.md).
''', encoding='utf-8')
print(pd.DataFrame(summaries).to_string(index=False))
print(out / f'{stem}.png')
