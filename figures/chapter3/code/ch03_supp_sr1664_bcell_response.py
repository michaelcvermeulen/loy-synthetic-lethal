"""Exploratory SR-1664 boxplot from the frozen within-subset screen."""
from pathlib import Path
import hashlib
import json
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[3]
FIG = ROOT / 'figures/chapter3'
DATA = FIG / 'source_data/sr1664_bcell_response'
STEM = 'ch03_supp_sr1664_bcell_response'
sys.path.insert(0, str(ROOT / 'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
font = apply_style(9)
points = pd.read_csv(DATA / 'plotted_points.csv')
hit = pd.read_csv(DATA / 'saved_test.csv').iloc[0]
assert len(points) == 24 and not points.ModelID.duplicated().any()
groups = [('Y-positive', '#3C5488', 16), ('LOY', '#E64B35', 8)]
fig, ax = plt.subplots(figsize=(4.0, 4.5))
fig.subplots_adjust(left=.23, right=.97, top=.97, bottom=.23)
rng = np.random.default_rng(20261008)
for x, (group, colour, n) in enumerate(groups):
    vals = points.loc[points.y_status.eq(group), 'response'].to_numpy()
    assert len(vals) == n
    ax.boxplot([vals], positions=[x], widths=.46, manage_ticks=False,
               showfliers=False, patch_artist=True,
               boxprops=dict(facecolor=colour, edgecolor=colour, alpha=.20, lw=.8),
               medianprops=dict(color=colour, lw=1.5),
               whiskerprops=dict(color=colour, lw=.8), capprops=dict(color=colour, lw=.8))
    ax.scatter(x + rng.uniform(-.12, .12, len(vals)), vals, s=27,
               color=colour, edgecolor='white', linewidth=.4, alpha=.9, zorder=3)
low, high = min(points.response.min(), 0), max(points.response.max(), 0)
spread = high - low
ax.set_ylim(low - .08 * spread, high + .48 * spread)
ax.set_xlim(-.55, 1.55)
ax.axhline(0, color='#999999', linestyle=(0, (3, 3)), linewidth=.7, zorder=0)
ax.set_xticks([0, 1], ['Y-positive\nn = 16', 'LOY\nn = 8'])
ax.set_ylabel('SR-1664 response (log₂ fold change)')
ax.set_xlabel('Mature B-cell neoplasms\nPRISM REP.1M', labelpad=9)
trans = ax.get_xaxis_transform()
ax.plot([0, 0, 1, 1], [.77, .79, .79, .77], transform=trans, color='#222222', linewidth=.8)
ax.text(.5, .985, f'Δ = {hit.inside:+.3f}', transform=trans, ha='center', va='top')
ax.text(.5, .925, r'P = 6.50 × 10$^{-7}$', transform=trans, ha='center', va='top')
ax.text(.5, .865, f'FDR = {hit.inside_fdr_class:.3f}', transform=trans, ha='center', va='top')
fig.text(.60, .02, 'Lower response = greater drug sensitivity', ha='center', fontsize=8)
for ext in ['pdf', 'svg', 'png']:
    fig.savefig(FIG / f'{STEM}.{ext}', dpi=300, facecolor='white', bbox_inches='tight', pad_inches=.045)
plt.close(fig)
derived = ROOT / 'research/depmap_loy_reanalysis_v1/loy_reanalysis_v1_export_2026-10-07/derived/loy_reanalysis_v1'
inputs = [derived / p for p in ['cohort.rds', 'context_matrix.rds', 'matrices/prism_REP1M.rds']]
inputs += [ROOT / 'tables/chapter3/depmap_drug_within_subsets/within_subset_hits.csv',
           DATA / 'plotted_points.csv', DATA / 'saved_test.csv']
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
provenance = dict(status='Used in CH3-S3.3.3-F02 panel e; standalone preview retained',
                  compound=hit.compound_name, feature=hit.feature, dataset=hit.dataset,
                  context=hit.label, models=24, n_Y_positive=16, n_LOY=8,
                  response='PRISM log2 fold change; lower means greater sensitivity',
                  effect=float(hit.inside), p=float(hit.inside_p), fdr=float(hit.inside_fdr_class),
                  inference='Saved within-subset OLS P value on responses rank-normalised across eligible dataset models; original BH FDR across the disease class.',
                  models_refitted=False, raw_effect_matches_mean_difference=True, font=font,
                  input_sha256={p.relative_to(ROOT).as_posix():sha(p) for p in inputs},
                  exports={ext:sha(FIG / f'{STEM}.{ext}') for ext in ['pdf','svg','png']})
(DATA / 'provenance.json').write_text(json.dumps(provenance, indent=2)+'\n', encoding='utf-8')
(DATA / 'README.md').write_text(
    '# SR-1664 response in mature B-cell neoplasm models\n\n'
    'Exploratory display of the saved PRISM REP.1M within-subset association, requested 8 October 2026. '
    'Added to CH3-S3.3.3-F02 panel e on 8 October 2026; standalone preview retained.\n\n'
    'Points show the original log2 fold-change measurements for 16 Y-positive and eight LOY male models; '
    'boxes show medians and interquartile ranges and whiskers the extreme observations within 1.5 IQR. '
    'Lower responses indicate greater sensitivity. The positive effect indicates lower sensitivity in LOY. '
    'The raw group-mean difference reproduces the saved effect (+0.6672401518). '
    'P = 6.497567969e-7 is the saved OLS test of rank-normalised responses, not a new test on the plotted raw values. '
    'FDR = 0.03108631443 is the original BH value across compounds, datasets and subsets in the disease class. '
    'There is no interaction result in this plot; this comparison does not establish disease specificity.\n\n'
    'Reproduce with ch03_supp_sr1664_bcell_extract.R (Rscript --vanilla), then '
    'ch03_supp_sr1664_bcell_response.py in figures/chapter3/code. No models are refitted.\n', encoding='utf-8')
print(FIG / f'{STEM}.png')
