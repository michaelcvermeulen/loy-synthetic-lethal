"""Boxplot view of the same adjusted CRISPR values shown in panel A."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Patch

root = Path(__file__).resolve().parents[3]
out = root / 'figures/chapter3'
source = out / 'source_data/arhgef2_cin_exploratory/points.csv'
data = out / 'source_data/arhgef2_cin_boxplot'
data.mkdir(parents=True, exist_ok=True)
p = pd.read_csv(source).query("assay == 'CRISPR'").copy()
assert len(p) == 473 and not p.adjusted_dependency.isna().any()
assert p.groupby('y_status').size().to_dict() == {'Y-negative': 183, 'Y-positive': 290}
cuts = p.cin.quantile([1/3, 2/3]).to_numpy()
p['cin_group'] = pd.cut(p.cin, [-np.inf, *cuts, np.inf], labels=['Low CIN', 'Middle CIN', 'High CIN'])
assert p.cin_group.notna().all()
for font in ['Regular', 'SemiBold']:
    font_manager.fontManager.addfont(str(root / f'skills/thesis-figure-style/assets/fonts/OpenSans-{font}.ttf'))
plt.rcParams.update({'font.family': 'Open Sans', 'font.size': 9,
    'axes.spines.top': False, 'axes.spines.right': False, 'axes.linewidth': .65,
    'pdf.fonttype': 42, 'svg.fonttype': 'none', 'text.color': '#222222',
    'axes.labelcolor': '#222222'})
colours = {'Y-positive': '#3C5488', 'Y-negative': '#E64B35'}
rng = np.random.default_rng(20261007)
fig = plt.figure(figsize=(9, 4.7))
left = fig.add_axes([.105, .205, .24, .65])
right = fig.add_axes([.455, .205, .515, .65], sharey=left)
summary = []

def draw_box(ax, subset, group, x, width, context):
    vals = subset.adjusted_dependency.to_numpy()
    col = colours[group]
    ax.scatter(x + rng.uniform(-width*.39, width*.39, len(vals)), vals,
               s=10, c=col, alpha=.36, edgecolors='none', zorder=2)
    box = ax.boxplot([vals], positions=[x], widths=width, patch_artist=True,
        showfliers=False, manage_ticks=False,
        boxprops={'facecolor': col, 'edgecolor': col, 'alpha': .22, 'linewidth': 1},
        medianprops={'color': col, 'linewidth': 1.7},
        whiskerprops={'color': col, 'linewidth': 1},
        capprops={'color': col, 'linewidth': 1}, zorder=3)
    ax.text(x, 3.05, f'n = {len(vals)}', ha='center', fontsize=8, color=col)
    summary.append({'cin_group': context, 'y_status': group, 'n': len(vals),
        'q1': float(np.quantile(vals, .25)), 'median': float(np.median(vals)),
        'q3': float(np.quantile(vals, .75))})

for x, group in enumerate(colours):
    draw_box(left, p[p.y_status == group], group, x, .53, 'All CIN values')
left.set(xlim=(-.6, 1.6), ylim=(-3.2, 3.4), xticks=[0, 1],
         xticklabels=['Y-positive', 'Y-negative'], yticks=np.arange(-3, 4),
         ylabel='Adjusted ARHGEF2 gene effect\n(rank-normalised)')
left.text(.5, 1.055, 'All CIN values', transform=left.transAxes, ha='center', weight='semibold')
for i, context in enumerate(['Low CIN', 'Middle CIN', 'High CIN']):
    for j, group in enumerate(colours):
        draw_box(right, p[(p.cin_group == context) & (p.y_status == group)],
                 group, i + (-.19 if j == 0 else .19), .29, context)
right.set(xlim=(-.55, 2.55), xticks=[0, 1, 2], xticklabels=[
    f'Low CIN\n≤{cuts[0]:.3f}', f'Middle CIN\n>{cuts[0]:.3f}–{cuts[1]:.3f}', f'High CIN\n>{cuts[1]:.3f}'])
right.tick_params(axis='y', labelleft=False)
right.text(.5, 1.055, 'CIN thirds across the same 473 models', transform=right.transAxes,
           ha='center', weight='semibold')
fig.legend(handles=[Patch(facecolor=c, edgecolor=c, label=g, alpha=.65) for g,c in colours.items()],
           loc='upper center', bbox_to_anchor=(.60, 1.0), ncol=2, frameon=False,
           fontsize=8.5, handlelength=1.2, columnspacing=2)
fig.text(.025, .88, 'a', weight='semibold', fontsize=12)
fig.text(.385, .88, 'b', weight='semibold', fontsize=12)
fig.text(.105, .065, 'Lower values indicate stronger dependency. Boxes: median and IQR; whiskers: 1.5 × IQR; dots: individual models.', fontsize=8)
fig.text(.105, .027, 'CIN thirds are for display only; the original interaction test treated CIN as continuous.', fontsize=8, color='#555555')
stem = 'ch03_supp_arhgef2_cin_boxplot'
for ext in ['png', 'pdf', 'svg']:
    fig.savefig(out / f'{stem}.{ext}', dpi=250, facecolor='white')
plt.close(fig)
p.to_csv(data / 'plotted_points.csv', index=False)
pd.DataFrame(summary).to_csv(data / 'box_summaries.csv', index=False)
(data / 'provenance.json').write_text(json.dumps({
    'source': str(source), 'source_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
    'n': len(p), 'cin_tertile_cutoffs': cuts.tolist(),
    'y_variable': 'Identical adjusted rank-normalised ARHGEF2 CRISPR values from original exploratory panel A',
    'adjustment': 'Lineage and screening library contributions removed, centred at assay-wide covariate means; CIN, Y and Y:CIN retained',
    'grouping': 'CIN tertiles computed over all 473 CRISPR models, shared by both Y groups; descriptive only',
    'boxplot': 'Median, 25th–75th percentiles; whiskers to most extreme observations within 1.5 IQR; all observations shown as points',
    'inference': 'No new hypothesis tests; the original continuous interaction P is not a pooled or between-bin comparison P',
    'status': 'Exploratory preview; original figure and Results unchanged',
    'font': 'Bundled Open Sans'
}, indent=2) + '\n', encoding='utf-8')
(data / 'README.md').write_text('''# ARHGEF2 CRISPR boxplot preview

Uses exactly the adjusted rank-normalised measurements from panel A of the preceding exploratory display. The left panel pools CIN values; the right panel groups models into shared CIN tertiles to show the pattern across CIN. These groups are descriptive and do not replace the original continuous interaction model. Boxes show medians and interquartile ranges, whiskers reach the most extreme observations within 1.5 IQR, and every cell line is plotted. Counts appear above each box. No P value from the continuous model is reassigned to a boxplot comparison.

Source measurements, group assignments, summaries and exact cutoffs are saved here. Original exploratory model definitions and validation are in ../arhgef2_cin_exploratory/. No Results prose or selected figure was changed.
''', encoding='utf-8')
print(pd.DataFrame(summary).to_string(index=False))
print('CIN tertile cutoffs:', cuts)
print(out / f'{stem}.png')
