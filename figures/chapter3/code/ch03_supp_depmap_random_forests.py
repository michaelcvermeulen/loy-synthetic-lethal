"""Render S04 from frozen panel tables with embedded thesis fonts."""
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
DATA = FIG / 'source_data/depmap_random_forests'
sys.path.insert(0, str(ROOT / 'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
font = apply_style(9)
names = ['fig8a_nested_accuracy', 'fig8b_confounder_check', 'fig8c_y_status_forest']
accuracy, confounder, selected = [pd.read_csv(DATA / f'{name}.csv') for name in names]
assert len(confounder.query("role == 'hit'")) == 5
assert len(confounder.query("role != 'hit'")) == 200
assert selected.loc[selected.dataset == 'crispr_effect', 'stable'].item() == 4
assert selected.loc[selected.dataset != 'crispr_effect', 'stable'].sum() == 0

RED, GREY = '#E64B35', '#B8B8B8'
colours = {'DDX3X': '#222222', 'EIF1AX': RED, 'ZFX': '#3C5488',
           'RPS4X': '#00A087', 'COMPOSITE': '#8172B3', 'other': GREY}
fig = plt.figure(figsize=(9.2, 3.2))
a = fig.add_axes([.063, .30, .239, .64])
b = fig.add_axes([.364, .30, .225, .64])
c = fig.add_axes([.747, .30, .238, .64])
for ax in [a, b, c]:
    ax.tick_params(labelsize=7, length=3, pad=2)
for x, label in [(0.013, 'a'), (.312, 'b'), (.614, 'c')]:
    fig.text(x, .985, label, va='top', weight='semibold', fontsize=10)

order = ['m0', 'm1', 'm2', 'm3']
named = ['DDX3X', 'EIF1AX', 'ZFX', 'RPS4X', 'COMPOSITE']
for dependency, frame in sorted(accuracy.groupby('dependency'), key=lambda item: item[0] in named):
    frame = frame.set_index('model_short').loc[order]
    colour = colours.get(dependency, GREY)
    a.plot(range(4), frame.r2_mean, color=colour, lw=.8, marker='o', ms=2.6)
offsets = {'DDX3X': -.022, 'EIF1AX': .011, 'ZFX': 0, 'RPS4X': 0, 'COMPOSITE': .012}
for name in named:
    y = accuracy.loc[(accuracy.dependency == name) & (accuracy.model_short == 'm3'), 'r2_mean'].item()
    a.annotate(name, (3, y), xytext=(3.27, y + offsets[name]), color=colours[name],
               ha='left', va='center', fontsize=6.4,
               arrowprops={'arrowstyle': '-', 'lw': .5, 'color': colours[name]},
               annotation_clip=False)
a.axhline(0, color='#222222', lw=.6, ls=(0, (1, 2)))
a.set_xlim(-.3, 4.45)
a.set_ylim(-.13, .68)
a.set_yticks([0, .2, .4, .6])
a.set_xticks(range(4), ['lineage\n+ library', '+ Y\nstatus', '+ 7 Y\ngenes', '+ all\nexpression'])
a.set_ylabel('Out-of-fold R²', labelpad=3)

# ggplot2's original histogram uses 0.05-wide right-closed bins at boundary 0.
# Use its saved bin geometry so exact bin-edge observations remain unchanged.
hist = pd.read_csv(DATA / 'panel_b_histogram_geometry.csv')
assert int(hist['count'].sum()) == 200
b.bar(hist.xmin, hist['count'], width=hist.xmax-hist.xmin, align='edge',
      color=GREY, edgecolor='white', linewidth=.25)
b.axvline(.8, color='#222222', lw=.7, ls=(0, (4, 3)))
hits = confounder[confounder.role == 'hit']
b.scatter(hits.success, np.zeros(len(hits)), color=RED, s=16, zorder=3)
b.text(.98, 35, '\n'.join(hits.symbol), ha='right', va='top', color=RED,
       fontsize=6.5, linespacing=1.5)
b.set_xlim(-.05, 1.05)
b.set_ylim(-1.7, 37)
b.set_xticks([0, .25, .5, .75, 1], ['0.00', '0.25', '0.50', '0.75', '1.00'])
b.set_yticks([0, 10, 20, 30])
b.set_xlabel('Share of subsamples won by real Y status', fontsize=8, labelpad=6)
b.set_ylabel('Random non-hit genes', labelpad=3)

datasets = ['CRISPR dependencies', 'REP.PRIMARY', 'REP.1M', 'REP.300',
            'PRISM secondary', 'GDSC1', 'GDSC2', 'CTRP v2']
frame = selected.set_index('dataset_label').loc[datasets]
ys = np.arange(len(frame))
for column, shift, colour, label in [
    ('stable', .175, '#222222', 'stable in both seeds'),
    ('mean_null_confirmed', -.175, GREY, 'confirmed with shuffled Y (mean)')]:
    c.barh(ys+shift, frame[column], height=.325, color=colour, label=label)
    for y, value in zip(ys+shift, frame[column]):
        c.text(value+.035, y, f'{value:.1f}'.removesuffix('.0'), va='center', fontsize=6.2)
c.set_yticks(ys, datasets)
c.set_ylim(len(frame)-.4, -.6)
c.set_xlim(0, 5)
c.set_xticks(range(6))
c.set_xlabel('Features confirmed', labelpad=6)
c.legend(loc='upper left', bbox_to_anchor=(.05, -.25), frameon=False, fontsize=6.5,
         handlelength=1.2, handleheight=1.2, labelspacing=.8, borderaxespad=0)

stem = 'ch03_supp_depmap_random_forests'
for ext in ['pdf', 'svg', 'png']:
    fig.savefig(FIG / f'{stem}.{ext}', dpi=300, facecolor='white')
plt.close(fig)
(DATA / 'provenance.json').write_text(json.dumps({
    'source_tables': {name + '.csv': hashlib.sha256((DATA / f'{name}.csv').read_bytes()).hexdigest()
                      for name in names},
    'histogram_geometry': 'Original ggplot2 bin counts and edges, extracted without refitting',
    'font': font, 'figure_inches': [9.2, 3.2], 'statistics_refitted': False,
    'removed': ['panel titles', 'on-figure method footnotes'],
    'formal_document_caption_retained': True,
}, indent=2), encoding='utf-8')
print(f'Saved {FIG / (stem + ".pdf")}')
