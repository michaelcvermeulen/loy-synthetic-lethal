"""ch03_supp_crossline_gene_sets: every shared and divergent gene set, with NES per line and screens."""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from crossline_common import (ct, CT, LINES, THEME_LABEL, BLUE, RED, DARK, MM, panel, export)

STEM = 'ch03_supp_crossline_gene_sets'
SIGNED = LinearSegmentedColormap.from_list('signed', [BLUE, '#F7F7F7', RED])
gs = ct('gene_sets_shared_divergent')
theme_order = ['interferon_inflammatory', 'GRHL2_epithelial', 'mitotic', 'ribosomal/translation', 'polycomb_targets',
               'histone cluster (6p22)', 'positional', 'other']
gs['theme_rank'] = gs.theme.map({t: i for i, t in enumerate(theme_order)})
gs['class_rank'] = gs['class'].map({'shared': 0, 'divergent': 1})
gs = gs.sort_values(['class_rank', 'theme_rank', 'NES_DU145']).reset_index(drop=True)


def short(p, n=52):
    return p if len(p) <= n else p[:n - 1] + '…'


n = len(gs)
fig = plt.figure(figsize=(183 * MM, 240 * MM))
TOP, BOT = 0.905, 0.04
h = fig.add_axes([0.40, BOT, 0.15, TOP - BOT])        # NES heat, three lines
p = fig.add_axes([0.565, BOT, 0.055, TOP - BOT])      # HT1080 period signs
t = fig.add_axes([0.63, BOT, 0.30, TOP - BOT])        # theme and flags
cax = fig.add_axes([0.94, 0.70, 0.012, 0.18])         # colour bar
nes = gs[[f'NES_{L}' for L in LINES]].to_numpy()
bh = gs[[f'BH_{L}' for L in LINES]].to_numpy()
im = h.imshow(nes, cmap=SIGNED, norm=TwoSlopeNorm(0, -2.6, 2.6), aspect='auto', interpolation='nearest')
for i in range(n):
    for j in range(3):
        if bh[i, j] < .05:
            h.text(j, i, '•', ha='center', va='center', fontsize=6, color='#222222')
h.set_xticks(range(3), LINES, fontsize=6.4, rotation=45, ha='left'); h.xaxis.tick_top()
h.set_yticks(range(n), [short(x) for x in gs.pathway], fontsize=4.9)
per = np.column_stack([np.sign(gs['HT1080_NES_AUG_2025']), np.sign(gs['HT1080_NES_MAY_2026'])])
p.imshow(per, cmap=SIGNED, norm=TwoSlopeNorm(0, -1.5, 1.5), aspect='auto', interpolation='nearest')
p.set_xticks([0, 1], ['HT1080 AUG', 'HT1080 MAY'], fontsize=5.8, rotation=45, ha='left'); p.xaxis.tick_top()
p.set_yticks([])
t.set_xlim(0, 1); t.set_ylim(n - .5, -.5); t.axis('off')
for i, r in gs.iterrows():
    tags = []
    if r.positionally_concentrated:
        tags.append('positional')
    if r.BH_DU145 < .05 and not bool(r.DU145_survives_CN_screen):
        tags.append('DU145 CN-sensitive')
    t.text(0, i, THEME_LABEL.get(r.theme, r.theme) + ('  [' + ', '.join(tags) + ']' if tags else ''), va='center', fontsize=5.2, color=DARK)
for ax in (h, p):
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
split = int((gs['class'] == 'shared').sum())
for ax in (h, p):
    ax.axhline(split - .5, color='#222222', lw=1.0)
t.axhline(split - .5, color='#222222', lw=1.0, xmax=.6)
fig.text(.02, .985, f"Above the line: {split} shared sets (BH < 0.05 in at least two lines, same sign in all three).", fontsize=6.8, weight='semibold')
fig.text(.02, .972, f"Below the line: {n - split} divergent sets (significant in at least two lines with opposite signs).", fontsize=6.8, weight='semibold')
cb = fig.colorbar(im, cax=cax); cb.set_label('NES', fontsize=6.5); cb.ax.tick_params(labelsize=6)
fig.text(.40, .012, 'Dots: BH < 0.05 in that line (all eligible sets). HT1080 columns: sign in each sampling period.', fontsize=5.8, color=DARK)

export(fig, STEM,
       {'shared_and_divergent_sets': gs.drop(columns=['theme_rank', 'class_rank'])},
       {'inputs': [CT / 'gene_sets_shared_divergent.csv'],
        'notes': ('All 46 shared and 23 divergent gene sets from the cross-line analysis, ordered by class, leading-edge theme and DU145 NES. '
                  'Names longer than 52 characters are truncated; full names, leading edges and statistics are in the source table.')})
