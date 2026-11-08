"""Compact methods version of the 44-library RNA-seq QC figure.

Reads the frozen plotted table written by make_rna_all_samples_qc.py; no
metric is recalculated and no analysis is rerun. Styled with the
thesis-methods-figure skill (7 pt, no titles, tight panels).
"""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'skills/thesis-methods-figure/scripts'))
import methods_style as ms
ms.apply_methods_style(grid=True)
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

OUT = ROOT / 'figures/chapter3'
SRC = OUT / 'source_data/rna_all_samples_qc'
d = pd.read_csv(SRC / 'plotted_metrics.tsv', sep='\t').sort_values('plot_x').reset_index(drop=True)
assert len(d) == 44 and d.qc_id.is_unique and d.plot_x.is_unique and d.display_condition.notna().all()

colours = {'PARENTAL': '#8C8C8C', 'CTRL': '#3C5488', 'sgCTRL': '#3C5488', 'LOY': '#E64B35',
           'LOY_LXHA': '#C77C52', 'LOY_UTY_AB': '#00A087', 'LOY_KDM5D_AB': '#8491B4', 'LOY_KDM5D_UTY_AB': '#7E6148'}
legend_labels = {'PARENTAL': 'Parental', 'sgCTRL': 'sgCTRL', 'LOY': 'Y-depleted', 'LOY_LXHA': 'Empty vector',
                 'LOY_UTY_AB': 'UTY', 'LOY_KDM5D_AB': 'KDM5D', 'LOY_KDM5D_UTY_AB': 'UTY + KDM5D'}
groups = {'A549': 'A549 · Aug 2025', 'HT1080': 'HT1080 · Aug 2025',
          'HT1080 ADDBACK': 'HT1080 · May 2026', 'DU145': 'DU145 · May 2026'}
shaded = ['HT1080', 'DU145']
x = d.plot_x.to_numpy()
bar_colours = d.display_condition.map(colours)
assert bar_colours.notna().all()
spans = {g: (t.plot_x.min() - .5, t.plot_x.max() + .5) for g, t in d.groupby('group')}

fig, axs = plt.subplots(6, 1, figsize=(6.5, 6.3), sharex=True, height_ratios=[1, 1, 1, 1, 1, 1])
for ax in axs:
    for g in shaded:
        ax.axvspan(*spans[g], color='#F3F4F5', zorder=0, lw=0)
    ax.spines['bottom'].set_visible(False)
    ax.tick_params(axis='x', length=0)
    ax.xaxis.grid(False)
    ax.set_xlim(-1, x[-1] + 1)

def band_legend(ax, handles, labels):
    """Legend in the band above the panel, right-aligned, level with the letter."""
    ax.legend(handles, labels, loc='lower right', bbox_to_anchor=(1, 1), ncol=len(labels),
              borderaxespad=0, frameon=False)

a, b, c, e_dup, f_mt, g_bias = axs
a.bar(x, d.star_input_pairs / 1e6, color=bar_colours, width=.78)
a.set_ylabel('Read pairs (M)')
a.set_ylim(0, 52)
for g, label in groups.items():
    ms.strip(a, label, span=(spans[g][0] - .35, spans[g][1] + .35))  # reach into the gaps so labels fit

bottom = np.zeros(len(d))
parts = [(d.star_unique_pct, 'Unique', '#3C5488'), (d.star_multi_pct, 'Multiple loci', '#91B8B5'),
         (100 - d.star_unique_pct - d.star_multi_pct, 'Unmapped', '#D5D5D5')]
for v, label, col in parts:
    assert np.isfinite(v).all() and (v >= -1e-5).all()
    b.bar(x, v, bottom=bottom, width=.78, color=col, label=label)
    bottom += v
b.set_ylabel('Pairs (%)')
b.set_ylim(0, 100)
band_legend(b, *b.get_legend_handles_labels())

bottom = np.zeros(len(d))
for key, label, col in [('qm_exonic_pct', 'Exonic', '#3C5488'), ('qm_intronic_pct', 'Intronic', '#E9A36F'),
                        ('qm_intergenic_pct', 'Intergenic', '#B8B3C9')]:
    c.bar(x, d[key], bottom=bottom, width=.78, color=col, label=label)
    bottom += d[key]
assert np.allclose(bottom, 100)
c.set_ylabel('Reads (%)')
c.set_ylim(0, 100)
band_legend(c, *c.get_legend_handles_labels())

e_dup.bar(x, d.dup_pct, color=bar_colours, width=.78)
e_dup.set_ylabel('Duplication (%)')
e_dup.set_ylim(0, 40)

f_mt.scatter(x - .15, d.chrM_pct, s=5, color='#3C5488', label='chrM alignments', zorder=3)
f_mt.scatter(x + .15, 100 * d.pic_ribosomal, s=6, marker='D', color='#E69F00', label='rRNA bases', zorder=3)
f_mt.set_ylabel('Abundance (%)')
f_mt.set_ylim(0, 20)
band_legend(f_mt, *f_mt.get_legend_handles_labels())

g_bias.scatter(x, d.qm_bias_53, s=7, c=bar_colours, zorder=3)
ms.reference_line(g_bias, 1)
g_bias.set_ylabel('5′/3′ coverage')
g_bias.set_ylim(.95, 1.4)
g_bias.set_xticks(x, d.display_label, rotation=90)

handles = [Line2D([], [], marker='s', ls='', markersize=4.5, color=colours[k]) for k in legend_labels]
fig.legend(handles, list(legend_labels.values()), loc='outside lower center', ncol=len(legend_labels),
           handletextpad=.2, columnspacing=1)
fig.align_ylabels(axs)
ms.label_panels(axs)
stem = 'ch03_methods_rna_library_qc'
ms.save_methods_figure(fig, stem, OUT)
print('Rendered', len(d), 'libraries:', OUT / stem)
