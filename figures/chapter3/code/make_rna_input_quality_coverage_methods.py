"""Compact methods version of the read-quality, input-RNA and coverage QC figure.

Reads the frozen plotted tables written by make_rna_quality_coverage.py; no
metric is recalculated. Per-library panels share one x axis and one row of
sample labels, with gene-body coverage as a final row of small multiples.
Styled with the thesis-methods-figure skill.
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
matplotlib.rcParams['figure.constrained_layout.h_pad'] = .2 / 72
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
from matplotlib.lines import Line2D

OUT = ROOT / 'figures/chapter3'
SRC = OUT / 'source_data/rna_quality_coverage'
d = pd.read_csv(SRC / 'plotted_metrics_with_provider.tsv', sep='\t').sort_values('plot_x').reset_index(drop=True)
cov = pd.read_csv(SRC / 'plotted_coverage.tsv', sep='\t')
assert len(d) == 44 and d.qc_id.is_unique and d.plot_x.is_unique
assert set(cov.qc_id) == set(d.qc_id) and (cov.groupby('qc_id').pos.nunique() == 100).all()
for col in ['q30_raw', 'q30_trim', 'adapter_trimmed_pct', 'fastp_pass_pct']:
    assert d[col].between(0, 100).all()
assert d.integrity.notna().all() and d.size_percent.notna().all()

colours = {'PARENTAL': '#8C8C8C', 'CTRL': '#3C5488', 'sgCTRL': '#3C5488', 'LOY': '#E64B35',
           'LOY_LXHA': '#C77C52', 'LOY_UTY_AB': '#00A087', 'LOY_KDM5D_AB': '#8491B4', 'LOY_KDM5D_UTY_AB': '#7E6148'}
legend_labels = {'PARENTAL': 'Parental', 'sgCTRL': 'sgCTRL', 'LOY': 'Y-depleted', 'LOY_LXHA': 'Empty vector',
                 'LOY_UTY_AB': 'UTY', 'LOY_KDM5D_AB': 'KDM5D', 'LOY_KDM5D_UTY_AB': 'UTY + KDM5D'}
groups = {'A549': 'A549 · Aug 2025', 'HT1080': 'HT1080 · Aug 2025',
          'HT1080 ADDBACK': 'HT1080 · May 2026', 'DU145': 'DU145 · May 2026'}
integrity_label = {'A549': 'RINe', 'HT1080': 'RINe', 'HT1080 ADDBACK': 'RQN', 'DU145': 'RIN (reported)'}
size_label = {'A549': 'DV200', 'HT1080': 'DV200', 'HT1080 ADDBACK': '200–10,000 nt', 'DU145': 'DV200'}
shaded = ['HT1080', 'DU145']
x = d.plot_x.to_numpy()
point_colours = d.display_condition.map(colours)
assert point_colours.notna().all()
spans = {g: (t.plot_x.min() - .5, t.plot_x.max() + .5) for g, t in d.groupby('group')}

fig = plt.figure(figsize=(6.5, 7.0))
top, bottom = fig.subfigures(2, 1, height_ratios=[5.0, 1.5])
axs = top.subplots(5, 1, sharex=True)
for ax in axs:
    for g in shaded:
        ax.axvspan(*spans[g], color='#F3F4F5', zorder=0, lw=0)
    ax.spines['bottom'].set_visible(False)
    ax.tick_params(axis='x', length=0)
    ax.xaxis.grid(False)
    ax.set_xlim(-1, x[-1] + 1)
a, b, c, e_int, f_size = axs


def band_notes(ax, labels):
    """Small grey notes centred over each group, in the band above the panel."""
    shift = mtransforms.ScaledTranslation(0, rcParams_band() / 2 / 72, fig.dpi_scale_trans)
    trans = mtransforms.blended_transform_factory(ax.transData, ax.transAxes) + shift
    for g, text in labels.items():
        ax.text(sum(spans[g]) / 2, 1, text, transform=trans, ha='center', va='center',
                fontsize=6, color='#555555', clip_on=False)


def rcParams_band():
    return matplotlib.rcParams['font.size'] * 1.7


before = a.bar(x - .19, d.q30_raw, width=.36, color='#999999', zorder=3)
after = a.bar(x + .19, d.q30_trim, width=.36, color='#3C5488', zorder=3)
a.set_ylabel('Q30 bases (%)')
a.set_ylim(0, 100)
b.legend([before, after], ['Before filtering', 'After filtering'], loc='lower right', bbox_to_anchor=(1, 1), ncol=2, borderaxespad=0)
for g, label in groups.items():
    ms.strip(a, label, span=(spans[g][0] - .35, spans[g][1] + .35))

b.bar(x, d.adapter_trimmed_pct, color=point_colours, width=.78)
b.set_ylabel('Trimmed (%)')
b.set_ylim(0, 35)

c.bar(x, d.fastp_pass_pct, color=point_colours, width=.78, zorder=3)
c.set_ylabel('Retained (%)')
c.set_ylim(0, 100)

e_int.bar(x, d.integrity, color=point_colours, width=.78, zorder=3)
e_int.set_ylabel('RNA integrity')
e_int.set_ylim(0, 10.5)
band_notes(e_int, integrity_label)

f_size.bar(x, d.size_percent, color=point_colours, width=.78, zorder=3)
f_size.set_ylabel('RNA (%)')
f_size.set_ylim(0, 100)
band_notes(f_size, size_label)
f_size.set_xticks(x, d.display_label, rotation=90)

cov_axes = bottom.subplots(1, 4, sharey=True)
ymax = np.ceil(cov.cov_rel.max() * 10) / 10
for ax, (g, label) in zip(cov_axes, groups.items()):
    for r in d.loc[d['group'] == g].itertuples():
        t = cov.loc[cov.qc_id == r.qc_id].sort_values('pos')
        ax.plot(t.pos, t.cov_rel, color=colours[r.display_condition], lw=.6, alpha=.75)
    ms.reference_line(ax, 1)
    ax.set_xlim(0, 99)
    ax.set_ylim(0, ymax)
    ax.set_xticks([0, 49.5, 99], ['5′', 'Mid', '3′'])
    for tick, ha in zip(ax.get_xticklabels(), ['left', 'center', 'right']):
        tick.set_ha(ha)  # keep 3′ and the next panel's 5′ apart
    ax.xaxis.grid(False)
    ms.strip(ax, label)
    ax.label_outer()
cov_axes[0].set_ylabel('Coverage / mean')

handles = [Line2D([], [], marker='s', ls='', markersize=4.5, color=colours[k]) for k in legend_labels]
fig.legend(handles, list(legend_labels.values()), loc='outside lower center', ncol=len(legend_labels),
           handletextpad=.2, columnspacing=1)
top.align_ylabels(axs)
ms.label_panels(list(axs) + [cov_axes[0]])
stem = 'ch03_methods_rna_input_quality_coverage'
ms.save_methods_figure(fig, stem, OUT)
print('Rendered', len(d), 'libraries and', len(cov), 'coverage values:', OUT / stem)
