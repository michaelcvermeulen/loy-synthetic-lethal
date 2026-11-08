"""Draft comparison of saved original and CN-exclusion NES; leave S03 unchanged."""
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from du145_rna_common import FIG, T, MM, export

STEM = 'ch03_supp_du145_gsea_original_vs_cn_scatter'
INPUT = FIG / 'source_data/du145_rna_gene_set_themes/panel_a_all_significant_sets.csv'
THEMES = {
    'interferon_inflammatory': ('Interferon / inflammatory', '#00A087'),
    'GRHL2_epithelial': ('Epithelial', '#E39D69'),
    'mitotic': ('G2/M', '#B89139'),
    'ribosomal/translation': ('Ribosomal', '#4DBBD5'),
    'histone cluster (6p22)': ('Histone-associated', '#8491B4'),
    'other': ('Other', '#909090'),
    'positional': ('Chromosomal-position sets', '#7E6148'),
}


def main():
    data = pd.read_csv(INPUT)
    original = pd.read_csv(T / 'fgsea_all_collections.csv')
    extra = original[['pathway', 'BH_CN_stable_only', 'size_CN_stable_only']]
    data = data.merge(extra, on='pathway', validate='one_to_one')
    assert len(data) == 154 and data.survives_CN_screen.sum() == 139
    assert data.BH_across_collections.lt(.05).all()
    data['paired_NES_available'] = np.isfinite(data.NES) & np.isfinite(data.NES_CN_stable_only)
    data['NES_change'] = data.NES_CN_stable_only - data.NES
    paired = data[data.paired_NES_available].copy()
    missing = data[~data.paired_NES_available].copy()
    assert len(paired) == 150 and len(missing) == 4
    assert set(missing.pathway) == {'chr20p11', 'chr20p12', 'chr20p13', 'chr11q24'}
    assert (paired.survives_CN_screen == ((paired.BH_CN_stable_only < .05) &
            (np.sign(paired.NES) == np.sign(paired.NES_CN_stable_only)))).all()

    fig = plt.figure(figsize=(183 * MM, 120 * MM))
    ax = fig.add_axes([.095, .145, 101/183, 101/120])
    ax.set(xlim=(-2.5, 2.5), ylim=(-2.5, 2.5), aspect='equal',
           xlabel='Original normalised enrichment score (NES)',
           ylabel='NES after copy-number exclusion')
    ax.set_xticks([-2, -1, 0, 1, 2]); ax.set_yticks([-2, -1, 0, 1, 2])
    ax.axhline(0, color='#E1E1E1', lw=.6, zorder=0)
    ax.axvline(0, color='#E1E1E1', lw=.6, zorder=0)
    ax.plot([-2.5, 2.5], [-2.5, 2.5], color='#A0A0A0', ls='--', lw=.8, zorder=1)
    for theme, (label, colour) in THEMES.items():
        for retained in [True, False]:
            subset = paired[(paired.theme == theme) & (paired.survives_CN_screen == retained)]
            ax.scatter(subset.NES, subset.NES_CN_stable_only,
                       s=24 if retained else 39, facecolors=colour if retained else 'white',
                       edgecolors='white' if retained else colour, linewidths=.35 if retained else 1.1,
                       alpha=.9 if retained else 1, zorder=3 if retained else 4)
    # Label the sole paired score changing by more than 0.2 NES units.
    labelled = paired[paired.NES_change.abs() >= .2]
    assert list(labelled.pathway) == ['chr11q23']
    point = labelled.iloc[0]
    ax.annotate('chr11q23', (point.NES, point.NES_CN_stable_only),
                xytext=(-2.3, 1.90), fontsize=8, ha='left', va='center',
                arrowprops=dict(arrowstyle='-', color='#666666', lw=.7), color='#333333')
    ax.text(.04, .94, '150 paired scores', transform=ax.transAxes, fontsize=8, va='top')

    right = .70
    fig.text(right, .965, 'Gene-set theme', fontsize=8.4, weight='semibold', va='top')
    theme_handles = [Line2D([], [], marker='o', ls='', ms=5, markeredgewidth=0,
                           color=colour, label=label) for label, colour in THEMES.values()]
    fig.legend(handles=theme_handles, loc='upper left', bbox_to_anchor=(right-.015, .93),
               frameon=False, fontsize=7.5, handletextpad=.45, labelspacing=.65)
    status_handles = [
        Line2D([], [], marker='o', ls='', ms=5, color='#555555',
               label='Significant; same direction'),
        Line2D([], [], marker='o', ls='', ms=5, markerfacecolor='white', markeredgecolor='#555555',
               label='Not retained'),
        Line2D([], [], color='#A0A0A0', ls='--', lw=.8, label='Unchanged NES'),
    ]
    fig.text(right, .60, 'After CN exclusion', fontsize=8.4, weight='semibold', va='top')
    fig.legend(handles=status_handles, loc='upper left', bbox_to_anchor=(right-.015, .568),
               frameon=False, fontsize=7.2, handletextpad=.45, labelspacing=.65)
    fig.text(right, .40, '139 / 154 retain significance\nand enrichment direction',
             fontsize=8.2, weight='semibold', va='top', linespacing=1.5)
    fig.text(right, .285, 'Filtered NES unavailable (4):\nchr20p11, chr20p12,\nchr20p13, chr11q24',
             fontsize=7.6, color='#555555', va='top', linespacing=1.5)
    fig.text(.095, .025, 'CN exclusion removes genes in copy-number-altered intervals.',
             fontsize=7.5, color='#555555', va='bottom')

    export(fig, STEM,
           {'all_originally_significant_sets': data, 'plotted_paired_scores': paired,
            'filtered_NES_unavailable': missing, 'labelled_sets': labelled},
           {'inputs': [INPUT, T / 'fgsea_all_collections.csv'],
            'status': 'Draft alternative for S03a; selected S03 unchanged',
            'created_date': '2026-10-05',
            'notes': ('Show 150 paired NES values from the 154 originally significant sets. '
                      'Four sets without a filtered NES are listed separately and not plotted as zero. '
                      'Colour is the existing leading-edge theme. Filled markers retain FDR <0.05 '
                      'and NES direction after CN exclusion; open markers do not. The sole absolute '
                      'NES change >=0.2 is labelled. Equal linear axes and an identity line are used. '
                      'No jitter, correlation test, enrichment rerun or gene-set selection change.'),
            'summary': {'original_sets': 154, 'paired_scores': 150, 'unavailable': 4,
                        'retained_significance_and_direction': 139,
                        'paired_absolute_NES_change_below_0_2': int((paired.NES_change.abs() < .2).sum())},
            'theme_colours': THEMES})


if __name__ == '__main__':
    main()
