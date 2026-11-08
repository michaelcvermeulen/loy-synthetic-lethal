"""DU145 S03: significant gene-set themes and passage-specific leading-edge expression."""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from du145_leading_edge_genes import select_leading_edge_genes, BIO_THEMES
from du145_rna_common import (DARK, GREY, BLUE, RED, MM, panel, export, THEMES, THEME_ORDER, THEME_LABEL,
                              STATUS_ORDER, STATUS_COLOUR, STATUS_LABEL, FOLLOW)

from du145_cn_only_pathways import cn_only_themes, STATUS_ORDER, STATUS_COLOUR, STATUS_LABEL, THEME_LABEL

STEM = 'ch03_supp_du145_rna_gene_set_themes'
th = cn_only_themes()
rng = np.random.default_rng(7)

WIDTH, HEIGHT = 183, 158
fig = plt.figure(figsize=(WIDTH * MM, HEIGHT * MM))
def mm_axes(x, y, width, height):
    return fig.add_axes([x/WIDTH, y/HEIGHT, width/WIDTH, height/HEIGHT])
def mm_text(x, y, text, **kwargs):
    return fig.text(x/WIDTH, y/HEIGHT, text, **kwargs)
a = mm_axes(29.3, 103, 146.4, 43)

# ---- a: every significant set, positioned by NES within its theme.
ya = {t: i for i, t in enumerate(THEME_ORDER[::-1])}
for st in STATUS_ORDER[::-1]:
    d = th[th.status.eq(st)]
    a.scatter(d.NES, d.theme.map(ya) + rng.uniform(-.28, .28, len(d)), s=11, color=STATUS_COLOUR[st], lw=0, alpha=.9,
              label=STATUS_LABEL[st], zorder=3 if st == 'retained' else 2)
a.axvline(0, lw=.6, color='#BBBBBB')
a.set_yticks(list(ya.values()), [THEME_LABEL[t] for t in ya], fontsize=6.6)
a.set_ylim(-.6, len(ya) - .4); a.set_xlim(-2.9, 2.9)
a.set_xlabel('Normalised enrichment score')
for t, y_ in ya.items():
    d = th[th.theme.eq(t)]
    a.text(2.85, y_, f'{(d.status == "retained").sum()} of {len(d)}', ha='right', va='center', fontsize=6, color=DARK)
a.text(2.85, len(ya) - .45, 'CN retained / total', ha='right', va='bottom', fontsize=5.8, color=DARK)
handles = [Line2D([], [], ls='', marker='o', ms=4, color=STATUS_COLOUR[s], label=STATUS_LABEL[s]) for s in STATUS_ORDER]
a.legend(handles=handles, fontsize=6, loc='center', borderaxespad=.3, handletextpad=.2, labelspacing=.25)
panel(fig, a, 'a', '', dx=-.14)

# ---- b: recurrent leading-edge genes, with actual passage-specific effects.
selected_genes, candidates, memberships, gene_summary = select_leading_edge_genes(th)
assert len(selected_genes) == 30
cmap = LinearSegmentedColormap.from_list('DU145_gene_effect', [BLUE, '#FAFAFA', RED])
norm = TwoSlopeNorm(vmin=-1.5, vcenter=0, vmax=1.5)
positions = [(17,59),(78,59),(139,59),(47.5,21),(108.5,21)]
short_titles = ['Interferon', 'Epithelial', 'G2/M', 'Ribosomal', 'Histone-associated']
for theme, title, (x,y) in zip(BIO_THEMES, short_titles, positions):
    rows = selected_genes[selected_genes.theme.eq(theme)]
    matrix = rows[['pair_P12','pair_P13','pair_P14']].to_numpy()
    b = mm_axes(x,y,24,24)
    im = b.imshow(matrix, aspect='auto', cmap=cmap, norm=norm, interpolation='nearest')
    labels = [r.gene_name + ('*' if r.gene_FDR_below_005 else '') for r in rows.itertuples()]
    b.set_yticks(range(len(rows)), labels, fontsize=6.7)
    b.set_xticks(range(3), ['P12','P13','P14'], fontsize=6.4)
    b.xaxis.tick_top()
    b.tick_params(length=0, pad=2)
    b.spines[:].set_visible(False)
    for i in range(len(rows)):
        for j in range(3):
            value = matrix[i,j]
            label = f'{value:.1f}'
            if label == '-0.0':
                label = '0.0'
            b.text(j,i,label,ha='center',va='center',fontsize=6.3,
                   color='white' if abs(value)>.9 else '#333333')
    frequency = mm_axes(x+26.5,y,11.5,24)
    frequency.set(xlim=(0,1),ylim=(len(rows)-.5,-.5))
    frequency.axis('off')
    frequency.text(.5,-.85,'Sets',ha='center',va='center',fontsize=6.3,color=DARK)
    for i,r in enumerate(rows.itertuples()):
        frequency.text(.5,i,f'{r.leading_edge_sets}/{r.retained_sets}',ha='center',va='center',
                       fontsize=6.3,color=DARK)
    mm_text(x-14,y+31,f'{title} ({int(rows.retained_sets.iloc[0])} sets)',fontsize=7.2,
            weight='semibold',color=BLUE if rows.mean_difference.mean()<0 else '#B5422B')
mm_text(3.5,94,'b',fontsize=12,weight='semibold')
# The shared colour scale is capped for large negative effects; all printed effects are untruncated.
cbar_ax = mm_axes(50,9,66,2.3)
cb = fig.colorbar(im,cax=cbar_ax,orientation='horizontal',extend='min',ticks=[-1.5,-.75,0,.75,1.5])
cb.ax.tick_params(labelsize=6.4)
cb.ax.set_title('Y-depleted minus sgCTRL log$_2$CPM difference',fontsize=6.8,pad=5)
mm_text(122,10,'* Gene-level FDR < 0.05',fontsize=6.3,color=DARK)
mm_text(122,5,'Colours capped; values printed in full range',fontsize=5.8,color=DARK)

export(fig, STEM,
       {'panel_a_all_significant_sets': th.drop(columns=['leadingEdge']),
        'panel_b_selected_genes': selected_genes,
        'panel_b_all_candidate_genes': candidates,
        'panel_b_leading_edge_memberships': memberships,
        'panel_b_theme_summary': gene_summary,
        'all_significant_sets_with_leading_edges': th},
       {'inputs': [__import__('du145_rna_common').T / 'fgsea_all_collections.csv', THEMES, FOLLOW / 'gene_set_themes/classify_themes.R',
                   *[__import__('du145_rna_common').T / f'{name}.csv' for name in
                     ['pathway_TF_effect_ranking','DU145_LOY_vs_sgCTRL_primary_effects','DE_paired_voom_primary']]],
        'notes': ('CN-only sensitivity: retained means significant with matching NES direction after exclusion of CN-changed genes; no positional exclusion. Themes are assigned from leading-edge composition, not set names, using gene groups fixed in advance: '
                  'Hallmark IFN-gamma, IFN-alpha and inflammatory response; Hallmark G2/M checkpoint, mitotic spindle and '
                  'E2F; the DepMap GRHL2 top-200 co-expression neighbourhood; cytosolic and mitochondrial ribosomal proteins '
                  '(>= 40% of the leading edge); replication-dependent histones (>= 30%). Otherwise the largest reference '
                  'overlap of at least 20% assigns the theme. Full names, statistics and leading edges for all 154 sets are '
                  'in all_significant_sets_with_leading_edges.csv. Former panel c was removed at the author\'s request '
                  'on 9 October 2026; its saved target-link sensitivity data remain historical. Panel b shows six '
                  'genes per coherent biological theme, ranked by frequency in the original leading edges of '
                  'CN-retained sets, then absolute mean expression difference, then symbol. No gene-level '
                  'FDR or consistency filter is used for selection. Colours show saved log2CPM differences '
                  'for each passage pair, capped at +/-1.5; all numerical labels show the full values to '
                  'one decimal. An asterisk on the gene label denotes paired-model FDR <0.05 across '
                  '18,767 genes. Counts are gene-set membership frequencies, not independent replicates. '
                  'The exact frozen pathway symbol-to-gene mapping is retained. Other and positional '
                  'groups are omitted from panel b because they are not coherent biological themes.')})
