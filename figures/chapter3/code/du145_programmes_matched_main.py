"""Selected DU145 F04: compact pathways, leading edges, CN sensitivity and markers."""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from du145_rna_common import (read, T, FIG, THEMES, FOLLOW, BLUE, RED, DARK,
                              MM, export, THEME_ORDER)
from du145_cn_only_pathways import (cn_only_themes, THEME_LABEL, STATUS_ORDER,
                                   STATUS_COLOUR, STATUS_LABEL)


def render_main():
    stem = 'ch03_main_du145_rna_programmes_regulators'
    source = FIG/'source_data/du145_rna_pathways/pathway_display.csv'
    display = pd.read_csv(source)
    fg = read('fgsea_all_collections').set_index('pathway')
    cam = read('camera_all_collections').set_index('pathway')
    assert len(fg) == len(cam) == 23252 and len(display) == 19
    assert cam.correlation_mode.eq('fixed_0.01').all()
    for row in display.itertuples():
        assert np.isclose(row.NES, fg.loc[row.pathway, 'NES'], rtol=1e-12)
        assert np.isclose(row.fgsea_FDR, fg.loc[row.pathway, 'BH_across_collections'], rtol=1e-12, atol=1e-15)
        assert np.isclose(row.CAMERA_FDR, cam.loc[row.pathway, 'BH_across_collections'], rtol=1e-12, atol=1e-25)
    assert display.NES.abs().max() <= 2.6
    # Display selection only: exact frozen symbol/ranking values used by fgsea.
    ranking = read('pathway_TF_effect_ranking').set_index('gene_name', verify_integrity=True)
    leading_records = []
    for row in display.itertuples():
        genes = row.leadingEdge.split(';')
        assert len(genes) == len(set(genes)) and len(genes) >= 3
        ranked = ranking.loc[genes].reset_index()
        ranked['abs_effect'] = ranked.LOY_vs_sgCTRL.abs()
        ranked = ranked.sort_values(['abs_effect', 'gene_name'], ascending=[False, True], kind='stable')
        chosen = ranked.head(3)
        assert (np.sign(chosen.LOY_vs_sgCTRL) == np.sign(row.NES)).all()
        label = ', '.join(chosen.gene_name)
        display.loc[display.pathway == row.pathway, 'leading_genes_display'] = label
        for rank, gene in enumerate(ranked.itertuples(), 1):
            leading_records.append(dict(pathway=row.pathway, gene_name=gene.gene_name,
                                        gene_id=gene.gene_id, ranking_effect=gene.LOY_vs_sgCTRL,
                                        abs_effect=gene.abs_effect, within_leading_edge_rank=rank,
                                        displayed=rank <= 3))
    leading_table = pd.DataFrame(leading_records)
    assert leading_table.displayed.sum() == 57
    nrows = int(display.row.max()) + 1
    heat = np.full((nrows, 1), np.nan)
    heat[display.row.to_numpy(), 0] = display.NES
    sig = fg.BH_across_collections < .05

    width, height = 183, 190
    fig = plt.figure(figsize=(width*MM, height*MM))
    def axes(x, y, w, h):
        return fig.add_axes([x/width, y/height, w/width, h/height])
    def text(x, y, label, **kwargs):
        return fig.text(x/width, y/height, label, **kwargs)

    # a, b: narrow NES strip, row-matched leading genes and compact FDR columns.
    a = axes(77, 86, 8.5, 87)
    genes_ax = axes(90, 86, 52, 87)
    b = axes(146, 86, 33, 87)
    cmap = LinearSegmentedColormap.from_list('DU145_NES', [BLUE, '#FAFAFA', RED])
    im = a.imshow(heat, aspect='auto', cmap=cmap,
                  norm=TwoSlopeNorm(vmin=-2.6, vcenter=0, vmax=2.6), interpolation='nearest')
    a.set_yticks(display.row, display.display_label, fontsize=6.8)
    a.set_xticks([])
    a.tick_params(length=0)
    a.spines[:].set_visible(False)
    for theme, rows in display.groupby('theme', sort=False):
        a.text(-.63, rows.row.min()-1, theme, ha='right', va='center',
               fontsize=6.9, weight='semibold', color=BLUE if rows.NES.mean()<0 else '#B5422B',
               clip_on=False)
    genes_ax.set(xlim=(0, 1), ylim=(nrows-.5, -.5))
    genes_ax.axis('off')
    b.set(xlim=(-.5, 1.5), ylim=(nrows-.5, -.5))
    b.axis('off')
    for r in display.itertuples():
        genes_ax.text(0, r.row, r.leading_genes_display, ha='left', va='center',
                      fontsize=6.8, color='#333333')
        if r.fgsea_FDR < .05:
            a.text(0, r.row, '*', ha='center', va='center', fontsize=10,
                   color='white' if abs(r.NES)>2.6*.55 else '#242424')
        for col, val in enumerate([r.fgsea_FDR, r.CAMERA_FDR]):
            b.text(col, r.row, f'{val:.1e}', ha='center', va='center', fontsize=6.8,
                   weight='semibold' if val<.05 else 'normal', color=DARK)
    for col, label in enumerate(['fgsea', 'CAMERA']):
        b.text(col, -.95, label, ha='center', va='bottom', fontsize=7)
    text(4, 184, 'a', fontsize=12, weight='semibold')
    text(13, 184, f'{len(fg):,} sets tested; {sig.sum()} significant '
         f'({(sig & (fg.NES<0)).sum()} lower, {(sig & (fg.NES>0)).sum()} higher)', fontsize=7.4)
    text(145, 184, 'b', fontsize=12, weight='semibold')
    text(153, 184, 'Global FDR', fontsize=8.5, weight='semibold')
    text(81.25, 177, 'NES', fontsize=7.1, ha='center')
    text(90, 177, 'Leading-edge genes', fontsize=7.1)
    cax = axes(40, 76, 53, 2.1)
    cb = fig.colorbar(im, cax=cax, orientation='horizontal', ticks=[-2,-1,0,1,2])
    cb.ax.tick_params(labelsize=6.6, pad=1.5)
    cb.ax.text(1.07, .5, 'NES', transform=cb.ax.transAxes, fontsize=7, va='center')

    # c: all significant sets; proportions and screen definitions are unchanged.
    c = axes(30, 18, 56, 45)
    th = cn_only_themes()
    counts = th.groupby(['theme','status']).size().unstack(fill_value=0).reindex(
        index=THEME_ORDER, columns=STATUS_ORDER, fill_value=0)
    yc = np.arange(len(THEME_ORDER))[::-1]
    left = np.zeros(len(THEME_ORDER))
    for status in STATUS_ORDER:
        c.barh(yc, counts[status], left=left, color=STATUS_COLOUR[status],
               height=.68, label=STATUS_LABEL[status])
        left += counts[status].to_numpy()
    direction = th.groupby('theme').NES.agg(
        lambda v: 'all down' if (v<0).all() else ('all up' if (v>0).all() else 'mixed'))
    for y, theme, number in zip(yc, THEME_ORDER, left):
        c.text(number+.8, y, f'{int(number)}, {direction[theme]}',
               va='center', fontsize=6.1, color=DARK)
    c.set_yticks(yc, [THEME_LABEL[t] for t in THEME_ORDER], fontsize=6.6)
    c.set(xlim=(0,66), xlabel='Significant gene sets')
    c.set_xticks([0,20,40,60])
    c.tick_params(axis='x', labelsize=7)
    c.xaxis.label.set_fontsize(7.3)
    c.legend(loc='lower left', bbox_to_anchor=(-.02,1.03), ncol=1, fontsize=6.2,
             handlelength=.9, handletextpad=.4, labelspacing=.2, borderpad=0)
    text(4, 69, 'c', fontsize=12, weight='semibold')

    # d: the same ten marker genes and saved passage-pair differences.
    d = axes(112, 18, 49, 48)
    selected = ['TACSTD2','AP1M2','CLDN7','MAL2','GRHL2','BST2','CXCL1','LTB','CXCL10','B2M']
    effects, de = read('DU145_LOY_vs_sgCTRL_primary_effects'), read('DE_paired_voom_primary')
    mk = effects.set_index('gene_name').loc[selected].reset_index()
    mk = mk.merge(de[['gene_id','logFC','FDR']], on='gene_id', validate='one_to_one')
    assert mk.gene_name.tolist() == selected
    mk['marker_group'] = ['Epithelial']*5 + ['Interferon/inflammatory']*5
    mk['sig'] = mk.FDR < .05
    ys = np.arange(len(mk))[::-1].astype(float)
    ys[:5] += .6
    for passage, marker, offset in [('P12','o',.13),('P13','s',0),('P14','^',-.13)]:
        for i, row in enumerate(mk.itertuples()):
            d.scatter(getattr(row,'pair_'+passage), ys[i]+offset, marker=marker, s=17,
                       facecolors=RED if row.sig else 'white', edgecolors=RED, lw=.6, zorder=3)
        d.plot([], [], ls='', marker=marker, color=RED, ms=3.5, label=passage)
    for i, row in enumerate(mk.itertuples()):
        vals = [row.pair_P12,row.pair_P13,row.pair_P14]
        d.hlines(ys[i], min(vals), max(vals), color='#CCCCCC', lw=.6, zorder=1)
        d.text(1.24, ys[i], f'{row.CTRL_mean_CPM:.1f}', transform=d.get_yaxis_transform(),
               fontsize=6.5, ha='center', va='center', color=DARK)
    d.text(1.24,1.025,'sgCTRL\nCPM', transform=d.transAxes, fontsize=6.3,
           ha='center', va='bottom', color=DARK)
    d.axhline(4.8, color='#DDDDDD', lw=.6)
    d.axvline(0, color='#AAAAAA', lw=.7)
    d.set(ylim=(-.6,10.2), xlim=(-4.6,.9), xlabel='Passage-pair log$_2$CPM difference')
    d.set_xticks([-4,-2,0])
    d.set_yticks(ys, mk.gene_name, fontsize=6.8)
    d.tick_params(axis='x', labelsize=7)
    d.xaxis.label.set_fontsize(7.1)
    d.legend(loc='lower left', bbox_to_anchor=(-.04,1.045), ncol=3, fontsize=6.3,
             columnspacing=.7, handletextpad=.2, borderpad=0, handlelength=.8)
    text(97, 69, 'd', fontsize=12, weight='semibold')
    text(7, 5, '* fgsea FDR < 0.05. Bold: FDR < 0.05. Leading-edge genes: three largest absolute ranking effects per set.', fontsize=6.2)
    text(7, 1, 'H, Hallmark; R, Reactome; CGP, chemical/genetic perturbations; MF/CC, GO molecular function/cellular component.', fontsize=6.1)

    export(fig, stem, {
        'panel_a_b_pathways': display,
        'panel_a_leading_edge_memberships': leading_table,
        'panel_c_theme_status_counts': counts.reset_index(),
        'panel_c_all_significant_sets': th,
        'panel_c_collection_summary': read('enrichment_collection_summary'),
        'panel_d_markers': mk,
    }, {
        'inputs': [source, THEMES] + [T/f'{n}.csv' for n in [
            'fgsea_all_collections','camera_all_collections','enrichment_collection_summary',
            'DU145_LOY_vs_sgCTRL_primary_effects','DE_paired_voom_primary','pathway_TF_effect_ranking']],
        'notes': ('Selected F04 revision: compact 183 x 190 mm layout; panel a adds three leading-edge genes '
                  'per pathway, ranked by absolute original mean log2CPM difference with alphabetical '
                  'tie-breaking. The exact frozen pathway symbol-to-gene mapping is retained. '
                  'No gene-level significance or passage-consistency filter is applied to these labels. '
                  'a/b retain the 19-set A549-format DU145 pathway display, '
                  'c retains the original all-set CN sensitivity summary, and d retains the same '
                  'ten passage-specific marker genes. The former full Hallmark display is archived; '
                  'all four former epithelial bars are represented in the heatmap. '
                  'Effect-ranked fgsea and paired CAMERA with fixed correlation 0.01 use global '
                  'BH families of 23,252 sets. No analysis, threshold, gene selection or statistic '
                  'has been refitted. CN exclusion uses absolute lpWGS segment difference >=0.2; '
                  '139/154 sets retain significance and direction. The two nonsignificant fgsea '
                  'Hallmark cell-cycle sets are retained in a/b for context. '
                  'Historical panel_a_hallmark, panel_b_theme_status_counts and panel_c_epithelial_sets '
                  'tables in this folder correspond to the archived layout; current tables are '
                  'listed in README.md.'),
    })
    return display, counts, mk


if __name__ == '__main__':
    render_main()
