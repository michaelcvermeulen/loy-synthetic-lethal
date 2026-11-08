"""Y/top-hit expression display and the unchanged RNA correlation panel."""
import hashlib
import json
import shutil
import numpy as np
import pandas as pd
from matplotlib.colors import to_rgba
from scipy.cluster.hierarchy import dendrogram, linkage, leaves_list
from scipy.spatial.distance import pdist

COL = {'CTRL': '#3C5488', 'LOY': '#E64B35', 'PARENTAL': '#888888'}


def draw_y_heatmap(fig, figure_root):
    source = figure_root / 'source_data/ht1080_rna_expression_changes'
    path = source / 'panel_b_Y_sample_expression.csv'
    ydata = pd.read_csv(path)
    ydata['heatmap_group'] = 'Y-linked'
    hit_path = source / 'panel_b_hit_sample_expression.csv'
    hits = pd.read_csv(hit_path)
    data = pd.concat([ydata, hits], ignore_index=True)
    meta = data[['sample_key', 'short', 'treatment', 'batch', 'passage']].drop_duplicates()
    meta['group_order'] = meta.treatment.map({'PARENTAL': 0, 'CTRL': 1, 'LOY': 2})
    meta = meta.sort_values(['batch', 'group_order', 'passage'])
    ygenes = ydata[['gene_id', 'gene_name', 'start', 'heatmap_group']].drop_duplicates().sort_values('start')
    hgenes = hits[['gene_id', 'gene_name', 'start', 'heatmap_group', 'heatmap_rank', 'logFC', 'FDR']].drop_duplicates()
    hgenes['group_order'] = hgenes.heatmap_group.map({'Higher': 0, 'Lower': 1})
    hgenes = hgenes.sort_values(['group_order', 'heatmap_rank'])
    genes = pd.concat([ygenes, hgenes], ignore_index=True)
    matrix = data.pivot(index='gene_id', columns='sample_key', values='log2_CPM_plus1').loc[genes.gene_id, meta.sample_key]
    assert matrix.shape[1] == 15 and np.isfinite(matrix).all().all()
    # Cluster the complete displayed matrix; do not constrain by gene class or population.
    row_distance = pdist(matrix.to_numpy(), metric='correlation')
    column_distance = pdist(matrix.to_numpy().T, metric='correlation')
    assert np.isfinite(row_distance).all() and np.isfinite(column_distance).all()
    row_tree = linkage(row_distance, method='average', optimal_ordering=True)
    column_tree = linkage(column_distance, method='average', optimal_ordering=True)
    row_order, column_order = leaves_list(row_tree), leaves_list(column_tree)
    genes = genes.iloc[row_order].reset_index(drop=True)
    meta = meta.iloc[column_order].reset_index(drop=True)
    matrix = matrix.iloc[row_order, column_order]
    ax = fig.add_axes([.407, .635, .193, .255])
    maximum = float(np.ceil(matrix.to_numpy().max()))
    im = ax.imshow(matrix, cmap='magma', vmin=0, vmax=maximum, aspect='auto')
    ax.set_yticks(range(len(genes)), genes.gene_name, fontsize=6)
    ax.set_xticks(range(len(meta)), meta.short, fontsize=6.2, rotation=90)
    ax.yaxis.tick_right()
    ax.tick_params(length=0, pad=2)
    ax.spines[:].set_visible(False)
    group_colours = {'Y-linked': '#555555', 'Higher': COL['LOY'], 'Lower': COL['CTRL']}
    for tick, group in zip(ax.get_yticklabels(), genes.heatmap_group):
        tick.set_color(group_colours[group])
    row_ax = fig.add_axes([.357, .635, .044, .255])
    row_dendro = dendrogram(row_tree, orientation='left', no_labels=True, ax=row_ax,
                           color_threshold=0, above_threshold_color='#555555')
    row_ax.set_ylim(10 * len(genes), 0)
    row_ax.axis('off')
    column_ax = fig.add_axes([.407, .914, .193, .042])
    column_dendro = dendrogram(column_tree, no_labels=True, ax=column_ax,
                              color_threshold=0, above_threshold_color='#555555')
    column_ax.set_xlim(0, 10 * len(meta))
    column_ax.axis('off')
    for tree_ax in [row_ax, column_ax]:
        for collection in tree_ax.collections:
            collection.set_linewidth(.6)
    assert row_dendro['leaves'] == row_order.tolist()
    assert column_dendro['leaves'] == column_order.tolist()
    for y, values in [(.894, [to_rgba(COL[t]) for t in meta.treatment]),
                      (.904, [to_rgba('#DADADA' if b == 'MAY_2026' else '#EEEEEE') for b in meta.batch])]:
        strip = fig.add_axes([.407, y, .193, .007])
        strip.imshow(np.array(values).reshape(1, 15, 4), aspect='auto')
        strip.axis('off')
    cb = fig.colorbar(im, cax=fig.add_axes([.407, .575, .193, .012]), orientation='horizontal')
    cb.set_ticks(np.arange(0, maximum + .1, 2))
    cb.ax.tick_params(labelsize=6.5, length=2)
    cb.set_label('log₂(CPM + 1)', fontsize=7, labelpad=1)
    cb.outline.set_visible(False)
    fig.text(.355, .925, 'b', fontsize=12, weight='semibold')
    matrix.to_csv(source / 'panel_b_expression_heatmap_matrix.csv')
    genes.to_csv(source / 'panel_b_heatmap_gene_order.csv', index=False)
    meta.to_csv(source / 'panel_b_Y_sample_order.csv', index=False)
    for name, tree in [('gene', row_tree), ('sample', column_tree)]:
        pd.DataFrame(tree, columns=['left', 'right', 'distance', 'count']).to_csv(
            source / f'panel_b_{name}_linkage.csv', index=False)
    # Linkage leaf IDs refer to the initial, unclustered order in this table.
    pd.DataFrame({'leaf': range(len(row_order)), 'gene_id': matrix.index[np.argsort(row_order)]}).to_csv(
        source / 'panel_b_gene_linkage_input_order.csv', index=False)
    pd.DataFrame({'leaf': range(len(column_order)), 'sample_key': matrix.columns[np.argsort(column_order)]}).to_csv(
        source / 'panel_b_sample_linkage_input_order.csv', index=False)
    paths = [path, hit_path, source/'panel_b_hit_selection.csv',
             figure_root/'source_data/ht1080_recurrent_gene_expression/all15_context_normalization.csv',
             figure_root/'source_data/shared_sex_chromosome/HT1080_effects.csv']
    (source/'panel_b_expression_provenance.json').write_text(json.dumps({
        'date': '2026-10-09', 'gene_count': len(genes), 'sample_count': len(meta),
        'genes': genes.gene_name.tolist(), 'sample_order': meta.short.tolist(),
        'Y_selection': 'Named non-PAR chrY genes with mean sgCTRL CPM >=1 in either period',
        'hit_selection': 'Three lowest-FDR genes per direction among FDR <0.05, absolute log2FC >=1 and matching directions in both periods; ties by gene ID',
        'normalisation': 'Saved all-15 contextual TMM effective library sizes',
        'scale': 'log2(CPM + 1), common scale; no row scaling, no period adjustment',
        'colour_limits': [0, maximum],
        'clustering': {'axes': ['genes', 'samples'], 'distance': '1 - Pearson correlation',
                       'linkage': 'average', 'optimal_leaf_ordering': True,
                       'input': 'All 20 displayed genes x 15 samples, unadjusted log2(CPM + 1)',
                       'groups_constrained': False},
        'row_order': 'Hierarchical clustering leaves, top to bottom',
        'column_order': 'Hierarchical clustering leaves, left to right',
        'analysis_refitted': False,
        'inputs_sha256': {str(p.relative_to(figure_root)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
    }, indent=2)+'\n', encoding='utf-8')


def draw_saved_correlation(fig, figure_root):
    source = figure_root / 'source_data/ht1080_rna_global_structure'
    out = figure_root / 'source_data/ht1080_pca_parental_combined'
    files = ['correlation_cluster_order.csv', 'Fig01_design_QC_c.csv',
             'correlation_linkage.csv', 'Fig01_design_QC_d.csv']
    hashes = {}
    for name in files:
        p = source / name
        hashes[name] = hashlib.sha256(p.read_bytes()).hexdigest()
        shutil.copy2(p, out / name)
    order = pd.read_csv(source/'correlation_cluster_order.csv')['sample'].tolist()
    cr = pd.read_csv(source/'Fig01_design_QC_c.csv')
    matrix = cr.pivot(index='sample_2', columns='sample_1', values='r').loc[order, order]
    tree = pd.read_csv(source/'correlation_linkage.csv').to_numpy()
    design = pd.read_csv(source/'Fig01_design_QC_d.csv')
    meta = design.set_index('short').loc[order]
    ax = fig.add_axes([.755, .212, .21, .313])
    im = ax.imshow(matrix, cmap='magma', vmin=matrix.to_numpy().min(), vmax=1, aspect='auto')
    ax.set_xticks(range(12), order, rotation=90, fontsize=6.5)
    ax.set_yticks(range(12), order, fontsize=6.5)
    ax.yaxis.tick_right(); ax.tick_params(length=0, pad=2); ax.spines[:].set_visible(False)
    da = fig.add_axes([.696, .212, .034, .313])
    dendrogram(tree, orientation='left', no_labels=True, ax=da,
               color_threshold=0, above_threshold_color='#666666')
    da.set_ylim(0, 120); da.axis('off')
    for x, values in [(.734, [to_rgba('#DADADA' if q == 'MAY_2026' else '#EEEEEE') for q in meta.batch]),
                      (.745, [to_rgba(COL[q]) for q in meta.treatment])]:
        strip = fig.add_axes([x, .212, .007, .313])
        strip.imshow(np.array(values).reshape(12, 1, 4), aspect='auto'); strip.axis('off')
    cb = fig.colorbar(im, cax=fig.add_axes([.755, .144, .21, .012]), orientation='horizontal')
    cb.set_ticks([.94, .96, .98, 1]); cb.ax.tick_params(labelsize=6.5, length=2)
    cb.set_label('Pearson r', fontsize=7, labelpad=1); cb.outline.set_visible(False)
    fig.text(.685, .55, 'g', fontsize=12, weight='semibold')
    fig.text(.755, .55, 'Expression correlations', fontsize=8.5)
    (out/'correlation_panel_g_provenance.json').write_text(json.dumps({
        'date': '2026-10-09', 'source_panel': 'Former F01b', 'new_panel': 'S01g',
        'sample_count': 12, 'sample_order': order, 'correlations_and_linkage_unchanged': True,
        'input_sha256': hashes,
    }, indent=2)+'\n', encoding='utf-8')
