"""Combine selected A549 F01 and F03 using their frozen numerical inputs."""
from pathlib import Path
import hashlib
import json
import shutil
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / 'figures/chapter3'
NUCLEAR = OUT / 'source_data/a549_primary_nuclear_RINe'
SOURCE = OUT / 'source_data/a549_rna_thesis'
DATA = OUT / 'source_data/a549_global_pathways_combined'
sys.path.insert(0, str(ROOT / 'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
from a549_matched_gsea import THEMES

BLUE, RED, TEAL, BROWN = '#3C5488', '#E64B35', '#00A087', '#7E6148'
STEM = 'ch03_main_a549_rna_global_structure'
GENES = ['KDM5D', 'RPS4Y1', 'CALB2', 'PDE3A', 'FGA', 'KRT19', 'PLCXD1', 'CSF2RA']


def render():
    font = apply_style(8)
    DATA.mkdir(parents=True, exist_ok=True)
    inputs = {
        'DE_nuclear_RINe_adjusted.csv': NUCLEAR / 'DE_nuclear_RINe_adjusted.csv',
        'DE_nuclear_unadjusted.csv': NUCLEAR / 'DE_nuclear_unadjusted.csv',
        'PCA_core.csv': SOURCE / 'PCA_core.csv',
        'PCA_variance_core.csv': SOURCE / 'PCA_variance_core.csv',
        'main_pathway_display.csv': SOURCE / 'main_pathway_display.csv',
        'camera_all_collections.csv': SOURCE / 'camera_all_collections.csv',
    }
    for name, source in inputs.items():
        shutil.copy2(source, DATA / name)
    de = pd.read_csv(inputs['DE_nuclear_RINe_adjusted.csv'])
    un = pd.read_csv(inputs['DE_nuclear_unadjusted.csv'])
    pca = pd.read_csv(inputs['PCA_core.csv'])
    variance = pd.read_csv(inputs['PCA_variance_core.csv'])
    pathways = pd.read_csv(inputs['main_pathway_display.csv']).set_index('pathway')
    camera = pd.read_csv(inputs['camera_all_collections.csv']).set_index('pathway')
    assert len(de) == 22843 and (de.FDR < .05).sum() == 27
    assert len(pca) == 8 and len(pathways) == 17
    assert pathways.NES.abs().max() <= 2.6

    fig = plt.figure(figsize=(7.2, 7.6))
    ax_a = fig.add_axes([.10, .665, .34, .285])
    ax_b = fig.add_axes([.565, .665, .405, .285])
    for group, colour in [('CTRL', BLUE), ('LOY', RED)]:
        group_data = pca[pca.treatment.eq(group)]
        ax_a.scatter(group_data.PC1, group_data.PC2, c=colour, s=33,
                     edgecolor='white', lw=.4)
        for row in group_data.itertuples():
            ax_a.annotate(str(row.sample_index), (row.PC1, row.PC2),
                          xytext=(5, 4), textcoords='offset points', fontsize=7, color=colour)
    ax_a.set(xlabel=f'PC1 ({100 * variance.core_variance.iloc[0]:.1f}%)',
             ylabel=f'PC2 ({100 * variance.core_variance.iloc[1]:.1f}%)')
    ax_a.margins(.18)
    ax_a.text(.02, .02, '8 libraries; 22,792 genes', transform=ax_a.transAxes,
              fontsize=6.8, color='#555555')
    ax_a.legend(handles=[Line2D([], [], marker='o', ls='', color=c, label=l)
                         for c, l in [(BLUE, 'sgCTRL'), (RED, 'Y-depleted')]],
                loc='lower left', bbox_to_anchor=(-.02, 1.015), ncol=2,
                frameon=False, fontsize=7, handletextpad=.3, columnspacing=.8, borderaxespad=0)

    y_gene = de.chrom.eq('chrY') & ~de.PAR
    other = ~y_gene & ~de.PAR
    significant = de.FDR < .05
    layers = [(other & ~significant, '#C8C8C8', 3),
              (other & significant & (de.logFC < 0), BLUE, 13),
              (other & significant & (de.logFC > 0), RED, 13),
              (de.PAR, TEAL, 17), (y_gene, BROWN, 18)]
    for mask, colour, size in layers:
        d = de[mask]
        ax_b.scatter(d.logFC, -np.log10(d.PValue), s=size, c=colour,
                     alpha=.6 if size == 3 else 1, linewidths=.25,
                     edgecolors='white' if size > 3 else 'none',
                     rasterized=size == 3, zorder=2 if size == 3 else 4)
    xlim = (min(-10, np.floor(min(de.logFC.min(), un.logFC.min())) - 1),
            max(5, np.ceil(max(de.logFC.max(), un.logFC.max())) + 1))
    ylim = (-.15, max(-np.log10(de.PValue.min()), -np.log10(un.PValue.min())) + 1.1)
    ax_b.set(xlim=xlim, ylim=ylim, xlabel='log₂FC (Y-depleted / sgCTRL)',
             ylabel='−log₁₀(nominal P)')
    for x in [-1, 1]:
        ax_b.axvline(x, c='#AAAAAA', ls=':', lw=.5)
    cutoff = de.loc[significant, 'PValue'].max()
    ax_b.axhline(-np.log10(cutoff), c='#666666', ls='--', lw=.6)
    label_positions = {
        'KDM5D': (-10.0, 5.10), 'RPS4Y1': (-5.2, 10.65),
        'CALB2': (-3.7, 5.85), 'PDE3A': (2.6, 5.65),
        'FGA': (3.1, 4.80), 'KRT19': (-3.9, 4.75),
        'PLCXD1': (.05, 3.80), 'CSF2RA': (-3.7, 2.80),
    }
    for gene in GENES:
        row = de[de.gene_name.eq(gene)].iloc[0]
        ax_b.annotate(gene, (row.logFC, -np.log10(row.PValue)),
                      xytext=label_positions[gene], textcoords='data', fontsize=6.5,
                      arrowprops=dict(arrowstyle='-', color='#888888', lw=.45), zorder=6)
    ax_b.text(.98, .98, '27 significant\n10 non-Y/non-PAR', transform=ax_b.transAxes,
              ha='right', va='top', fontsize=7)
    ax_b.legend(handles=[Line2D([], [], marker='o', ls='', color=c, label=l)
                         for c, l in [(BROWN, 'Y non-PAR'), (TEAL, 'PAR'),
                                      (BLUE, 'Other: lower'), (RED, 'Other: higher')]],
                loc='lower left', bbox_to_anchor=(-.015, 1.015), ncol=2,
                fontsize=6.6, frameon=False, handletextpad=.3, columnspacing=.8, borderaxespad=0)

    # One integrated pathway panel preserves all rows, themes, NES and both FDRs.
    rows, headings, y_position = [], [], 0
    for theme, keys in THEMES:
        keys = [key for key in keys if key in pathways.index]
        headings.append((y_position, theme))
        y_position += .9
        for key in keys:
            rows.append((y_position, key))
            y_position += 1
    assert [key for _, key in rows] == list(pathways.index)
    ax_c = fig.add_axes([.10, .085, .87, .465])
    ax_c.set(xlim=(0, 1), ylim=(y_position - .35, -1.45))
    ax_c.axis('off')
    cmap = LinearSegmentedColormap.from_list('NES', [BLUE, '#FAFAFA', RED])
    norm = TwoSlopeNorm(vmin=-2.6, vcenter=0, vmax=2.6)
    tile_x, tile_w, fg_x, camera_x = .535, .12, .765, .935
    for y, theme in headings:
        ax_c.text(0, y, theme, va='center', fontsize=7.6, weight='semibold', color='#414141')
    plotted = []
    for y, key in rows:
        row = pathways.loc[key]
        nes = float(row.NES)
        fg_fdr = float(row.BH_across_collections)
        camera_fdr = float(camera.loc[key, 'BH_across_collections'])
        ax_c.text(.012, y, row.display_label, va='center', fontsize=7.25)
        ax_c.add_patch(Rectangle((tile_x, y - .48), tile_w, .96,
                                 facecolor=cmap(norm(nes)), edgecolor='none'))
        # Exact scores in the tiles also make the colour scale directly readable.
        ax_c.text(tile_x + tile_w / 2, y, f'{nes:.2f}', va='center', ha='center',
                  fontsize=7.2, color='white' if abs(nes) > 1.65 else '#242424')
        for x, value in [(fg_x, fg_fdr), (camera_x, camera_fdr)]:
            ax_c.text(x, y, f'{value:.1e}', va='center', ha='center', fontsize=7.1,
                      weight='semibold' if value < .05 else 'normal', color='#333333')
        plotted.append({'pathway': key, 'display_label': row.display_label, 'NES': nes,
                        'fgsea_FDR': fg_fdr, 'CAMERA_FDR': camera_fdr})
    ax_c.text(tile_x + tile_w / 2, -.80, 'NES', ha='center', va='center', fontsize=7.4)
    ax_c.text(fg_x, -.80, 'fgsea FDR', ha='center', va='center', fontsize=7.4)
    ax_c.text(camera_x, -.80, 'CAMERA FDR', ha='center', va='center', fontsize=7.4)
    for letter, x, y in [('a', .025, .965), ('b', .505, .965), ('c', .025, .555)]:
        fig.text(x, y, letter, fontsize=12, weight='semibold')
    colour_axis = fig.add_axes([.535, .040, .32, .012])
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=colour_axis,
                      orientation='horizontal', ticks=[-2, -1, 0, 1, 2])
    cb.ax.tick_params(labelsize=6.5, length=2, pad=1)
    fig.text(.10, .046, 'H, Hallmark; R, Reactome; GO, GO Biological Process.', fontsize=6.4)
    fig.text(.10, .024, 'Bold FDR: <0.05. Pathway analyses use unadjusted expression.', fontsize=6.4)
    pd.DataFrame(plotted).to_csv(DATA / 'panel_c_display.csv', index=False)
    de.loc[de.gene_name.isin(GENES), ['gene_id', 'gene_name', 'logFC', 'PValue', 'FDR']].to_csv(
        DATA / 'panel_b_label_values.csv', index=False)
    exports = {}
    for ext in ['pdf', 'svg', 'png']:
        path = OUT / f'{STEM}.{ext}'
        fig.savefig(path, dpi=300, facecolor='white')
        exports[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
    plt.close(fig)
    provenance = {
        'scope': 'F01 PCA and adjusted volcano plus all 17 F03 pathways in panel c; no model refitting.',
        'font': font, 'size_inches': [7.2, 7.6],
        'inputs': {name: {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
                   for name, path in inputs.items()},
        'gene_tests': 22843, 'significant_genes': 27, 'pathways_displayed': 17,
        'pathway_testing_family': 23457, 'camera_inter_gene_correlation': .01,
        'volcano_limits': {'x': xlim, 'y': ylim}, 'exports': exports,
    }
    (DATA / 'provenance.json').write_text(json.dumps(provenance, indent=2), encoding='utf-8')
    print(f'Wrote {STEM}: 3 panels; 17 pathway rows; no statistical changes.')


if __name__ == '__main__':
    render()
