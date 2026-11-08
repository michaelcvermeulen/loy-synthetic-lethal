"""DU145 RNA supplementary figures: per-pair gene detail, and copy-number / model sensitivity."""
import numpy as np
import pandas as pd
import argparse
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from adjustText import adjust_text
import matplotlib.patheffects as pe
from du145_rna_common import (read, T, FOLLOW, COND, COND_LABEL, PASSAGE_MARKER, DARK, GREY, BLUE, RED, MM, panel, export)

from du145_reporting import annotate, LABEL
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--gene-detail-only', action='store_true', help='Export S01 without regenerating S02.')
args = parser.parse_args()
de = annotate(read('DE_paired_voom_primary'))
eff = read('DU145_LOY_vs_sgCTRL_primary_effects')
unp = read('DE_unpaired_voom_sensitivity')[['gene_id', 'logFC', 'FDR']].rename(columns={'logFC': 'logFC_unpaired', 'FDR': 'FDR_unpaired'})
cn = read('CN_gene_level')
win = read('CN_window_10Mb_RNA_vs_DNA')
arm = read('CN_arm_level_RNA_vs_DNA')
byarm = read('CN_FDR05_genes_by_arm')
cpm = read('core_CPM_all_genes')
PAIRS = ['pair_P12', 'pair_P13', 'pair_P14']

# ================= Supplementary 1: per-pair detail for the thesis gene set =================
ts = de[de.reporting_subset].merge(eff[['gene_id'] + PAIRS], on='gene_id')
ts = ts.sort_values('logFC')
from pathlib import Path
FIG_SOURCE=Path(__file__).resolve().parents[1]/'source_data/du145_rna_gene_detail'
FIG_SOURCE.mkdir(exist_ok=True,parents=True)
fig = plt.figure(figsize=(183 * MM, 250 * MM))
# Panel a: row-scaled expression of the complete reporting subset in all six libraries.
from matplotlib.colors import LinearSegmentedColormap
logcpm = read('core_logCPM_all_genes')
sample_order=[f'LNCaP-{condition}-P{passage}' for condition in ['sgCTRL','LOY'] for passage in [12,13,14]]
raw=logcpm.pivot(index='gene_id',columns='sample_key',values='log2_CPM').loc[ts.gene_id,sample_order]
assert raw.shape==(114,6) and not raw.isna().any().any()
row_sd=raw.std(axis=1,ddof=1)
assert (row_sd>0).all()
z=raw.sub(raw.mean(axis=1),axis=0).div(row_sd,axis=0)
assert np.allclose(z.mean(axis=1),0,atol=1e-10)
heat_axes=[fig.add_axes([left,.36,.18,.56]) for left in [.12,.44,.76]]
b=fig.add_axes([.12,.055,.80,.21])
cmap=LinearSegmentedColormap.from_list('expression_z',[BLUE,'#F7F7F7',RED])
for ax,indices in zip(heat_axes,np.array_split(np.arange(len(ts)),3)):
    chunk=ts.iloc[indices]
    im=ax.imshow(z.iloc[indices],aspect='auto',cmap=cmap,vmin=-2,vmax=2,interpolation='nearest')
    ax.set_yticks(range(len(chunk)),chunk.gene_name,fontsize=6)
    ax.set_xticks(range(6),['C12','C13','C14','Y12','Y13','Y14'],rotation=90,fontsize=6.5)
    for tick,col in zip(ax.get_xticklabels(),[BLUE]*3+[RED]*3):tick.set_color(col)
    ax.tick_params(length=0)
    ax.axvline(2.5,color='white',lw=1.2)
    for spine in ax.spines.values():spine.set_visible(False)
heat_axes[0].text(-.45,1.035,'a',transform=heat_axes[0].transAxes,weight='semibold',fontsize=12)
cax=fig.add_axes([.44,.948,.18,.009])
cb=fig.colorbar(im,cax=cax,orientation='horizontal',ticks=[-2,0,2]);cb.ax.tick_params(labelsize=6)
cb.ax.set_title('Row z-score',fontsize=7,pad=3)
fig.text(.5,.31,'C, sgCTRL; Y, Y-depleted; numbers indicate passage',ha='center',fontsize=7)
raw.to_csv(FIG_SOURCE / 'heatmap_log2CPM.csv')
z.to_csv(FIG_SOURCE / 'heatmap_row_z.csv')

# b: abundance versus effect size, showing the fold-change rule favours low-abundance genes.
sig = de.FDR < .05
thesis = de.reporting_subset
b.scatter(de.mean_CPM[~sig], de.logFC[~sig], s=3, color='#DDDDDD', lw=0, rasterized=True)
b.scatter(de.mean_CPM[sig & ~thesis], de.logFC[sig & ~thesis], s=10, color='#F0B27A', lw=0, label='FDR < 0.05')
b.scatter(de.mean_CPM[thesis], de.logFC[thesis], s=8, color=RED, lw=0, label=LABEL)
b.set_xscale('symlog', linthresh=.1)
for yv in (-.5, .5):
    b.axhline(yv, ls='--', lw=.6, color='#CCCCCC')
b.axvline(5, ls=':', lw=.7, color=BLUE)
b.text(4.3, 2.2, '5 CPM', fontsize=6.2, color=BLUE, ha='right')
lowc = int((ts.mean_CPM < 5).sum())
b.text(.97, .04, f'{lowc} of {len(ts)} selected genes below 5 CPM', transform=b.transAxes, ha='right', fontsize=6.4, color=DARK)
key = ts[ts.gene_name.isin(['GRHL2', 'CCR9', 'TACSTD2', 'SFTPA1', 'TAF9', 'THBD'])]
texts = [b.text(r.mean_CPM, r.logFC, r.gene_name, fontsize=6.3, color=DARK) for r in key.itertuples()]
adjust_text(texts, ax=b, only_move={'text': 'xy'}, arrowprops=dict(arrowstyle='-', color='#BBBBBB', lw=.45))
b.set_xlabel('Mean CPM (log scale)'); b.set_ylabel('log$_2$FC (Y-depleted / sgCTRL)')
b.legend(fontsize=6.3, loc='upper right', borderaxespad=.3)
panel(fig, b, 'b', '', dx=-.085)

export(fig, 'ch03_supp_du145_rna_gene_detail',
       {'panel_a_thesis_genes_by_pair': ts[['gene_name', 'chrom', 'logFC', 'FDR', 'mean_CPM'] + PAIRS],
        'panel_a_heatmap_row_z': z.reset_index(),
        'panel_a_heatmap_log2CPM': raw.reset_index(),
        'panel_b_abundance_effect': de[['gene_id', 'gene_name', 'mean_CPM', 'logFC', 'FDR', 'reporting_subset']]},
       {'inputs': [T / 'DE_paired_voom_primary.csv', T / 'DU145_LOY_vs_sgCTRL_primary_effects.csv', T / 'core_logCPM_all_genes.csv'],
        'notes': ('Panel a: frozen TMM log2 CPM (prior count 2), row z-scores across all six libraries using sample SD. Rows sorted by model log2FC, split into three consecutive blocks. Columns group sgCTRL P12/P13/P14 followed by LOY P12/P13/P14, with a divider between conditions. Colours clipped at +/-2. No passage correction of expression. '
                  'CCR9 is retained here for completeness although it is excluded from interpretation (10-16 reads per sgCTRL '
                  'library, none in Y-depleted libraries). Panel b reports abundance against effect size for all genes. '
                  'Highlighted genes meet FDR <0.05, |log2FC| >=0.5 and passage-direction agreement.')})

if args.gene_detail_only:
    raise SystemExit(0)

# S02 now includes the former F01 supporting panels; retain the current combined layout.
from make_du145_global_expression_combined import draw_supp
draw_supp()
