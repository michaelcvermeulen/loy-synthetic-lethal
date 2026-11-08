"""Extract saved results for an exploratory ZMYND8 review; no model refitting."""
from pathlib import Path
import hashlib
import json
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import LinearSegmentedColormap

sys.stdout.reconfigure(encoding='utf-8')
THESIS = Path(__file__).resolve().parents[3]
PROJECTS = THESIS.parent
FIG = THESIS / 'figures/chapter3'
SRC = FIG / 'source_data/zmynd8_exploratory_2026-10-05'
SRC.mkdir(parents=True, exist_ok=True)
NOTE = THESIS / 'research/molecular_characterisation/zmynd8_exploratory_2026-10-05'
NOTE.mkdir(parents=True, exist_ok=True)
ANALYSES = {
    'DU145': PROJECTS / 'RNAseq/analysis/LNCAP/SEPT_23_2026_REANALYZE',
    'HT1080': PROJECTS / 'RNAseq/analysis/HT1080/SEPT_11_2026_REANALYZE',
    'A549': PROJECTS / 'RNAseq/analysis/A549/SEPT_18_2026_REANALYZE',
}
DE = {
    'DU145': ANALYSES['DU145'] / 'outputs/tables/DE_paired_voom_primary.csv',
    'HT1080': FIG / 'source_data/ht1080_batch_aware_de/batchDE_batch_condition__combined.csv',
    'A549': FIG / 'source_data/a549_primary_nuclear_RINe/DE_nuclear_RINe_adjusted.csv',
}
MODELS = {'DU145': 'Passage-blocked limma-voom; 18,767 non-Y/non-PAR genes',
          'HT1080': 'Batch + condition limma-voom; 25,526 non-Y/non-PAR genes',
          'A549': 'RINe-adjusted limma-voom; 22,843 nuclear genes including Y/PAR'}
TARGETS = ['SNAI2', 'CD44', 'VEGFA', 'EGFR', 'MMP1', 'MMP3']
CONTEXT = ['ZMYND8', 'KDM5A', 'KDM5B', 'KDM5C', 'KDM5D', 'CHD4', 'GATAD2A', 'HDAC1', 'HDAC2', 'MMP2', 'MMP7']
GENES = CONTEXT + TARGETS
inputs = set()

def read(path):
    inputs.add(path)
    return pd.read_csv(path, low_memory=False)

rows = []
for line, path in DE.items():
    df = read(path)
    assert df.gene_id.is_unique
    for gene in GENES:
        subset = df[df.gene_name.eq(gene)]
        assert len(subset) <= 1, (line, gene)
        row = {'cell_line': line, 'gene_name': gene, 'model': MODELS[line],
               'source': str(path), 'tested': len(subset) == 1,
               'category': 'Li2016 illustrated corepressed target' if gene in TARGETS else 'context'}
        if len(subset):
            r = subset.iloc[0]
            row.update(gene_id=r.gene_id, logFC=r.logFC, FDR=r.FDR,
                       CI_low=r.get('CI_low', r.get('CI.L', np.nan)),
                       CI_high=r.get('CI_high', r.get('CI.R', np.nan)))
        rows.append(row)
rna = pd.DataFrame(rows)
rna.to_csv(SRC / 'RNA_primary_models.csv', index=False)
z = rna[rna.gene_name.eq('ZMYND8')].set_index('cell_line').loc[list(DE)]
assert z.tested.all()
z.assign(fold_change=2**z.logFC).to_csv(SRC / 'ZMYND8_primary_RNA.csv')

periods = []
for period in ['AUG_2025', 'MAY_2026']:
    d = read(FIG / f'source_data/ht1080_batch_aware_de/batchDE_passage_period__{period}.csv')
    selected = d[d.gene_name.isin(['ZMYND8'] + TARGETS)].copy()
    selected['period'] = period
    periods.append(selected)
pd.concat(periods).to_csv(SRC / 'HT1080_period_specific_RNA.csv', index=False)

du = read(ANALYSES['DU145'] / 'outputs/tables/DU145_LOY_vs_sgCTRL_all_gene_effects.csv')
du[du.gene_name.isin(GENES)].to_csv(SRC / 'DU145_expression_and_passage_effects.csv', index=False)
cpm = read(ANALYSES['DU145'] / 'outputs/tables/core_CPM_all_genes.csv')
du_samples = cpm[cpm.gene_id.eq(z.loc['DU145', 'gene_id'])].copy()
assert len(du_samples) == 6
du_samples['condition'] = np.where(du_samples.sample_key.str.contains('sgCTRL'), 'sgCTRL', 'Y-depleted')
du_samples['passage'] = du_samples.sample_key.str.extract(r'(P\d+)$')
du_samples.to_csv(SRC / 'DU145_ZMYND8_sample_CPM.csv', index=False)
adu = read(ANALYSES['DU145'] / 'outputs/tables/Y_PAR_gene_table.csv')
adu[adu.gene_name.eq('KDM5D')].to_csv(SRC / 'DU145_KDM5D_separate_Y_PAR_model.csv', index=False)
aa = read(ANALYSES['A549'] / 'outputs/tables/A549_LOY_vs_sgCTRL_all_gene_effects.csv')
aa[aa.gene_name.isin(GENES)].to_csv(SRC / 'A549_descriptive_expression_and_filter_status.csv', index=False)

tf = []
for line, analysis in ANALYSES.items():
    d = read(analysis / 'outputs/tf_modules_v1/tables/TF_robustness.csv')
    d = d[d.source.eq('ZMYND8')].copy()
    d.insert(0, 'cell_line', line)
    d['analysis_scope'] = ('Saved unadjusted target-expression analysis; not the RINe-adjusted primary gene model'
                           if line == 'A549' else 'Saved CollecTRI target-expression analysis')
    tf.append(d)
pd.concat(tf).to_csv(SRC / 'ZMYND8_saved_regulator_inference.csv', index=False)

atac = read(FIG / 'source_data/ht1080_rna_atac_peak_scatter/promoter_peak_links.csv')
atac = atac[atac.gene_name.isin(GENES)].copy()
atac.to_csv(SRC / 'HT1080_promoter_peak_results.csv', index=False)
atac.groupby('gene_name').agg(n_tested_peaks=('peak_id', 'nunique'),
                             median_logFC=('logFC', 'median'), min_FDR=('FDR', 'min')).to_csv(
                                 SRC / 'HT1080_promoter_summary.csv')

cn_rows = []
f = ANALYSES['DU145'] / 'outputs/tables/CN_gene_level.csv'
d = read(f); g = d[d.gene_name.eq('ZMYND8')].iloc[0]
cn_rows.append({'cell_line': 'DU145', 'arm': '20q', 'delta_log2_copy_ratio': g['CN_delta_LOY2X'],
                'scope': 'Gene-overlapping segment; LOY-2X minus sgCTRL', 'source': str(f)})
f = FIG / 'source_data/a549_du145_rna_dna/A549_gene_mapping_and_effects.csv'
d = read(f); g = d[d.gene_name.eq('ZMYND8')].iloc[0]
cn_rows.append({'cell_line': 'A549', 'arm': '20q', 'delta_log2_copy_ratio': g.DNA,
                'scope': 'Saved gene-mapped DNA effect; RNA columns here are legacy and not used', 'source': str(f)})
f = FIG / 'source_data/ht1080_atac_copy_number/gene_CN_context.csv'
d = read(f); g = d[d.gene_name.eq('ZMYND8')].iloc[0]
cn_rows.append({'cell_line': 'HT1080', 'arm': '20q', 'delta_log2_copy_ratio': g.CN_promoter_delta,
                'scope': 'Promoter-mapped earlier D0 DNA effect', 'source': str(f)})
pd.DataFrame(cn_rows).to_csv(SRC / 'ZMYND8_local_copy_number_context.csv', index=False)

for f in (THESIS / 'skills/thesis-figure-style/assets/fonts').glob('OpenSans-*.ttf'):
    if 'variable' not in f.name:
        font_manager.fontManager.addfont(str(f))
plt.rcParams.update({'font.family': 'Open Sans', 'font.size': 8, 'axes.labelsize': 8,
                     'axes.titlesize': 9, 'pdf.fonttype': 42, 'svg.fonttype': 'none',
                     'axes.spines.top': False, 'axes.spines.right': False,
                     'axes.linewidth': .65, 'xtick.major.width': .65, 'ytick.major.width': .65})
fig = plt.figure(figsize=(183/25.4, 120/25.4))
gs = fig.add_gridspec(2, 2, width_ratios=[1.04, 1.14], height_ratios=[1, 1.05],
                      left=.11, right=.98, bottom=.16, top=.93, wspace=.52, hspace=.76)
ax = fig.add_subplot(gs[0, 0])
for i, (line, r) in enumerate(z.iterrows()):
    ax.errorbar(r.logFC, i, xerr=[[r.logFC-r.CI_low], [r.CI_high-r.logFC]],
                fmt='o', color='#343434', capsize=2, markersize=4, linewidth=1)
    ax.text(.99, i, f'{r.logFC:+.2f}; {r.FDR:.3f}', ha='right', va='center',
            transform=ax.get_yaxis_transform(), fontsize=7)
ax.axvline(0, color='.7', linestyle='--', linewidth=.7, zorder=0)
ax.set(yticks=range(3), yticklabels=list(z.index), ylim=(2.65, -.7), xlim=(-.55, .99),
       xticks=[-.5, 0, .5], xlabel='ZMYND8 RNA log₂ fold change')
ax.set_title('ZMYND8 expression', loc='left', pad=13)
ax.text(.99, 1.03, 'log₂FC; FDR', transform=ax.transAxes, fontsize=7, ha='right')
ax.text(-.25, 1.11, 'a', transform=ax.transAxes, fontsize=12, fontweight='semibold')

ax = fig.add_subplot(gs[1, 0])
colors = {'sgCTRL': '#3C5488', 'Y-depleted': '#E64B35'}
for passage, sub in du_samples.groupby('passage'):
    sub = sub.set_index('condition')
    ax.plot([0, 1], [sub.loc[c, 'CPM'] for c in colors], color='.8', lw=.8, zorder=1)
for i, (condition, color) in enumerate(colors.items()):
    sub = du_samples[du_samples.condition.eq(condition)]
    ax.scatter(np.full(len(sub), i), sub.CPM, color=color, s=24, zorder=3)
    ax.plot([i-.15, i+.15], [sub.CPM.mean()]*2, color=color, lw=1.6)
ax.set(xlim=(-.45, 1.45), ylim=(0, max(du_samples.CPM.max()*1.14, 100)),
       xticks=[0,1], xticklabels=list(colors), ylabel='ZMYND8 expression (CPM)')
ax.set_title('DU145 passage pairs', loc='left', pad=8)
ax.text(-.25, 1.09, 'b', transform=ax.transAxes, fontsize=12, fontweight='semibold')

ax = fig.add_subplot(gs[:, 1])
eff = rna[rna.gene_name.isin(TARGETS)].pivot(index='gene_name', columns='cell_line', values='logFC').reindex(index=TARGETS, columns=list(DE))
fdr = rna[rna.gene_name.isin(TARGETS)].pivot(index='gene_name', columns='cell_line', values='FDR').reindex_like(eff)
cmap = LinearSegmentedColormap.from_list('blue_white_red', ['#3C5488', '#FFFFFF', '#E64B35'])
cmap.set_bad('#E6E6E6')
im = ax.imshow(eff.values, cmap=cmap, vmin=-5.5, vmax=5.5, aspect='auto')
for i in range(len(TARGETS)):
    for j in range(3):
        v, q = eff.iloc[i,j], fdr.iloc[i,j]
        label = 'Not\ntested' if pd.isna(v) else f'{v:+.2f}' + ('*' if q < .05 else '')
        ax.text(j, i, label, ha='center', va='center', fontsize=8,
                color='white' if pd.notna(v) and abs(v)>3.6 else '#242424')
ax.set(xticks=range(3), xticklabels=list(DE), yticks=range(len(TARGETS)), yticklabels=TARGETS)
ax.tick_params(length=0, pad=5)
ax.set_title('Published ZMYND8–KDM5D targets', loc='left', pad=13)
for spine in ax.spines.values(): spine.set_visible(False)
ax.text(-.31, 1.045, 'c', transform=ax.transAxes, fontsize=12, fontweight='semibold')
cbax = fig.add_axes([.625, .078, .29, .022])
cb = fig.colorbar(im, cax=cbax, orientation='horizontal', ticks=[-5, 0, 5])
cb.set_label('RNA log₂ fold change (Y-depleted / sgCTRL)', fontsize=7, labelpad=3)
cb.ax.tick_params(labelsize=7, length=2)
fig.text(.11, .063, '* FDR < 0.05 in the original full model.\nMMP1: 0.000367; MMP3: 0.0471 (HT1080).', fontsize=7)
stem = FIG / 'ch03_supp_zmynd8_exploratory'
for suffix in ['pdf', 'svg', 'png']:
    fig.savefig(stem.with_suffix('.'+suffix), dpi=300, facecolor='white')
plt.close(fig)

provenance = {'date': '2026-10-05', 'status': 'Exploratory; unselected; not inserted into thesis',
              'comparison_order': list(DE), 'models': MODELS,
              'target_selection': {'genes': TARGETS, 'source': 'Li et al. 2016, Figure 4 and Figure S2',
                                   'doi': '10.1016/j.molcel.2016.06.035',
                                   'scope': 'Six illustrated corepressed loci, not the full published signature'},
              'new_statistical_tests': False, 'intervals': 'Saved limma 95% confidence intervals',
              'inputs': [{'path': str(f), 'sha256': hashlib.sha256(f.read_bytes()).hexdigest()} for f in sorted(inputs)],
              'figure': str(stem.with_suffix('.pdf'))}
(SRC / 'provenance.json').write_text(json.dumps(provenance, indent=2), encoding='utf-8')
(SRC / 'README.md').write_text('Exploratory ZMYND8 review, 5 October 2026.\n\n'
    'Reproduce with `figures/chapter3/code/make_zmynd8_exploratory.py`. No models refitted. '
    'RNA_primary_models.csv preserves the current per-line DE model and original full-family FDR. '
    'Absent genes are untested, not zero. KDM5D is excluded from DU145/HT1080 non-Y primary models; '
    'the separate DU145 Y/PAR model is saved explicitly. '
    'A549 regulator inference predates the RINe-adjusted primary model. '
    'The three local DNA values use their saved mappings (gene/segment versus promoter), not matched RNA specimens.\n\n'
    'The unselected figure is not an active chapter figure. Its exploratory legend and interpretation are in '
    '`research/molecular_characterisation/zmynd8_exploratory_2026-10-05/ZMYND8_exploratory_review.md`.\n', encoding='utf-8')
print('Saved', stem.with_suffix('.pdf'))
print(z[['logFC','FDR','CI_low','CI_high']].to_string())
print(du_samples[['sample_key','CPM']].to_string(index=False))
