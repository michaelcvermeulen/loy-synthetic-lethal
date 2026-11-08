"""Combine saved Y-call, dependency and replication panels without refitting."""
from pathlib import Path
import hashlib
import json
import shutil
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator

ROOT = Path(__file__).resolve().parents[3]
FIG = ROOT / 'figures/chapter3'
SOURCE = ROOT / 'research/depmap_loy_reanalysis_v1/loy_reanalysis_v1_export_2026-10-07/results/loy_reanalysis_v1/figure_source_data'
DATA = FIG / 'source_data/depmap_y_status_dependencies'
sys.path.insert(0, str(ROOT / 'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
font = apply_style(8)
DATA.mkdir(parents=True, exist_ok=True)
names = ['fig1a_rna_score_by_sex', 'fig1b_rna_vs_dna', 'fig1c_call_tiers', 'fig2a_gene_volcano', 'fig2b_hit_boxplots']
tables = {}
hashes = {}
for name in names:
    source = SOURCE / f'{name}.csv'
    target = DATA / source.name
    shutil.copy2(source, target)
    tables[name] = pd.read_csv(target)
    hashes[source.name] = hashlib.sha256(target.read_bytes()).hexdigest()
rna, dna, tiers, genes, box = [tables[name] for name in names]
assert len(genes) == 17205
assert set(genes.loc[genes.fdr < .05, 'symbol']) == {'DDX3X', 'EIF1AX', 'RPS4X', 'ZFX'}
assert tiers.N.sum() == 1143
assert rna.groupby('Sex').size().to_dict() == {'Female': 711, 'Male': 905}
assert box.groupby(['symbol', 'status']).size().unstack().eq([230, 338]).all().all()

# Female models are a descriptive reference, from the same frozen matrix.
female_path = DATA / 'panel_e_female_reference.csv'
female = pd.read_csv(female_path)
assert female.groupby('symbol').size().eq(478).all()
assert set(female.symbol) == {'DDX3X', 'EIF1AX', 'RPS4X', 'ZFX'}
assert female.Sex.eq('Female').all() and female.status.eq('Female').all()
assert female.OncotreePrimaryDisease.ne('Non-Cancerous').all()
assert female.EngineeredModel.fillna('').eq('').all()
assert not female.ModelID.isin(box.ModelID).any()
panel_e = pd.concat([box, female[box.columns]], ignore_index=True)
panel_e.to_csv(DATA / 'panel_e_all_points.csv', index=False)
hashes[female_path.name] = hashlib.sha256(female_path.read_bytes()).hexdigest()


# The frozen analysis cohort is male-derived, cancer and not engineered.
# Prevalence is descriptive across all eligible models, before assay restriction.
cohort_path = SOURCE.parents[2] / 'derived/loy_reanalysis_v1/cohort.csv'
cohort = pd.read_csv(cohort_path, low_memory=False)
assert len(cohort) == 1076 and cohort.ModelID.is_unique
assert cohort.Sex.eq('Male').all()
assert cohort.OncotreePrimaryDisease.ne('Non-Cancerous').all()
assert cohort.EngineeredModel.fillna('').eq('').all()
cohort[['ModelID', 'PatientID', 'OncotreeLineage', 'y_call']].to_csv(DATA / 'lineage_prevalence_models.csv', index=False)
prevalence = cohort.groupby('OncotreeLineage', dropna=False).agg(
    eligible=('ModelID', 'size'),
    y_negative=('y_call', lambda s: s.eq('Y_negative').sum()),
    y_positive=('y_call', lambda s: s.eq('Y_positive').sum())).reset_index()
prevalence['classified'] = prevalence.y_negative + prevalence.y_positive
prevalence['set_aside'] = prevalence.eligible - prevalence.classified
prevalence['loy_percent'] = 100 * prevalence.y_negative / prevalence.classified
prevalence['displayed'] = prevalence.classified.ge(10)
prevalence = prevalence.sort_values(['loy_percent', 'OncotreeLineage'], ascending=[False, True])
prevalence.to_csv(DATA / 'lineage_prevalence_all.csv', index=False)
shown = prevalence[prevalence.displayed].copy()
assert prevalence.classified.sum() == 910 and prevalence.y_negative.sum() == 377
assert len(shown) == 17 and shown.classified.sum() == 875
shown.to_csv(DATA / 'lineage_prevalence_displayed.csv', index=False)

BLUE, RED, GREY, TEAL = '#3C5488', '#E64B35', '#B8B8B8', '#00A087'
colours = {'Y-positive': BLUE, 'Y-negative': RED, 'set aside': GREY}
BASE_HEIGHT, REPLICATION_HEIGHT = 4.9, 3.3
HEIGHT = BASE_HEIGHT + REPLICATION_HEIGHT
fig = plt.figure(figsize=(9.6, HEIGHT))
def axis(rect):
    # Compact a-e above the aligned replication and within-lineage panels.
    x, y, width, height = rect
    ax = fig.add_axes([x, (REPLICATION_HEIGHT + y * BASE_HEIGHT) / HEIGHT,
                       width, height * BASE_HEIGHT / HEIGHT])
    ax.tick_params(length=3, pad=2, labelsize=7.5)
    return ax
def tag(x, y, text):
    fig.text(x, (REPLICATION_HEIGHT + y * BASE_HEIGHT) / HEIGHT, text,
             fontsize=11, weight='semibold', ha='left', va='top')

# a: Preserve the two sexes and their independent frequency scales.
af = axis([.060, .795, .208, .160])
am = axis([.060, .565, .208, .178])
bins = np.arange(0, np.ceil(rna.rna_y_mean.max()/.2)*.2+.201, .2)
for ax, sex, colour in [(af, 'Female', '#7E6148'), (am, 'Male', '#555555')]:
    values = rna.loc[rna.Sex == sex, 'rna_y_mean']
    ax.hist(values, bins=bins, color=colour, edgecolor='white', linewidth=.12)
    for cutoff in [1., 2.5]:
        ax.axvline(cutoff, color='#333333', ls=(0, (3, 2)), lw=.65)
    ax.text(.98, .91, f'{sex}\nn = {len(values):,}', transform=ax.transAxes, ha='right', va='top', fontsize=7.5)
    ax.set_xlim(-.08, bins[-1])
    ax.yaxis.set_major_locator(MaxNLocator(nbins=3, integer=True))
    ax.set_ylabel('Cell lines', labelpad=3)
af.tick_params(labelbottom=False)
am.set_xlabel('Mean log₂(TPM + 1)\nof seven Y-linked genes', labelpad=3)
tag(.010, .98, 'a')

# b: The original points, categories and relative copy-ratio measurement.
b = axis([.367, .565, .253, .390])
for status in ['set aside', 'Y-positive', 'Y-negative']:
    part = dna[dna.status == status]
    b.scatter(part.dna_y_ratio, part.rna_y_mean, c=colours[status], s=8, alpha=.70, linewidths=0)
b.set_xlabel('Y copy ratio\n(median of seven genes)', labelpad=3)
b.set_ylabel('Mean log₂(TPM + 1)\nof seven Y-linked genes', labelpad=3)
b.xaxis.set_major_locator(MaxNLocator(nbins=4))
b.yaxis.set_major_locator(MaxNLocator(nbins=5))
b.legend(handles=[Line2D([], [], marker='o', color='none', markerfacecolor=colours[s], markeredgewidth=0, markersize=4, label=s) for s in colours], loc='lower right', fontsize=7, handletextpad=.4, labelspacing=.25, borderaxespad=.4)
tag(.300, .98, 'b')

# c: Y-negative percentage among classified cancer models, with denominators.
c = axis([.817, .565, .126, .390])
display_names = {
    'Esophagus/Stomach': 'Oesophagus/stomach', 'Head and Neck': 'Head and neck',
    'Bladder/Urinary Tract': 'Bladder/urinary tract', 'Biliary Tract': 'Biliary tract',
    'Soft Tissue': 'Soft tissue', 'Peripheral Nervous System': 'Peripheral nervous system',
    'CNS/Brain': 'CNS/brain'}
ys = np.arange(len(shown))
c.barh(ys, shown.loy_percent, color=RED, height=.68)
c.set_yticks(ys, [display_names.get(v, v) for v in shown.OncotreeLineage], fontsize=6.5)
c.set_ylim(len(shown)-.4, -.6)
c.set_xlim(0, 115)  # Leave room for labels beyond the longest bar.
c.set_xticks([0, 50, 100])
c.set_xlabel('Y-negative (%)', labelpad=3)
c.tick_params(axis='y', length=0, pad=4)
c.spines['left'].set_visible(False)
for y, row in zip(ys, shown.itertuples()):
    c.annotate(f'{row.loy_percent:.1f}%', (row.loy_percent, y),
               xytext=(2, 0), textcoords='offset points',
               ha='left', va='center', fontsize=6.5, clip_on=False)
    c.text(1.045, y, f'{row.y_negative}/{row.classified}', transform=c.get_yaxis_transform(),
           ha='left', va='center', fontsize=6.5, clip_on=False)
c.text(1.045, 1.02, 'Y− / n', transform=c.transAxes, fontsize=6.5, ha='left', va='bottom')
tag(.655, .98, 'c')

# The former c panel is retained in a standalone supplementary figure.
tier_order = ['RNA and DNA agree', 'DNA only', 'RNA only', 'RNA, DNA partial', 'RNA intermediate', 'DNA partial, no RNA', 'no data', 'RNA and DNA disagree']
labels = ['RNA and DNA agree', 'DNA only', 'RNA only', 'RNA; DNA partial', 'RNA intermediate', 'DNA partial; no RNA', 'No data', 'RNA and DNA disagree']

# d: Preserve raw-P vertical axis, BH-FDR categories and primary effect sizes.
d = axis([.060, .090, .280, .342])
classes = [('other', GREY, 3), ('X-Y partner, not a hit', TEAL, 13), ('hit', RED, 20)]
for kind, colour, size in classes:
    part = genes[genes['class'] == kind]
    d.scatter(part.estimate, -np.log10(part.p), s=size, color=colour, alpha=.8, linewidths=0)
d.axvline(0, color='#BBBBBB', lw=.5, zorder=0)
d.set_xlabel('Chronos effect: Y− minus Y+', labelpad=3)
d.set_ylabel('−log₁₀ P', labelpad=3)
d.set_xlim(genes.estimate.min()-.12, genes.estimate.max()+.12)
d.set_ylim(-1.5, 74)
d.yaxis.set_major_locator(MaxNLocator(nbins=5))
offsets = {'DDX3X': (9, 3), 'EIF1AX': (7, -8), 'ZFX': (12, 7), 'RPS4X': (-39, 5), 'KDM5C': (-45, 27), 'USP9X': (-44, 5), 'PRKX': (9, 5)}
annotate = genes[(genes.fdr < .05) | ((genes['class'] == 'X-Y partner, not a hit') & (genes.p < .01))]
for row in annotate.itertuples():
    d.annotate(row.symbol, (row.estimate, -np.log10(row.p)), xytext=offsets.get(row.symbol, (9, 8)), textcoords='offset points', fontsize=7.5, fontstyle='italic', arrowprops=dict(arrowstyle='-', lw=.45, color='#666666', shrinkA=1, shrinkB=3))
d.legend(handles=[Line2D([], [], marker='o', color='none', markerfacecolor=colour, markeredgewidth=0, markersize=4, label=label) for label, colour in [('FDR <0.05', RED), ('Other X–Y paralog', TEAL), ('Other gene', GREY)]], loc='upper right', fontsize=6.8, handletextpad=.35, labelspacing=.2, borderaxespad=.2)
tag(.010, .465, 'd')

# e: Male comparisons plus female reference values; saved male statistics are unchanged.
rng = np.random.default_rng(20261005)
female_rng = np.random.default_rng(20261008)
e_colours = {**colours, 'Female': '#7E6148'}
gene_order = ['DDX3X', 'EIF1AX', 'ZFX', 'RPS4X']
for i, gene in enumerate(gene_order):
    e = axis([.437+i*.140, .090, .117, .342])
    frame = panel_e[panel_e.symbol == gene]
    values = []
    for x, status in enumerate(['Y-positive', 'Y-negative', 'Female'], start=1):
        v = frame.loc[frame.status == status, 'chronos'].to_numpy()
        values.append(v)
        jitter_rng = female_rng if status == 'Female' else rng
        e.scatter(x+jitter_rng.uniform(-.19, .19, len(v)), v, s=3.5, alpha=.40, color=e_colours[status], linewidths=0)
    e.boxplot(values, positions=[1, 2, 3], widths=.45, showfliers=False, whis=1.5,
              boxprops=dict(lw=.65), medianprops=dict(lw=.9, color='#222222'),
              whiskerprops=dict(lw=.65), capprops=dict(lw=.65))
    row = genes[genes.symbol == gene].iloc[0]
    lo, hi = frame.chronos.min(), frame.chronos.max()
    e.set_ylim(lo-.07*(hi-lo), hi+.39*(hi-lo))
    e.set_xlim(.48, 3.60)
    # Bracket only the two male groups to identify the saved comparison.
    e.plot([1, 1, 2, 2], [.78, .80, .80, .78], transform=e.get_xaxis_transform(),
           color='#222222', lw=.75, solid_capstyle='butt', zorder=5)
    e.set_xticks([1, 2, 3], ['Male\nY+', 'Male\nY−', 'Female'])
    e.tick_params(axis='x', labelsize=6.7)
    e.yaxis.set_major_locator(MaxNLocator(nbins=4))
    e.text(.5, 1.025, gene, fontstyle='italic', weight='semibold', transform=e.transAxes, ha='center', va='bottom', fontsize=8.5)
    exp = int(np.floor(np.log10(row.fdr)))
    mantissa = row.fdr / 10**exp
    e.text(.5, .975, f'Δ = {row.estimate:.2f}\nFDR {mantissa:.1f} × 10$^{{{exp}}}$', transform=e.transAxes, ha='center', va='top', fontsize=6.8, linespacing=1.3)
    if i == 0:
        e.set_ylabel('Chronos gene effect', labelpad=3)
tag(.375, .465, 'e')

# f: The former F02, compacted into a two-by-two gene facet layout.
replication_path = FIG / 'source_data/depmap_core_replication/plotted_estimates.csv'
replication = pd.read_csv(replication_path)
assert len(replication) == 11
replication_genes = ['DDX3X', 'EIF1AX', 'ZFX', 'RPS4X']
replication_sources = ['Cell lines, CRISPR', 'Sanger organoids, CRISPR', 'NextGen, CRISPR']
replication_colours = ['#444444', '#009E87', '#8570B2']
replication_labels = ['Cell lines', 'Organoids', 'NextGen']
fig.text(.010, 3.26 / HEIGHT, 'f', fontsize=11, weight='semibold', ha='left', va='top')
for index, gene in enumerate(replication_genes):
    facet_row, facet_col = divmod(index, 2)
    f = fig.add_axes([.103 + facet_col * .164, (1.99 - facet_row * 1.31) / HEIGHT, .143, .95 / HEIGHT])
    f.axvline(0, color='#929292', lw=.65, ls=(0, (3, 3)), zorder=0)
    f.text(.31, 1.12, gene, transform=f.transAxes, fontstyle='italic',
           weight='semibold', fontsize=8.5, ha='center')
    f.text(.91, 1.12, 'FDR', transform=f.transAxes, ha='center', fontsize=6.8)
    for y, source, colour in zip([2, 1, 0], replication_sources, replication_colours):
        matches = replication[(replication.symbol == gene) & (replication.source == source)]
        if matches.empty:
            assert gene == 'EIF1AX' and source == 'Sanger organoids, CRISPR'
            f.text(-1.15, y, 'Not in library', va='center', ha='center', fontsize=6.6, color='#777777')
            f.text(.78, y, '—', va='center', ha='center', fontsize=6.8)
            continue
        assert len(matches) == 1
        row = matches.iloc[0]
        effect, low, high, fdr = [float(row[k]) for k in ('smd', 'smd_low', 'smd_high', 'display_fdr')]
        assert low <= effect <= high
        f.errorbar(effect, y, xerr=[[effect-low], [high-effect]], fmt='o', markersize=4.5,
                   color=colour, markerfacecolor=colour if fdr < .05 else 'white',
                   markeredgewidth=.9, elinewidth=.85, capsize=2, zorder=3)
        qtext = f'{fdr:.1e}' if fdr < .001 else f'{fdr:.3f}'
        f.text(.78, y, qtext, va='center', ha='center', fontsize=6.8)
    f.set_yticks([2, 1, 0], replication_labels if facet_col == 0 else ['', '', ''])
    f.tick_params(axis='y', length=0, pad=4, labelsize=7)
    f.tick_params(axis='x', length=3, pad=2, labelsize=7)
    f.set_ylim(-.55, 2.65)
    f.set_xlim(-2.3, 1.08)
    f.set_xticks([-2, -1, 0])
    f.spines['left'].set_visible(False)
    f.spines['bottom'].set_bounds(-2.3, .18)
fig.text(.215, .28 / HEIGHT, 'Standardised dependency difference\n(Y-negative − Y-positive)',
         ha='center', va='bottom', fontsize=7.5, linespacing=1.25)

# g: Saved forest plots from supplementary S01b; no model or interval refit.
lineage_path = SOURCE / 'fig3b_hits_within_lineage.csv'
shutil.copy2(lineage_path, DATA / lineage_path.name)
lineage = pd.read_csv(lineage_path)
assert len(lineage) == 68 and lineage.variant.nunique() == 17
assert set(lineage.symbol) == set(gene_order)
assert lineage.groupby('variant').size().eq(4).all()
assert lineage.n.ge(10).all() and lineage.n_negative.ge(3).all() and lineage.n_positive.ge(3).all()
assert (lineage.estimate < 0).sum() == 66
lineage['ci_low'] = lineage.estimate - 1.96 * lineage.se
lineage['ci_high'] = lineage.estimate + 1.96 * lineage.se
lineage.to_csv(DATA / 'lineage_dependency_plotted.csv', index=False)
lineage_order = lineage.groupby('variant', as_index=False).n.mean().sort_values(['n','variant'], ascending=[False,False]).variant.tolist()
fig.text(.430, 3.26 / HEIGHT, 'g', fontsize=11, weight='semibold', ha='left', va='top')
for index, gene in enumerate(gene_order):
    g = fig.add_axes([.565 + index * .11, .48 / HEIGHT, .093, 2.54 / HEIGHT])
    q = lineage[lineage.symbol == gene].set_index('variant').loc[lineage_order]
    positions = np.arange(len(q))
    g.axvline(0, color='#929292', lw=.6, ls=(0, (2, 2)), zorder=0)
    g.errorbar(q.estimate, positions, xerr=[q.estimate-q.ci_low, q.ci_high-q.estimate],
               fmt='o', color=RED, markersize=2.7, markeredgewidth=0,
               elinewidth=.65, capsize=0, zorder=3)
    lo, hi = min(q.ci_low.min(), 0), max(q.ci_high.max(), 0)
    pad = (hi-lo)*.10
    g.set_xlim(lo-pad, hi+pad)
    g.set_ylim(len(q)-.45, -.6)
    g.set_yticks(positions, [display_names.get(v,v) for v in lineage_order] if index == 0 else ['']*len(q))
    g.tick_params(axis='y', length=0, pad=4, labelsize=6.4)
    g.tick_params(axis='x', length=3, pad=2, labelsize=6.8)
    g.xaxis.set_major_locator(MaxNLocator(nbins=3))
    g.spines['left'].set_visible(index == 0)
    g.text(.5, 1.035, gene, transform=g.transAxes, ha='center', va='bottom',
           fontstyle='italic', weight='semibold', fontsize=8.5)
fig.text(.777, .085 / HEIGHT, 'Within-lineage Chronos effect (Y-negative − Y-positive; 95% CI)',
         ha='center', va='bottom', fontsize=7.5)

stem = 'ch03_main_depmap_y_status_dependencies'
for ext in ['pdf', 'svg', 'png']:
    fig.savefig(FIG / f'{stem}.{ext}', dpi=300, facecolor='white')
plt.close(fig)

# Call-evidence exports are independent; rebuild only when explicitly requested.
if '--with-call-evidence' in sys.argv:
    supp, s = plt.subplots(figsize=(6.0, 3.1))
    supp.subplots_adjust(left=.32, right=.96, bottom=.16, top=.90)
    positions = np.arange(len(tier_order))
    left = np.zeros(len(tier_order))
    for status in colours:
        values = tiers[tiers.status == status].groupby('tier').N.sum().reindex(tier_order, fill_value=0).to_numpy()
        s.barh(positions, values, left=left, height=.67, color=colours[status], edgecolor='white', linewidth=.35, label=status)
        for y, start, value in zip(positions, left, values):
            if value >= 40:
                s.text(start+value/2, y, str(value), va='center', ha='center',
                       color='white' if status != 'set aside' else '#333333', fontsize=8)
        left += values
    s.set_yticks(positions, labels, fontsize=8)
    s.invert_yaxis()
    s.set_xlabel('Male cell lines')
    s.set_xlim(0, 795)
    s.set_xticks([0, 200, 400, 600, 800])
    s.tick_params(axis='y', length=0, pad=5)
    s.spines['left'].set_visible(False)
    s.legend(loc='lower left', bbox_to_anchor=(-.04, 1.015), ncol=3, fontsize=8,
             handlelength=1.2, columnspacing=1.2, borderaxespad=0)
    supp_stem = 'ch03_supp_depmap_y_call_evidence'
    for ext in ['pdf', 'svg', 'png']:
        supp.savefig(FIG / f'{supp_stem}.{ext}', dpi=300, facecolor='white')
    plt.close(supp)
(DATA / 'provenance.json').write_text(json.dumps({
    'source': str(SOURCE), 'source_sha256': hashes, 'font': font,
    'panel_mapping': {'a': names[0], 'b': names[1], 'c': 'lineage_prevalence_displayed', 'd': names[3], 'e': 'panel_e_all_points.csv',
                      'f': 'depmap_core_replication/plotted_estimates.csv',
                      'g': 'fig3b_hits_within_lineage.csv'},
    'replication': {'source': str(replication_path.relative_to(ROOT)),
                    'source_sha256': hashlib.sha256(replication_path.read_bytes()).hexdigest(),
                    'source_manifest': 'figures/chapter3/source_data/depmap_core_replication/source_manifest.json',
                    'estimates': len(replication), 'missing': 'EIF1AX in Sanger library',
                    'gene_order': replication_genes, 'source_order': replication_sources,
                    'former_figure': 'CH3-S3.3.1-F02', 'layout': '2 x 2 facets'},
    'within_lineage': {'source': str(lineage_path),
                      'source_sha256': hashlib.sha256(lineage_path.read_bytes()).hexdigest(),
                      'former_panel': 'CH3-S3.3.2-S01b', 'estimates': 68, 'lineages': 17,
                      'minimum_lineage_models': 10, 'minimum_models_per_y_group': 3,
                      'interval': 'estimate +/- 1.96 * saved donor-clustered robust SE',
                      'lineage_order': lineage_order, 'gene_order': gene_order},
    'supplementary_panel': {'id': 'CH3-S3.4.1-S01', 'source': names[2]},
    'lineage_prevalence': {'source': str(cohort_path), 'source_sha256': hashlib.sha256(cohort_path.read_bytes()).hexdigest(),
                           'eligible_models': len(cohort), 'classified': 910, 'y_negative': 377, 'y_positive': 533,
                           'set_aside_excluded': 166, 'minimum_classified_per_displayed_lineage': 10,
                           'displayed_lineages': 17, 'displayed_models': 875, 'grouping': 'OncotreeLineage',
                           'definition': '100 * Y-negative / (Y-negative + Y-positive)'},
    'female_reference': {'source_table': 'panel_e_female_reference.csv',
                         'extraction_script': 'figures/chapter3/code/ch03_main_depmap_female_reference_extract.R',
                         'source_manifest': 'figures/chapter3/source_data/depmap_y_status_dependencies/female_reference_provenance.json',
                         'models_per_gene': 478, 'donors': 457,
                         'eligibility': 'Annotated female cancer models; not engineered; finite Chronos values',
                         'interpretation': 'Descriptive reference; effects and FDRs compare male Y-negative with male Y-positive'},
    'panel_e_comparison_brackets': 'Male Y-positive versus male Y-negative; no bracket to female reference',
    'figure_inches': [9.6, HEIGHT], 'statistics_refitted': False,
    'random_seed_for_visual_jitter': 20261005,
    'former_figures': ['CH3-S3.4.1-F01', 'CH3-S3.4.2-F01']
}, indent=2), encoding='utf-8')
print(json.dumps({'figure': str(FIG / f'{stem}.png'), 'font': font, 'panels': 7,
                  'replication_estimates': len(replication), 'source_counts': {k: len(v) for k, v in tables.items()}}))
