"""Display saved drug controls, pooled compounds and classes; do not refit tests."""
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
from matplotlib.colors import Normalize

ROOT = Path(__file__).resolve().parents[3]
FIG = ROOT / 'figures/chapter3'
EXPORT = ROOT / 'research/depmap_loy_reanalysis_v1/loy_reanalysis_v1_export_2026-10-07'
SOURCE = EXPORT / 'derived/loy_reanalysis_v1/final'
DATA = FIG / 'source_data/depmap_drug_controls_pooling'
STEM = 'ch03_supp_depmap_drug_controls_pooling'
sys.path.insert(0, str(ROOT / 'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
font = apply_style(8)
DATA.mkdir(parents=True, exist_ok=True)
input_names = ['pipeline_controls.csv', 'drug_feature_annotation.csv', 'drug_pan_consensus.csv',
               'drug_pan_mechanism.csv', 'forest_y_status_summary.csv', 'drug_pan_results.csv']
tables, manifest = {}, {}
for name in input_names:
    path = SOURCE / name
    tables[name] = pd.read_csv(path)
    shutil.copy2(path, DATA / name)
    manifest[name] = {'path': str(path.relative_to(ROOT)), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}

controls = tables['pipeline_controls.csv'].query("dataset != 'crispr_effect'").copy()
controls = controls.merge(tables['drug_feature_annotation.csv'][['dataset', 'feature', 'compound_name']],
                          on=['dataset', 'feature'], how='left', validate='one_to_one')
assert len(controls) == 12 and controls.compound_name.notna().all()
assert controls['pass'].all() and controls.p.lt(.01).all()
assert np.array_equal(np.sign(controls.estimate), controls.expected_sign)
pooled = tables['drug_pan_consensus.csv']
classes = tables['drug_pan_mechanism.csv']
rf = tables['forest_y_status_summary.csv']
assert len(pooled) == 352 and pooled.fdr.ge(.05).all()
assert pooled.families.ge(2).all()
assert classes.groupby('class_type').size().to_dict() == {'mechanism': 218, 'target': 385}
assert classes.fdr.ge(.05).all()
assert rf.query("dataset != 'crispr_effect'").stable.eq(0).all()
controls.to_csv(DATA / 'plotted_controls.csv', index=False)

order = ['prism_REPPRIMARY', 'prism_REP1M', 'prism_REP300', 'prism_secondary', 'gdsc1', 'gdsc2', 'ctrp']
labels = ['PRISM REP.PRIMARY', 'PRISM REP.1M', 'PRISM REP.300', 'PRISM secondary', 'GDSC1', 'GDSC2', 'CTRP v2']
types = ['BRAF altered -> BRAF inhibitor sensitivity', 'TP53 altered -> MDM2 inhibitor resistance']
matrix = np.full((7, 2), np.nan)
for row in controls.itertuples():
    matrix[order.index(row.dataset), types.index(row.control)] = -np.log10(row.p)

fig = plt.figure(figsize=(8.9, 6.35))
a = fig.add_axes([.165, .525, .35, .375])
colourmap = plt.get_cmap('magma_r').copy()
colourmap.set_bad('#eeeeee')
im = a.imshow(matrix, cmap=colourmap, norm=Normalize(0, 32), aspect='auto')
a.set_yticks(np.arange(7), labels, fontsize=7.5)
a.set_xticks([0, 1], ['BRAF alteration\nBRAF inhibitors', 'TP53 alteration\nMDM2 inhibitors'], fontsize=8)
a.xaxis.tick_top()
a.tick_params(length=0, pad=5)
for spine in a.spines.values():
    spine.set_visible(False)
a.set_xticks(np.arange(-.5, 2, 1), minor=True)
a.set_yticks(np.arange(-.5, 7, 1), minor=True)
a.grid(which='minor', color='white', linewidth=1.2)
a.tick_params(which='minor', length=0)
for i, dataset in enumerate(order):
    for j, control in enumerate(types):
        rows = controls[(controls.dataset == dataset) & (controls.control == control)]
        if rows.empty:
            a.text(j, i, 'Not tested', ha='center', va='center', fontsize=7, color='#777777')
            continue
        row = rows.iloc[0]
        colour = 'white' if matrix[i, j] > 17 else '#222222'
        name = row.compound_name.title().replace('Plx', 'PLX').replace('Nutlin-3A', 'Nutlin-3a')
        a.text(j, i, f'{name}\nP = {row.p:.1e}\nn = {row.n_exposed}/{row.n_unexposed}',
               ha='center', va='center', fontsize=6.8, linespacing=1.13, color=colour)
a.text(.25, -.058, 'Altered: more sensitive', transform=a.transAxes, ha='center', va='top', fontsize=7)
a.text(.75, -.058, 'Altered: less sensitive', transform=a.transAxes, ha='center', va='top', fontsize=7)
cax = fig.add_axes([.535, .575, .010, .245])
bar = fig.colorbar(im, cax=cax, ticks=[0, 10, 20, 30])
bar.set_label('−log₁₀ nominal P', fontsize=7, labelpad=4)
bar.ax.tick_params(labelsize=6.5, length=2)
bar.outline.set_visible(False)

b = fig.add_axes([.69, .55, .285, .33])
c = fig.add_axes([.105, .10, .365, .29])
d = fig.add_axes([.61, .10, .365, .29])

def association_panel(ax, frame, x, xlim, ticks, xlabel, family_label, colour):
    ax.scatter(frame[x], -np.log10(frame.p), s=13, color=colour, alpha=.65, linewidths=0)
    ax.axvline(0, color='#c9c9c9', lw=.6, zorder=0)
    nominal = -np.log10(.05)
    ax.axhline(nominal, color='#858585', lw=.7, ls=(0, (3, 3)), zorder=0)
    ax.text(.97, nominal + .08, 'Nominal P = 0.05', transform=ax.get_yaxis_transform(),
            fontsize=6.5, ha='right', color='#666666')
    ax.set(xlim=xlim, ylim=(-.05, 3.35), xticks=ticks, yticks=[0, 1, 2, 3], xlabel=xlabel, ylabel='−log₁₀ nominal P')
    ax.tick_params(length=3, pad=2, labelsize=7.5)
    ax.text(.5, 1.10, family_label, transform=ax.transAxes, ha='center', fontsize=8.5)
    ax.text(.5, 1.025, f'0/{len(frame):,} at FDR <0.05; minimum FDR = {frame.fdr.min():.3f}',
            transform=ax.transAxes, ha='center', fontsize=7)

association_panel(b, pooled, 'estimate', (-.43, .43), [-.4, -.2, 0, .2, .4],
                  'Pooled LOY-associated response difference (SD)', 'Shared compounds', '#555555')
association_panel(c, classes[classes.class_type.eq('mechanism')], 'mean_t', (-1.1, 1.75), [-1, 0, 1],
                  'Mean signed drug-response t statistic', 'Mechanism-of-action classes', '#517491')
association_panel(d, classes[classes.class_type.eq('target')], 'mean_t', (-1.1, 1.75), [-1, 0, 1],
                  'Mean signed drug-response t statistic', 'Target classes', '#517491')

# Label the four lowest nominal P values per panel, retaining cutoff ties.
# Coincident mechanism classes share a callout so their labels do not overlap.
label_layout = {
    'b': {
        'TOSEDOSTAT': ('Tosedostat', -.40, 3.16, 'left'),
        'ABT-737': ('ABT-737', .27, 2.86, 'center'),
        'GSK429286A': ('GSK429286A', .035, 2.04, 'center'),
        'LENALIDOMIDE': ('Lenalidomide', -.405, 2.26, 'left'),
        'BERZOSERTIB': ('Berzosertib', -.39, 1.08, 'left'),
        'ETP-46464': ('ETP-46464', -.40, .79, 'left'),
        'BEZ235': ('BEZ235', -.37, .32, 'left'),
        'VE 821': ('VE-821', .29, .54, 'center'),
        'CERALASERTIB': ('Ceralasertib', .29, .17, 'center'),
    },
    'c': {
        'TACHYKININ ANTAGONIST': ('Tachykinin\nantagonist', .60, 2.86, 'center'),
        'SODIUM CHANNEL BLOCKER': ('Sodium channel\nblocker', -.18, 2.19, 'center'),
        'FGFR2': ('FGFR2; inhibitor\nof FGFR1', 1.30, 2.41, 'center'),
        'INHIBITOR OF FGFR1': None,
    },
    'd': {
        'SLC12A1': ('SLC12A1', .98, 3.10, 'center'),
        'SLC12A2': ('SLC12A2', 1.49, 3.10, 'center'),
        'TRPV1': ('TRPV1', .37, 2.78, 'center'),
        'GPR35': ('GPR35', .54, 2.05, 'center'),
        'ATR': ('ATR', -.91, 1.13, 'center'),
    },
}
label_records = []
for panel, ax, frame, xcolumn, namecolumn in [
    ('b', b, pooled, 'estimate', 'compound_name'),
    ('c', c, classes[classes.class_type.eq('mechanism')], 'mean_t', 'class'),
    ('d', d, classes[classes.class_type.eq('target')], 'mean_t', 'class'),
]:
    ordered = frame.sort_values(['p', namecolumn])
    selected_labels = ordered[ordered.p.le(ordered.iloc[3].p)].copy()
    selected_labels['selection_reason'] = 'top nominal P'
    if panel == 'b':
        atr = ordered[ordered.targets.fillna('').str.contains(r'(?:^|[;,])\s*ATR\s*(?:$|[;,])')].copy()
        assert set(atr.compound_name) == {'BERZOSERTIB', 'ETP-46464', 'BEZ235', 'VE 821', 'CERALASERTIB'}
        assert not selected_labels[namecolumn].isin(atr[namecolumn]).any()
        atr['selection_reason'] = 'ATR target annotation; author-requested label'
        selected_labels = pd.concat([selected_labels, atr], ignore_index=True)
    if panel == 'd':
        atr = ordered[ordered[namecolumn].eq('ATR')].copy()
        assert len(atr) == 1 and not selected_labels[namecolumn].eq('ATR').any()
        atr['selection_reason'] = 'literature-motivated candidate; author-requested label'
        selected_labels = pd.concat([selected_labels, atr], ignore_index=True)
    assert set(selected_labels[namecolumn]) == set(label_layout[panel])
    if panel == 'c':
        overlapping = selected_labels[selected_labels['class'].isin(['FGFR2', 'INHIBITOR OF FGFR1'])]
        assert overlapping[['mean_t', 'p']].drop_duplicates().shape[0] == 1
    for _, row in selected_labels.iterrows():
        source_label = row[namecolumn]
        label_records.append({'panel': panel, 'source_label': source_label,
                              'x': row[xcolumn], 'p': row.p, 'fdr': row.fdr,
                              'selection_reason': row['selection_reason']})
        layout = label_layout[panel][source_label]
        if layout is None:
            continue
        text, tx, ty, ha = layout
        point = (row[xcolumn], -np.log10(row.p))
        ax.scatter(*point, s=20, color='#333333', linewidths=0, zorder=4)
        ax.annotate(text, point, xytext=(tx, ty), ha=ha, va='center', fontsize=7,
                    linespacing=1.12, color='#292929',
                    bbox=dict(facecolor='white', edgecolor='none', pad=.6),
                    arrowprops=dict(arrowstyle='-', lw=.6, color='#666666', shrinkA=2, shrinkB=3),
                    zorder=5)
pd.DataFrame(label_records).to_csv(DATA / 'labelled_associations.csv', index=False)

# ATR is absent from the tested mechanism-class table: preserve that absence.
# Reproduce only its eligibility counts, without combining classes or fitting tests.
drug_results = tables['drug_pan_results.csv']
atr_mechanism_rows = []
for _, row in drug_results.iterrows():
    for term in str(row.mechanism).split(','):
        term = term.strip()
        if 'ATR' in term.split():
            atr_mechanism_rows.append({'class': term,
                'compound': row.compound_id if pd.notna(row.compound_id) else row.compound_name})
atr_mechanism_counts = pd.DataFrame(atr_mechanism_rows).groupby('class').agg(
    compounds=('compound', 'nunique'), tests=('compound', 'size')).reset_index()
assert atr_mechanism_counts.set_index('class').compounds.to_dict() == {
    'ATR KINASE INHIBITOR': 3, 'INHIBITOR OF ATR': 2, 'INHIBITOR OF ATR AND MTOR': 1}
assert not classes[classes.class_type.eq('mechanism')]['class'].isin(atr_mechanism_counts['class']).any()
atr_mechanism_counts.to_csv(DATA / 'atr_mechanism_eligibility.csv', index=False)
c.text(-.98, 3.17, 'ATR mechanism groups:\n1–3 compounds each;\nnot tested (minimum 5)',
       ha='left', va='top', fontsize=7, linespacing=1.2, color='#555555')

for x, y, letter in [(.015, .958, 'a'), (.615, .958, 'b'), (.025, .445, 'c'), (.53, .445, 'd')]:
    fig.text(x, y, letter, fontsize=11, weight='semibold', va='top')
fig.text(.54, .017, 'Negative values in b–d indicate greater sensitivity in Y-negative models.', ha='center', fontsize=7.5)

for ext in ['pdf', 'svg', 'png']:
    fig.savefig(FIG / f'{STEM}.{ext}', dpi=300, facecolor='white')
plt.close(fig)
provenance = {
    'inputs': manifest, 'font': font, 'statistics_refitted': False, 'figure_inches': [8.9, 6.35],
    'control_selection': 'Saved minimum-P eligible matching inhibitor per genotype and dataset; nominal P, not FDR.',
    'control_n_order': 'genotype-altered / unaltered; no control in REP.300',
    'pooled_family': '352 compounds present in at least two screening families; PRISM variants form one family.',
    'pooled_test': 'Saved mixed-model likelihood-ratio P; standardised response, lineage and source fixed effects, ModelID random intercept.',
    'class_test': 'Saved two-sided within-lineage permutation P, 2,000 permutations; BH separately for mechanisms and targets; >=5 compounds/class.',
    'threshold_line': 'Nominal P=0.05 only. FDR results are reported separately above each panel.',
    'labels': 'Four lowest nominal P values per panel b-d, including ties at the fourth value, plus five ATR-annotated pooled compounds in b and ATR in d at the author request; coincident FGFR2 and INHIBITOR OF FGFR1 classes share a label. Highlighting denotes labelled points, not FDR significance.',
    'atr_mechanism_note': 'No ATR mechanism class was tested: the three distinct source annotations have 3, 2 and 1 compounds, below the existing minimum of 5. Panel c displays an explanatory note, not a new data point or pooled class test.',
    'counts': {'drug_controls': 12, 'pooled_compounds': 352, 'mechanism_classes': 218, 'target_classes': 385},
    'forest_reference': 'Existing CH3-S3.3.3-S04 panel c; drug stable counts all zero in current export.'
}
(DATA / 'provenance.json').write_text(json.dumps(provenance, indent=2), encoding='utf-8')
print(json.dumps({'figure': str(FIG / f'{STEM}.png'), **provenance['counts']}))
