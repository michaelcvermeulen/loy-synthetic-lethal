"""Direct WRN-MSI associations within cancer groups; exploratory follow-up."""
from pathlib import Path
import hashlib
import json
import sys
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[3]
FIG = ROOT / 'figures/chapter3'
DATA = FIG / 'source_data/wrn_msi_context_followup'
STEM = 'ch03_supp_wrn_msi_context_followup'
sys.path.insert(0, str(ROOT / 'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
apply_style(8)
points = pd.read_csv(DATA / 'wrn_msi_points.csv')
assert points.ModelID.is_unique and points.Sex.eq('Male').all()
assert np.array_equal(points.MSI.to_numpy(), (points.MSIScore >= 20).astype(int).to_numpy())
levels = {'Lineage': 'OncotreeLineage', 'Primary disease': 'OncotreePrimaryDisease', 'Subtype': 'OncotreeSubtype'}
broader = points[points.benchmark_eligible].copy()
main = broader[broader.y_called].copy()
assert len(main) == 568 and main.MSI.sum() == 26


def design(frame, pooled=False):
    # The MSI term is tested last; retain only independent nuisance columns.
    nuisance = pd.DataFrame({'intercept': np.ones(len(frame))}, index=frame.index)
    candidates = pd.DataFrame(index=frame.index)
    if pooled:
        candidates = pd.get_dummies(frame.OncotreeLineage, prefix='lineage', drop_first=True, dtype=float)
    candidates = pd.concat([candidates, frame[['LibraryKY', 'LibraryHumagneCD']].astype(float)], axis=1)
    for col in candidates:
        trial = np.column_stack([nuisance.to_numpy(), candidates[col].to_numpy()])
        if np.linalg.matrix_rank(trial) > nuisance.shape[1]:
            nuisance[col] = candidates[col]
    exposed = frame.MSI.to_numpy(float)
    residual = exposed - nuisance.to_numpy() @ np.linalg.lstsq(nuisance, exposed, rcond=None)[0]
    if residual @ residual <= 1e-8 * len(frame):
        return None
    nuisance['MSI'] = exposed
    return nuisance


def fit_group(frame, level, group, cohort):
    pos = frame[frame.MSI.eq(1)]
    ref = frame[frame.MSI.eq(0)]
    row = dict(cohort=cohort, level=level, group=group, n=len(frame), n_msi=len(pos), n_reference=len(ref),
               donors_msi=pos.PatientID.nunique(), donors_reference=ref.PatientID.nunique(),
               mean_msi=pos.WRN.mean(), mean_reference=ref.WRN.mean(),
               median_msi=pos.WRN.median(), median_reference=ref.WRN.median(),
               raw_difference=pos.WRN.mean() - ref.WRN.mean(),
               eligible_five_each=min(len(pos), len(ref)) >= 5,
               status='no MSI-positive models' if len(pos) == 0 else 'no reference models' if len(ref) == 0 else 'descriptive only')
    if min(len(pos), len(ref)) < 2:
        return row
    x = design(frame, pooled=level == 'Pan-cancer')
    if x is None:
        row['status'] = 'MSI confounded with covariates'
        return row
    if len(frame) - x.shape[1] < 3:
        row['status'] = 'fewer than three residual degrees of freedom'
        return row
    fit = sm.OLS(frame.WRN, x).fit()
    idx = x.columns.get_loc('MSI')
    robust = fit.get_robustcov_results(cov_type='HC3', use_t=True)
    ci = robust.conf_int(alpha=.05)[idx]
    ols_ci = fit.conf_int(alpha=.05).loc['MSI']
    row.update(status='fitted', effect=fit.params.MSI, se_ols=fit.bse.MSI,
               ci_low_ols=ols_ci.iloc[0], ci_high_ols=ols_ci.iloc[1], p_ols=fit.pvalues.MSI,
               se_hc3=robust.bse[idx], ci_low_hc3=ci[0], ci_high_hc3=ci[1], p_hc3=robust.pvalues[idx],
               residual_df=fit.df_resid, library_terms=';'.join(c for c in x.columns if c.startswith('Library')))
    return row


def all_groups(frame, cohort):
    rows = [fit_group(frame, 'Pan-cancer', 'All lineages', cohort)]
    for level, column in levels.items():
        for group, part in frame.groupby(column, dropna=False, sort=True):
            rows.append(fit_group(part, level, str(group), cohort))
    table = pd.DataFrame(rows)
    table['fdr_hc3_eligible_within_level'] = np.nan
    # One prespecified gene, corrected only across eligible contexts at each level.
    for level in levels:
        idx = table.index[table.level.eq(level) & table.eligible_five_each & table.status.eq('fitted')]
        if len(idx):
            table.loc[idx, 'fdr_hc3_eligible_within_level'] = multipletests(table.loc[idx, 'p_hc3'], method='fdr_bh')[1]
    return table


primary = all_groups(main, 'LOY-screen cohort')
sensitivity = all_groups(broader, 'Broader benchmark cohort')
reference = pd.read_csv(ROOT / 'research/depmap_loy_reanalysis_v1/loy_reanalysis_v1_export_2026-10-07/derived/loy_reanalysis_v1/final/power_benchmarks.csv').query("gene == 'WRN'").iloc[0]
pooled_check = sensitivity.iloc[0]
assert np.isclose(pooled_check.effect, reference.effect, atol=1e-12, rtol=0)
assert np.isclose(pooled_check.se_ols, reference.se, atol=1e-12, rtol=0)
assert np.isclose(pooled_check.p_ols, reference.p, rtol=1e-9, atol=0)
assert pooled_check.n_msi == 27 and pooled_check.n_reference == 617
primary.to_csv(DATA / 'wrn_msi_within_contexts.csv', index=False)
sensitivity.to_csv(DATA / 'wrn_msi_broader_cohort.csv', index=False)
# Check repeated models using a deterministic single model per donor.
dedup = main.sort_values('ModelID').drop_duplicates('PatientID')
donor_check = all_groups(dedup, 'One model per donor')
donor_check.to_csv(DATA / 'wrn_msi_one_model_per_donor.csv', index=False)

shown = primary[primary.status.eq('fitted') & primary.level.ne('Pan-cancer')].copy()
shown = shown.sort_values(['level', 'effect'])
short = {'Esophagus/Stomach': 'Oesophagus/stomach', 'Colorectal Adenocarcinoma': 'Colorectal adenocarcinoma',
         'Esophagogastric Adenocarcinoma': 'Esophagogastric adenocarcinoma',
         'T-Lymphoblastic Leukemia/Lymphoma': 'T-lymphoblastic leukaemia/lymphoma',
         'B-Cell Acute Lymphoblastic Leukemia': 'B-cell ALL', 'T-Cell Acute Lymphoblastic Leukemia': 'T-cell ALL',
         'Renal Clear Cell Carcinoma': 'Clear-cell renal carcinoma'}
fig, axes = plt.subplots(1, 3, figsize=(12.8, 4.15))
fig.subplots_adjust(left=.225, right=.99, bottom=.22, top=.86, wspace=1.25)
for i, (level, ax) in enumerate(zip(levels, axes)):
    tab = shown[shown.level.eq(level)].sort_values('effect')
    for y, row in enumerate(tab.itertuples()):
        color = '#3C5488' if row.eligible_five_each else '#888888'
        ax.errorbar(row.effect, y, xerr=[[row.effect-row.ci_low_hc3], [row.ci_high_hc3-row.effect]],
                    fmt='o', color=color, mfc=color if row.eligible_five_each else 'white',
                    mec=color, ms=4.6, lw=.9, capsize=2)
    ax.axvline(0, color='#BBBBBB', ls='--', lw=.7)
    ax.set_yticks(range(len(tab)), [f'{short.get(row.group, row.group)} ({row.n_msi}/{row.n_reference})' for row in tab.itertuples()], fontsize=7.4)
    ax.set_ylim(len(tab)-.5, -.65)
    ax.set_xlim(shown.ci_low_hc3.min()-.15, shown.ci_high_hc3.max()+.15)
    ax.tick_params(axis='y', length=0)
    ax.set_title(level, fontsize=9, fontweight='semibold', pad=10)
    ax.set_xlabel('WRN Chronos effect\n(MSI-positive − reference)', fontsize=8)
    ax.text(-.1, 1.08, 'abc'[i], transform=ax.transAxes, fontsize=11, weight='semibold')
fig.legend(handles=[Line2D([], [], color='#3C5488', marker='o', lw=.9, label='At least five models per group'),
                    Line2D([], [], color='#888888', marker='o', markerfacecolor='white', lw=.9,
                           label='Fewer than five in a group: exploratory')],
           loc='lower center', bbox_to_anchor=(.54,.025), ncol=2, frameon=False, fontsize=8)
for ext in ['pdf','svg','png']:
    fig.savefig(FIG / f'{STEM}.{ext}', dpi=300, bbox_inches='tight', pad_inches=.045)
plt.close(fig)

# Raw profiles make the small MSI groups and their heterogeneity visible.
contexts = [('Lineage', r.group) for r in primary[primary.level.eq('Lineage') & primary.n_msi.gt(0)]
            .sort_values('raw_difference').itertuples()]
fig, axes = plt.subplots(3, 3, figsize=(9.4,7.1), sharey=True)
axes = axes.ravel()
rng = np.random.default_rng(20261008)
for ax, (level, group) in zip(axes, contexts):
    part = main[main[levels[level]].eq(group)]
    for x, msi, color in [(0,0,'#777777'),(1,1,'#8172B3')]:
        values = part.loc[part.MSI.eq(msi),'WRN'].to_numpy()
        if len(values) >= 2:
            ax.boxplot([values], positions=[x], widths=.5, patch_artist=True, showfliers=False,
                       boxprops=dict(facecolor=color,edgecolor=color,alpha=.2),medianprops=dict(color=color,lw=1.2),
                       whiskerprops=dict(color=color), capprops=dict(color=color))
        ax.scatter(x+rng.uniform(-.15,.15,len(values)),values,s=16,color=color,alpha=.8,edgecolors='white',linewidths=.3)
    ax.axhline(0,color='#BBBBBB',lw=.7,ls='--')
    ax.set_xticks([0,1],[f'Ref.\nn = {sum(part.MSI.eq(0))}',f'MSI+\nn = {sum(part.MSI.eq(1))}'],fontsize=8)
    ax.set_title(short.get(group,group),fontsize=8.5)
    ax.set_xlim(-.5,1.5)
for ax in axes[::3]:
    ax.set_ylabel('WRN gene effect (Chronos)')
for ax in axes[len(contexts):]:
    ax.set_visible(False)
fig.tight_layout(w_pad=1,h_pad=1.5)
for ext in ['pdf','svg','png']:
    fig.savefig(FIG / f'{STEM}_profiles.{ext}',dpi=300,bbox_inches='tight',pad_inches=.045)
plt.close(fig)

report = ['# Direct WRN–MSI associations within cancer groups', '',
          'Exploratory follow-up, 8 October 2026. These are observed WRN–MSI associations, not simulated LOY effects.', '',
          'The primary cohort is the same 568 male models with resolved Y status used in the LOY searches (26 MSI-positive, 542 reference). MSI-positive means the frozen genomic__MSI call (MSIScore ≥20); Y status is not the exposure in these models. The broader 644-model benchmark cohort is retained as sensitivity.', '',
          'Within each lineage, primary disease and subtype: WRN Chronos ~ library indicators + MSI. The pooled model additionally adjusts for lineage. Negative effects indicate stronger WRN dependency in MSI-positive models. Effects are ordinary least-squares coefficients. Displayed intervals and primary exploratory P values use heteroskedasticity-consistent HC3 standard errors and a t reference with residual degrees of freedom; ordinary model-based intervals and P values are also saved for comparison with the original benchmark. No genome-wide scan was run.', '',
          'At least five MSI-positive and five reference models are required for the main comparisons. Effects with two to four models in either group are exploratory and shown with open symbols. Groups with fewer than two in either group, unestimable contrasts or fewer than three residual degrees of freedom retain descriptive summaries only. All groups remain in the tables. FDR is calculated for the single WRN candidate across eligible comparisons separately at each classification level. The nested groups reuse models and are not independent replications.', '',
          'A one-model-per-donor sensitivity uses the alphabetically first ModelID for each PatientID. The pooled broader-cohort effect, standard error and P value reproduce the saved benchmark exactly within numerical tolerance.', '',
          '| Level | Group | MSI / reference | Effect | 95% HC3 interval | Nominal HC3 P |',
          '|---|---|---:|---:|---|---:|']
for row in primary[primary.status.eq('fitted')].itertuples():
    report.append(f'| {row.level} | {row.group} | {row.n_msi}/{row.n_reference} | {row.effect:.3f} | {row.ci_low_hc3:.3f} to {row.ci_high_hc3:.3f} | {row.p_hc3:.3g} |')
report += ['', 'Groups with MSI-positive models that lack an adjusted fit:', '',
           '| Level | Group | MSI / reference | Raw mean difference | Reason |',
           '|---|---|---:|---:|---|']
for row in primary[primary.n_msi.gt(0) & primary.status.ne('fitted')].itertuples():
    value = f'{row.raw_difference:.3f}' if np.isfinite(row.raw_difference) else 'Not estimable'
    report.append(f'| {row.level} | {row.group} | {row.n_msi}/{row.n_reference} | {value} | {row.status} |')
report += ['', '[Forest plot](../../ch03_supp_wrn_msi_context_followup.pdf)',
           '[Raw lineage profiles](../../ch03_supp_wrn_msi_context_followup_profiles.pdf)', '',
           'This follow-up has not been assigned a thesis figure number or inserted into Results.']
(DATA/'README.md').write_text('\n'.join(report)+'\n',encoding='utf-8')
provenance = {'analysis_date':'2026-10-08','source_export':'2026-10-07','new_analysis':True,
              'primary_models':len(main),'primary_msi':int(main.MSI.sum()),'broader_models':len(broader),
              'single_donor_models':len(dedup),'benchmark_reproduced':True,
              'minimum_per_group_main':5,'minimum_per_group_exploratory':2,
              'interval':'HC3 heteroskedasticity-consistent, t with residual df',
              'fdr_family':'One WRN candidate across eligible groups separately within lineage/disease/subtype',
              'input_sha256':hashlib.sha256((DATA/'wrn_msi_points.csv').read_bytes()).hexdigest()}
provenance['result_sha256'] = {p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in DATA.glob('*.csv')}
provenance['figure_sha256'] = {p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                             for p in FIG.glob(f'{STEM}*') if p.suffix in ['.pdf','.svg','.png']}
(DATA/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n',encoding='utf-8')
print(primary[primary.status.eq('fitted')][['level','group','n_msi','n_reference','effect','ci_low_hc3','ci_high_hc3','p_hc3','eligible_five_each','fdr_hc3_eligible_within_level']].to_string(index=False))
print('Broader pooled benchmark reproduced. Tables, forest plot, raw profiles and report saved.')
