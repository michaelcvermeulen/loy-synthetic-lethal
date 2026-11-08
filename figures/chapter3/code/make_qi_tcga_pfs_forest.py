"""Restyle published Qi survival estimates; no patient-level models are refitted."""
from pathlib import Path
import sys
import json
import hashlib
import numpy as np
import pandas as pd
from scipy.stats import norm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FixedFormatter, NullLocator

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / 'figures/chapter3'
# Directory retains its initial task name; the confirmed source endpoint is PFS.
DATA = OUT / 'source_data/qi_tcga_pfi'
sys.path.insert(0, str(ROOT / 'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style, condition_colours
font = apply_style(9)
red = condition_colours()['Y-depleted']
src = DATA / 'Qi_Table_S4.xlsx'
d = pd.read_excel(src, sheet_name='Table S4')
assert len(d) == 26 and d.Tumor_type.is_unique
ranked = d.sort_values('LOY_Cox_p')
bh_sorted = np.minimum(1, np.minimum.accumulate(
    (ranked.LOY_Cox_p.to_numpy() * len(d) / np.arange(1, len(d)+1))[::-1])[::-1])
assert np.allclose(bh_sorted, ranked.LOY_q)
# The authors' notebook cells 21/23 obtain coefficient P values from
# lifelines cph.summary['p']['Y_coded']: two-sided Wald P, not a model LRT.
# S4 does not tabulate CI limits. Recover the same normal-Wald intervals
# algebraically from its full-precision HR and P; do not refit models.
d['log_hr'] = np.log(d.LOY_HR)
d['se_log_hr_recovered'] = np.abs(d.log_hr) / norm.isf(d.LOY_Cox_p / 2)
d['ci_lower'] = np.exp(d.log_hr - norm.ppf(.975) * d.se_log_hr_recovered)
d['ci_upper'] = np.exp(d.log_hr + norm.ppf(.975) * d.se_log_hr_recovered)
d['ci_source'] = 'Normal-Wald interval recovered from published HR and coefficient Wald P'
d['endpoint'] = 'PFS (TCGA-CDR ExtraEndpoints)'
d['model'] = 'LOY + age at diagnosis + autosomal chromosome arms lost'
d.to_csv(DATA / 'published_adjusted_estimates_with_recovered_ci.csv', index=False)
reference = pd.read_csv(OUT / 'source_data/qi_tcga_loy_prevalence/prevalence_by_tumour_type.csv')
display = reference[['Tumor_type', 'tumour_name']].merge(d, on='Tumor_type', how='left', validate='one_to_one')
assert display.Tumor_type.tolist() == reference.Tumor_type.tolist()
assert set(display.loc[display.LOY_HR.isna(), 'Tumor_type']) == {'LAML', 'PCPG'}
display['estimate_status'] = np.where(display.LOY_HR.notna(), 'Published adjusted estimate', 'Not reported in Table S4')
display.to_csv(DATA / 'forest_display_table.csv', index=False)
# These rounded numbers are printed in the published authors' notebook cell 12.
overall = dict(HR=1.17, ci_lower=1.06, ci_upper=1.28, n=4271, events=1956,
               model='Unadjusted LOY-only Cox model; pooled across cancer types',
               endpoint='PFS (TCGA-CDR ExtraEndpoints)',
               source='Authors notebook cell 12; HR also reported in article Fig. 4A',
               precision='HR and CI available to two decimal places in saved output',
               logrank_p_article='<0.00117')
(DATA / 'published_overall_estimate.json').write_text(json.dumps(overall, indent=2), encoding='utf-8')

fig, ax = plt.subplots(figsize=(183/25.4, 215/25.4))
fig.subplots_adjust(left=.475, right=.735, bottom=.12, top=.85)
ax.set_xscale('log')
ax.set_xlim(.025, 64)
ax.set_ylim(30.6, -1.2)
ticks = [.05, .25, 1, 4, 16, 64]
ax.xaxis.set_major_locator(FixedLocator(ticks))
ax.xaxis.set_major_formatter(FixedFormatter([str(x) for x in ticks]))
ax.xaxis.set_minor_locator(NullLocator())
ax.tick_params(axis='x', labelsize=7.5, length=3)
ax.set_yticks([])
ax.spines['left'].set_visible(False)
ax.axvline(1, color='#999999', lw=.7, ls=(0,(3,3)), zorder=0)
ax.set_xlabel('Hazard ratio (LOY versus WT)', fontsize=8.5, labelpad=8)
trans = ax.get_yaxis_transform()
def text_at(x, y, label, **kw):
    # x is in figure fractions and y is on the shared data row scale.
    ax.text((x-.475)/(.735-.475), y, label, transform=trans,
            va='center', clip_on=False, **kw)
def interval(y, hr, lo, hi, marker='o', colour=red):
    ax.errorbar(hr, y, xerr=[[hr-lo], [hi-hr]], fmt=marker, ms=4,
                color=colour, ecolor=colour, elinewidth=.85, capsize=2,
                markeredgewidth=.7, zorder=3)

text_at(.035, -.85, 'All cancers (unadjusted)', fontsize=8, weight='semibold')
text_at(.851, -.85, 'HR (95% CI)', ha='center', fontsize=8, weight='semibold')
text_at(.983, -.85, 'q', ha='right', fontsize=8.5, weight='semibold')
interval(.4, 1.17, 1.06, 1.28, marker='D', colour='#444444')
text_at(.035, .4, 'All tumours', fontsize=8, weight='semibold')
text_at(.851, .4, '1.17 (1.06–1.28)', fontsize=7.2, ha='center', weight='semibold')
text_at(.983, .4, '—', fontsize=8, ha='right', color='#666666')
text_at(.035, 2, 'Cancer-specific estimates (adjusted)', fontsize=8, weight='semibold')
for i, row in enumerate(display.itertuples()):
    y = i + 3
    sig = pd.notna(row.LOY_q) and row.LOY_q < .05
    text_at(.035, y, f'{row.tumour_name} ({row.Tumor_type})', fontsize=7.7,
            weight='semibold' if sig else 'normal')
    if pd.isna(row.LOY_HR):
        text_at(.851, y, 'Not reported', fontsize=7.2, ha='center', color='#777777')
        text_at(.983, y, '—', fontsize=7.5, ha='right', color='#777777')
        continue
    interval(y, row.LOY_HR, row.ci_lower, row.ci_upper)
    text_at(.851, y, f'{row.LOY_HR:.2f} ({row.ci_lower:.2f}–{row.ci_upper:.2f})',
            ha='center', fontsize=7.2, weight='semibold' if sig else 'normal')
    text_at(.983, y, f'{row.LOY_q:.3f}', fontsize=7.5, ha='right',
            weight='semibold' if sig else 'normal', color=red if sig else '#444444')

fig.text(.035, .975, 'b', fontsize=13, weight='semibold', va='top')
fig.text(.085, .975, 'Complete LOY and progression-free survival', fontsize=11, weight='semibold', va='top')
fig.text(.085, .938, 'Published estimates from Qi et al. (2023)', fontsize=9.5, va='top')
fig.text(.085, .909, 'Cancer-specific models: age and chromosome-arm loss adjustment', fontsize=8.5, va='top')
fig.text(.475, .058, 'Lower hazard', fontsize=7.5, color='#555555')
fig.text(.735, .058, 'Higher hazard', fontsize=7.5, color='#555555', ha='right')
fig.text(.035, .025, 'Source: Qi et al. (2023), Table S4 and published analysis notebook.', fontsize=7.3, color='#555555')
stem = 'ch03_main_intro_tcga_loy_pfs_forest'
for ext in ['pdf','svg','png']:
    fig.savefig(OUT / f'{stem}.{ext}', dpi=600)
fig.savefig(OUT / 'previews' / f'{stem}.png', dpi=180)
plt.close(fig)
meta = dict(figure=stem, endpoint=overall['endpoint'], font=font, size_mm=[183,215],
            source_url='https://ars.els-cdn.com/content/image/1-s2.0-S0092867423006463-mmc4.xlsx',
            notebook_url='https://github.com/rheinbaylab/LOY_in_primary_tumors/blob/master/Notebooks/LOY%20Survival%20analyses%20FINAL.ipynb',
            source_doi='10.1016/j.cell.2023.06.006', downloaded='2026-09-17',
            source_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),
            notebook_sha256=hashlib.sha256((DATA/'LOY_Survival_analyses_FINAL.ipynb').read_bytes()).hexdigest(),
            script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            model=d.model.iloc[0], adjusted_estimates=26, not_reported=['LAML','PCPG'],
            ci_method=d.ci_source.iloc[0], q_method='Published BH adjustment across 26 LOY coefficient Wald tests',
            overall=overall, interpretation='Published associations, not causal effects; no new models fitted.')
(DATA / 'render_provenance.json').write_text(json.dumps(meta, indent=2), encoding='utf-8')
print('Rendered published PFS estimates: 26 adjusted cancer types, 2 unreported rows, separate unadjusted overall estimate.')
