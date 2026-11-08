"""Focused, post hoc PSMD8 inspection; retain the original expression split."""
from pathlib import Path
import json
import sys
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pymupdf as fitz
from prepare_ddx3x_partner_adjustment import nuisance_design
from prepare_ddx3y_expression_dependency import bootstrap_r
from prepare_xy_ranked_partners_by_y_status import sha256

THESIS = Path(__file__).resolve().parents[3]
OUT = THESIS / 'figures/chapter4'
PREVIOUS = OUT / 'source_data/ddx3y_expression_dependency'
DATA = OUT / 'source_data/psmd8_ddx3y_followup'
STEM = 'ch04_psmd8_ddx3y_followup_draft'
sys.path.insert(0, str(THESIS / 'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style

GROUPS = ['XY_low_DDX3Y', 'XY_high_DDX3Y', 'LOY', 'Female']
LABELS = ['XY: lower DDX3Y', 'XY: higher DDX3Y', 'LOY', 'Female']


def residual_values(frame, method='primary'):
    z = nuisance_design(frame, 'raw' if method == 'raw' else 'lineage_efficacy').to_numpy(float)
    if method == 'add_DDX3Y_expression':
        z = np.column_stack([z, frame.DDX3Y_expression])
    if method == 'add_DDX3X_expression':
        z = np.column_stack([z, frame.DDX3X_expression])
    vals = frame[['DDX3X', 'PSMD8']].to_numpy(float)
    residual = vals - z @ np.linalg.lstsq(z, vals, rcond=None)[0]
    return residual, z, float(np.corrcoef(residual, rowvar=False)[0, 1])


def main():
    DATA.mkdir(exist_ok=True, parents=True)
    protocol = dict(
        date='2026-10-20', purpose='User-requested inspection of PSMD8 anticorrelation.',
        selection='Post hoc focus on a gene observed in the existing 31-partner heatmap.',
        fixed_input='Reuse the previous expression-covered models and median split. '
                    'Do not optimize threshold or remove outliers from primary results.',
        methods=['Raw Pearson', 'Primary adjusted Pearson',
                 'Spearman of primary OLS residuals (descriptive rank robustness)',
                 'Add DDX3Y expression to primary nuisance design',
                 'Add DDX3X expression to primary nuisance design'],
        intervals='1000 donor bootstrap resamples, refitting nuisance regressions; '
                  '95% percentile intervals, fixed seed 20261020.',
        influence='Refit primary adjustment after removing each donor in turn and '
                  'each lineage in turn. All removal results retained; no primary exclusions.',
        interpretation='Preserve the original continuous interaction P values and '
                       '31-test BH FDRs. No new significance claim from selected '
                       'stratum intervals, and no physical/suppressor interaction inference.')
    path = DATA / 'protocol.json'
    if path.exists():
        assert json.loads(path.read_text()) == protocol
    else:
        path.write_text(json.dumps(protocol, indent=2), encoding='utf-8')
    source = PREVIOUS / 'analysis_model_values.csv'
    frame = pd.read_csv(source, index_col='ModelID')
    saved = pd.read_csv(PREVIOUS / 'stratum_correlations.csv')
    saved = saved.loc[saved.gene.eq('PSMD8')].set_index('stratum')
    all_rows, all_residuals, influence, bootstrap = [], [], [], []
    for group in GROUPS:
        sub = frame.loc[frame.display_group.eq(group)]
        for method in ['raw', 'primary', 'add_DDX3Y_expression', 'add_DDX3X_expression']:
            residual, design, r = residual_values(sub, method)
            row = dict(group=group, analysis=method, n_models=len(sub),
                       n_patients=sub.PatientID.nunique(), pearson_r=r,
                       spearman_of_residuals=stats.spearmanr(residual[:, 0], residual[:, 1]).statistic)
            if method == 'primary':
                np.testing.assert_allclose(r, saved.loc[group, 'pearson_r'], atol=1e-12)
                boots = bootstrap_r(sub.PSMD8, sub.DDX3X, design, sub.PatientID)
                row['r_ci_low'], row['r_ci_high'] = np.quantile(boots, [.025, .975])
                bootstrap.extend(dict(group=group, iteration=i, pearson_r=b)
                                 for i, b in enumerate(boots))
                all_residuals.append(pd.DataFrame(dict(
                    ModelID=sub.index, group=group, CellLineName=sub.CellLineName,
                    PatientID=sub.PatientID, OncotreeLineage=sub.OncotreeLineage,
                    ddx3x_residual=residual[:, 0], psmd8_residual=residual[:, 1])))
                for field, kind in [('PatientID', 'donor'), ('OncotreeLineage', 'lineage')]:
                    for value in sub[field].unique():
                        rest = sub.loc[sub[field].ne(value)]
                        _, _, leave_r = residual_values(rest)
                        influence.append(dict(
                            group=group, removed_kind=kind, removed_value=value,
                            n_removed=len(sub) - len(rest), n_remaining=len(rest),
                            pearson_r=leave_r, delta_r=leave_r - r,
                            removed_models='; '.join(sub.loc[sub[field].eq(value), 'CellLineName'])))
            all_rows.append(row)
        print(group, 'complete', flush=True)
    stats_table = pd.DataFrame(all_rows)
    residual_table = pd.concat(all_residuals, ignore_index=True)
    influence_table = pd.DataFrame(influence)
    stats_table.to_csv(DATA / 'correlations_and_sensitivities.csv', index=False)
    residual_table.to_csv(DATA / 'plotted_residuals.csv', index=False)
    influence_table.to_csv(DATA / 'leave_one_out_influence.csv', index=False)
    pd.DataFrame(bootstrap).to_csv(DATA / 'bootstrap_correlations.csv.gz', index=False,
                                 compression={'method': 'gzip', 'mtime': 0})
    interactions = pd.read_csv(PREVIOUS / 'expression_partner_interactions.csv')
    interactions.loc[interactions.gene.eq('PSMD8')].to_csv(
        DATA / 'previous_psmd8_interaction_tests.csv', index=False)
    ranges = influence_table.groupby(['group', 'removed_kind']).pearson_r.agg(['min', 'max'])
    ranges.to_csv(DATA / 'influence_ranges.csv')
    manifest = dict(
        **protocol, sources={str(p): sha256(p) for p in [
            source, PREVIOUS / 'stratum_correlations.csv',
            PREVIOUS / 'expression_partner_interactions.csv', path, Path(__file__)]},
        primary_matches_saved=True,
        source_release='DepMap Public 25Q3; prior cohort and Chapter 3 eligibility unchanged.')
    (DATA / 'manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    print(stats_table.to_string(index=False))
    print(ranges.to_string())
    render(stats_table, residual_table)


def render(statistics, residuals):
    apply_style(9)
    fig, axes = plt.subplots(2, 2, figsize=(8.5, 7.1), sharex=True, sharey=True)
    fig.subplots_adjust(left=.105, right=.975, bottom=.115, top=.925,
                        wspace=.15, hspace=.32)
    colors = ['#3C5488', '#3C5488', '#E64B35', '#8A6D9E']
    x = residuals.ddx3x_residual
    y = residuals.psmd8_residual
    xlim = (float(x.min() - .1), float(x.max() + .1))
    ylim = (float(y.min() - .08), float(y.max() + .12))
    for i, (group, label, color) in enumerate(zip(GROUPS, LABELS, colors)):
        ax = axes.flat[i]
        values = residuals.loc[residuals.group.eq(group)]
        row = statistics.loc[statistics.group.eq(group) & statistics.analysis.eq('primary')].iloc[0]
        x, y = values.ddx3x_residual.to_numpy(), values.psmd8_residual.to_numpy()
        ax.scatter(x, y, color=color, s=14, alpha=.48, edgecolors='none')
        grid = np.linspace(x.min(), x.max(), 100)
        slope = np.dot(x, y) / np.dot(x, x)
        ax.plot(grid, slope * grid, color='#222222', lw=1.3)
        ax.axhline(0, color='#BBBBBB', ls='--', lw=.5, zorder=0)
        ax.axvline(0, color='#BBBBBB', ls='--', lw=.5, zorder=0)
        ax.set(xlim=xlim, ylim=ylim)
        ax.text(.03, .96,
                f'r = {row.pearson_r:.2f} [{row.r_ci_low:.2f}, {row.r_ci_high:.2f}]\n'
                f'n = {int(row.n_models)} models',
                transform=ax.transAxes, va='top', fontsize=8.5,
                bbox=dict(facecolor='white', alpha=.90, edgecolor='none', pad=2))
        ax.text(0, 1.045, label, transform=ax.transAxes,
                fontsize=9.5, fontweight='semibold')
        ax.text(-.18, 1.045, chr(ord('a') + i), transform=ax.transAxes,
                fontsize=12, fontweight='semibold')
        if i // 2 == 1:
            ax.set_xlabel('DDX3X gene-effect residual')
        if i % 2 == 0:
            ax.set_ylabel('PSMD8 gene-effect residual')
    fig.text(.54, .032, 'Adjusted for lineage, screen efficacy and library/platform',
             ha='center', fontsize=8.5, color='#555555')
    for ext in ['pdf', 'svg', 'png']:
        fig.savefig(OUT / f'{STEM}.{ext}', dpi=600 if ext == 'png' else 300)
    fig.savefig(OUT / 'review' / f'{STEM}_preview.png', dpi=160)
    plt.close(fig)
    with fitz.open(OUT / f'{STEM}.pdf') as doc:
        page = doc[0]
        assert all('OpenSans' in f[3] for f in page.get_fonts())
        spans = [s for b in page.get_text('dict')['blocks'] if 'lines' in b
                 for line in b['lines'] for s in line['spans'] if s['text'].strip()]
        outside = [s['text'] for s in spans if not page.rect.contains(fitz.Rect(s['bbox']))]
        overlaps = []
        for i, a in enumerate(spans):
            for b in spans[i + 1:]:
                ra, rb = fitz.Rect(a['bbox']), fitz.Rect(b['bbox'])
                area = (ra & rb).get_area()
                if area > 1 and area / min(ra.get_area(), rb.get_area()) > .1:
                    overlaps.append([a['text'], b['text']])
        assert not outside, outside
        assert not overlaps, overlaps
        page.get_pixmap(matrix=fitz.Matrix(1.5, 1.5), colorspace=fitz.csGRAY).save(
            OUT / 'review' / f'{STEM}_greyscale.png')
    (DATA / 'render_manifest.json').write_text(json.dumps(dict(
        stem=STEM, dimensions_inches=[8.5, 7.1], font='Open Sans',
        x_limits=xlim, y_limits=ylim, group_order=GROUPS,
        interval='Nominal 95% donor-bootstrap percentile interval; post hoc selected gene.',
        text_outside_page=outside, text_overlaps=overlaps), indent=2), encoding='utf-8')


if __name__ == '__main__':
    main()
