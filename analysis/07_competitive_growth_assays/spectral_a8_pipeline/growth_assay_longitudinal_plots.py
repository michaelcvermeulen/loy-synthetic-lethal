"""Render D0/D4 trajectories from the approved report, without reprocessing events.

Run with a Python environment containing Matplotlib:
    python -m spectral.a8_pipeline.growth_assay_longitudinal_plots --report-dir <report-directory>
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator

from .growth_assay_replicate_plots import render_replicates

REPO = Path(__file__).resolve().parents[2]
DEFAULT_REPORT = REPO / 'analysis/spectral_growth_assay_D4_20261020'
GROUPS = ('WT', 'WT_C1', 'DDX3X', 'DDX3X_C1', 'SREBF1', 'TADA2B')
COLORS = {'LOY_GFP': '#1f7a4d', 'LOY_mCherry': '#a8329b'}
TEMPLATE = Path(__file__).with_name('growth_assay_followup.html')


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def render(report: Path) -> None:
    report = report.resolve()
    results_path = report / 'results.json'
    results_hash = digest(results_path)
    data = json.loads(results_path.read_text(encoding='utf-8'))
    if data['identity_status'] != 'user_confirmed' or data['primary_negative_rule'] != 'broad':
        raise ValueError('Expected approved D4 identities and broad-negative analysis')
    wells = data['wells']
    summaries = {(s['group'], s['orientation']): s for s in data['orientation_summaries']}
    figures = report / 'figures'
    figures.mkdir(exist_ok=True)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                         'axes.titlesize': 13, 'axes.labelsize': 10,
                         'svg.fonttype': 'none', 'pdf.fonttype': 42,
                         'axes.spines.top': False, 'axes.spines.right': False})
    output_paths = []
    values = [w['delta_log_odds'] for w in wells if w['role'] == 'mixture' and w['delta_log_odds'] is not None]
    import math
    delta_limits = (min(-.5, math.floor(min(values)*2)/2), max(.5, math.ceil(max(values)*2)/2))
    for mode in ('fraction', 'change'):
        normalized = mode == 'change'
        fig, axes = plt.subplots(2, 3, figsize=(12, 8), sharex=True, sharey=True)
        fig.subplots_adjust(left=.075, right=.97, top=.81, bottom=.15, hspace=.28, wspace=.18)
        fig.suptitle('D0 to D4 · ' + ('change in LOY log-odds' if normalized else 'LOY fraction'),
                     x=.075, y=.98, ha='left', fontsize=21, fontweight='bold', color='#172d24')
        fig.text(.075, .922, 'Corrected sample labels · broad-negative exclusion · one panel per assay arm',
                 fontsize=11, color='#52635b')
        legend = [Line2D([], [], color=COLORS[o], lw=2, linestyle='-' if i == 0 else '--',
                         label='LOY = GFP' if i == 0 else 'LOY = mCherry') for i, o in enumerate(COLORS)]
        legend += [Line2D([], [], color='#52635b', marker='D', linestyle='none', markersize=5,
                          label='Pooled D0 baseline'),
                   Line2D([], [], color='#52635b', marker='D', markerfacecolor='white', linestyle='none',
                          markersize=6, label='D4 mean'),
                   Line2D([], [], color='#52635b', marker='o', linestyle='none', markersize=4,
                          label='Individual D4 wells')]
        fig.legend(handles=legend, loc='upper left', bbox_to_anchor=(.067, .9), ncol=5,
                   frameon=False, fontsize=9, handlelength=2)
        for ax, group in zip(axes.flat, GROUPS):
            ax.set_title(group.replace('_', ' '), loc='left', color='#172d24', fontweight='bold')
            ax.set_xlim(-.55, 4.8)
            ax.set_ylim(*(delta_limits if normalized else (0, 100)))
            ax.set_xticks([0, 4], ['D0', 'D4'])
            ax.tick_params(axis='x', labelbottom=True)
            ax.yaxis.set_major_locator(MultipleLocator(.5 if normalized else 20))
            ax.set_axisbelow(True)
            ax.grid(axis='y', color='#dfe7e1', linewidth=.7)
            ax.spines[['left', 'bottom']].set_color('#a9b9af')
            ax.tick_params(colors='#52635b', length=3)
            if normalized:
                ax.axhline(0, color='#84968b', linewidth=.8, linestyle=':')
            for j, (orientation, color) in enumerate(COLORS.items()):
                s = summaries[group, orientation]
                ws = sorted((w for w in wells if w['role'] == 'mixture' and w['group'] == group
                             and w['orientation'] == orientation), key=lambda w: w['replicate'])
                if len(ws) != 3:
                    raise ValueError(f'Expected three D4 wells: {group} {orientation}')
                p0 = 0 if normalized else 100*s['d0_loy_fraction']
                mean = s['delta_log_odds_mean'] if normalized else 100*s['loy_fraction_mean']
                x0 = (-.035 if j == 0 else .035) if normalized else 0
                ax.plot([x0, 4], [p0, mean], color=color, lw=1.8, linestyle='-' if j == 0 else '--')
                ax.plot(x0, p0, marker='D', color=color, markersize=5, linestyle='none', zorder=4)
                key = 'delta_log_odds' if normalized else 'loy_fraction'
                scale = 1 if normalized else 100
                for w in ws:
                    if w[key] is not None:
                        ax.scatter(4+(w['replicate']-2)*.13, scale*w[key], s=23, color=color,
                                   edgecolors='white', linewidths=.5, zorder=5)
                ax.plot(4, mean, marker='D', color=color, markerfacecolor='white', markersize=6,
                        markeredgewidth=1.4, linestyle='none', zorder=6)
        for ax in axes[:, 0]:
            ax.set_ylabel('Change in LOY log-odds from D0' if normalized else 'LOY among retained cells (%)')
        fig.text(.075, .072, 'D0: pooled counts from the original mixing batch. D4: equally weighted mean of three wells; dots are offset for visibility.',
                 fontsize=9, color='#52635b')
        fig.text(.075, .047, 'Lines connect the two measured timepoints. DDX3X arms share D0; WT arms use separate D0 mixtures. No independent D0–D4 pairing is implied.',
                 fontsize=8.5, color='#52635b')
        stem = 'longitudinal_log_odds' if normalized else 'longitudinal_loy_fraction'
        for suffix in ('png', 'svg', 'pdf'):
            path = figures / f'{stem}.{suffix}'
            fig.savefig(path, dpi=180, facecolor='white')
            output_paths.append(path)
        plt.close(fig)

    rows = []
    for (group, orientation), s in summaries.items():
        rows.append(dict(group=group, orientation=orientation, day=0, kind='pooled_baseline',
                         sample='', replicate='', loy_fraction=s['d0_loy_fraction'], delta_log_odds=0))
        rows.append(dict(group=group, orientation=orientation, day=4, kind='replicate_mean',
                         sample='', replicate='', loy_fraction=s['loy_fraction_mean'], delta_log_odds=s['delta_log_odds_mean']))
    for w in wells:
        if w['role'] == 'mixture':
            rows.append(dict(group=w['group'], orientation=w['orientation'], day=4, kind='individual_well',
                             sample=w['analysis_sample'], replicate=w['replicate'], loy_fraction=w['loy_fraction'],
                             delta_log_odds=w['delta_log_odds']))
    csv_path = figures / 'longitudinal_data.csv'
    with csv_path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    output_paths.append(csv_path)
    output_paths.extend(render_replicates(data, figures))
    if digest(results_path) != results_hash:
        raise ValueError('Analysis changed during plotting; HTML not updated')
    payload = json.dumps(data, ensure_ascii=False, allow_nan=False).replace('<', '\\u003c')
    html = TEMPLATE.read_text(encoding='utf-8').replace('__DATA__', payload)
    html = html.replace('id="trajectory-exports" hidden', 'id="trajectory-exports"')
    html = html.replace('id="replicate-exports" hidden', 'id="replicate-exports"')
    if (report / 'same_genotype/report.html').exists():
        html = html.replace('id="same-genotype-link" hidden', 'id="same-genotype-link"')
    temporary = report / 'report.html.updating'
    temporary.write_text(html, encoding='utf-8')
    temporary.replace(report / 'report.html')
    manifest_path = report / 'code_manifest.json'
    manifest = json.loads(manifest_path.read_text(encoding='utf-8'))
    for p in (Path(__file__), TEMPLATE, Path(__file__).with_name('growth_assay_replicate_plots.py')):
        manifest[str(p.relative_to(REPO))] = digest(p)
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    provenance = dict(generated_utc=datetime.now(timezone.utc).isoformat(), results_sha256=results_hash,
                      matplotlib_version=matplotlib.__version__, script_sha256=digest(Path(__file__)),
                      files={str(p.relative_to(report)): digest(p) for p in output_paths})
    (figures / 'longitudinal_manifest.json').write_text(json.dumps(provenance, indent=2), encoding='utf-8')
    print(json.dumps(dict(report=str(report / 'report.html'), figures=[str(p) for p in output_paths],
                         measurements_unchanged=True), indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--report-dir', type=Path, default=DEFAULT_REPORT)
    render(parser.parse_args().report_dir)
