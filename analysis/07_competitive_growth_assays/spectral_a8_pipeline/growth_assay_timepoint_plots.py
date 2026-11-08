"""Render D0 → previous → current trajectories from a growth_assay_timepoint report, without reprocessing events.

Run with a Python environment containing Matplotlib:
    python -m spectral.a8_pipeline.growth_assay_timepoint_plots --report-dir <report-directory>
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
from hashlib import sha256
import json
import math
from pathlib import Path
import textwrap

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator

REPO = Path(__file__).resolve().parents[2]
PAIRS = (('WT', 'WT_C1'), ('DDX3X', 'DDX3X_C1'), ('SREBF1', 'TADA2B'))
SETS = ('WT', 'WT_C1', 'DDX3X', 'DDX3X_C1')  # backgrounds with two replicate sets
COLORS = {'LOY_GFP': '#1f7a4d', 'LOY_mCherry': '#a8329b', 'same_genotype': '#6b7280'}
STYLES = {'LOY_GFP': '-', 'LOY_mCherry': '--'}
OFFSET = {'LOY_GFP': -.3, 'LOY_mCherry': .3}
INK, MUTED, GRID = '#172d24', '#52635b', '#dfe7e1'


def panel_title(group: str) -> str:
    """C1-labelled arms are untreated second replicate sets (no drug was added)."""
    background, _, label = group.partition('_')
    return f'{background} · second set ({label} label)' if label else f'{background} · first set'


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def logit(p: float | None) -> float | None:
    return math.log(p / (1 - p)) if p is not None and 0 < p < 1 else None


def day(label: str) -> int:
    return int(label[1:])


def flagged(current: dict) -> set[str]:
    return {w['well'] for w in current['wells'] if w['qc_flags']}


def series(current: dict, previous: dict) -> list[dict]:
    """Every plotted competition value: pooled D0 mixes, earlier wells and current biological replicates."""
    rows, bad = [], flagged(current)
    for s in current['orientation_summaries']:
        key = dict(group=s['group'], orientation=s['orientation'])
        rows.append(dict(**key, day=0, kind='pooled_d0', replicate='', wells=s['d0_wells'], flagged=False,
                         fraction=s['d0_loy_fraction'], delta_log_odds=0.0))
        for w in previous['wells']:
            if w['role'] == 'mixture' and (w['group'], w['orientation']) == (s['group'], s['orientation']):
                rows.append(dict(**key, day=day(previous['timepoint']), kind='well', replicate=w['replicate'],
                                 wells=w['well'], flagged=False, fraction=w['loy_fraction'], delta_log_odds=w['delta_log_odds']))
        rows.append(dict(**key, day=day(previous['timepoint']), kind='arm_mean', replicate='', wells='', flagged=False,
                         fraction=s['previous_loy_fraction_mean'], delta_log_odds=s['previous_delta_log_odds_mean']))
        for b in current['biological_replicates']:
            if (b['group'], b['orientation']) == (s['group'], s['orientation']):
                rows.append(dict(**key, day=day(current['timepoint']), kind='biological_replicate', replicate=b['replicate'],
                                 wells=b['wells'], flagged=bool(set(b['wells'].split(';')) & bad),
                                 fraction=b['loy_fraction'], delta_log_odds=b['delta_log_odds']))
        rows.append(dict(**key, day=day(current['timepoint']), kind='arm_mean', replicate='', wells='', flagged=False,
                         fraction=s['loy_fraction_mean'], delta_log_odds=s['delta_log_odds_mean']))
    return rows


def control_series(current: dict, previous: dict) -> list[dict]:
    """Same-genotype GFP fractions: pooled D0, earlier wells and current biological replicates."""
    rows, bad = [], flagged(current)
    before = {c['genotype']: c for c in previous['dye_controls']}
    for c in current['dye_controls']:
        g, p0 = c['genotype'], c['d0_gfp_fraction_pooled']
        if not math.isclose(p0, before[g]['d0_gfp_fraction_pooled'], abs_tol=1e-12):
            raise ValueError('Same-genotype D0 references differ between reports')
        rows.append(dict(genotype=g, day=0, kind='pooled_d0', replicate='', wells=c['d0_wells'], flagged=False, fraction=p0))
        rows += [dict(genotype=g, day=day(previous['timepoint']), kind='well', replicate=w['replicate'], wells=w['well'],
                      flagged=False, fraction=w['gfp_fraction'])
                 for w in previous['wells'] if w['role'] == 'dye_neutral' and w['genotype'] == g]
        rows.append(dict(genotype=g, day=day(previous['timepoint']), kind='arm_mean', replicate='', wells='', flagged=False,
                         fraction=c['previous_gfp_fraction_mean']))
        rows += [dict(genotype=g, day=day(current['timepoint']), kind='biological_replicate', replicate=b['replicate'],
                      wells=b['wells'], flagged=bool(set(b['wells'].split(';')) & bad), fraction=b['gfp_fraction'])
                 for b in current['biological_replicates'] if b['group'] == 'same_genotype' and b['orientation'] == g]
        rows.append(dict(genotype=g, day=day(current['timepoint']), kind='arm_mean', replicate='', wells='', flagged=False,
                         fraction=c['gfp_fraction_mean']))
    for r in rows:
        r['gfp_delta_log_odds'] = logit(r['fraction']) - logit(next(x['fraction'] for x in rows if
                                                                  x['genotype'] == r['genotype'] and x['kind'] == 'pooled_d0'))
    return rows


def style_axis(ax, days: list[int], limits, step: float) -> None:
    ax.set_xlim(-1.2, days[-1] + 1.6)
    ax.set_ylim(*limits)
    ax.set_xticks(days, [f'D{d}' for d in days])
    ax.tick_params(axis='x', labelbottom=True)
    ax.yaxis.set_major_locator(MultipleLocator(step))
    ax.set_axisbelow(True)
    ax.grid(axis='y', color=GRID, linewidth=.7)
    ax.spines[['left', 'bottom']].set_color('#a9b9af')
    ax.tick_params(colors=MUTED, length=3)


def footnote(fig, x: float, top: float, notes: list[str], width: int) -> None:
    """Wrapped caption lines, top-aligned, so long notes never run off the figure."""
    text = '\n'.join(textwrap.fill(n, width) for n in notes)
    fig.text(x, top, text, fontsize=8.8, color=MUTED, va='top', linespacing=1.55)


def save(fig, stem: Path, paths: list[Path]) -> None:
    for suffix in ('png', 'svg', 'pdf'):
        path = stem.with_suffix('.' + suffix)
        fig.savefig(path, dpi=180, facecolor='white')
        paths.append(path)
    plt.close(fig)


def trajectories(current: dict, previous: dict, rows: list[dict], figures: Path, paths: list[Path]) -> None:
    tp, prev = current['timepoint'], previous['timepoint']
    days = [0, day(prev), day(tp)]
    deltas = [r['delta_log_odds'] for r in rows if r['delta_log_odds'] is not None]
    delta_limits = (min(-.5, math.floor(min(deltas) * 2) / 2), max(.5, math.ceil(max(deltas) * 2) / 2))
    for mode in ('fraction', 'log_odds'):
        log_odds = mode == 'log_odds'
        key, scale = ('delta_log_odds', 1) if log_odds else ('fraction', 100)
        fig, axes = plt.subplots(3, 2, figsize=(11, 12.4), sharex=True, sharey=True)
        fig.subplots_adjust(left=.085, right=.97, top=.865, bottom=.125, hspace=.3, wspace=.12)
        fig.suptitle(f'D0 to {tp} · ' + ('change in LOY log-odds' if log_odds else 'LOY fraction'),
                     x=.085, y=.98, ha='left', fontsize=20, fontweight='bold', color=INK)
        fig.text(.085, .943, f'{current["identity_status"].replace("_", " ")} · primary gates · no drug in any arm · '
                             'first replicate set left, second set right', fontsize=11, color=MUTED)
        legend = [Line2D([], [], color=COLORS[o], lw=2, linestyle=STYLES[o],
                         label='LOY = GFP' if o == 'LOY_GFP' else 'LOY = mCherry') for o in STYLES]
        legend += [Line2D([], [], color=MUTED, marker='D', linestyle='none', markersize=5, label='Pooled D0 mix'),
                   Line2D([], [], color=MUTED, marker='o', linestyle='none', markersize=3.5, label=f'{prev} wells'),
                   Line2D([], [], color=MUTED, marker='o', linestyle='none', markersize=6,
                          label=f'{tp} biological replicates (mean of 2 wells)'),
                   Line2D([], [], color=MUTED, marker='o', markerfacecolor='white', linestyle='none', markersize=6,
                          label='Replicate with a QC-flagged well'),
                   Line2D([], [], color=MUTED, marker='D', markerfacecolor='white', linestyle='none', markersize=7,
                          label='Arm mean')]
        fig.legend(handles=legend, loc='upper left', bbox_to_anchor=(.077, .93), ncol=4, frameon=False, fontsize=9,
                   handlelength=2.2)
        for ax, group in zip(axes.flat, [g for pair in PAIRS for g in pair]):
            ax.set_title(panel_title(group) if group in SETS else group, loc='left', color=INK, fontweight='bold',
                         fontsize=13)
            style_axis(ax, days, delta_limits if log_odds else (0, 100), .5 if log_odds else 20)
            if log_odds:
                ax.axhline(0, color='#a9b9af', linewidth=.9)
            for orientation, color in ((o, COLORS[o]) for o in STYLES):
                mine = [r for r in rows if (r['group'], r['orientation']) == (group, orientation)]
                means = sorted((r for r in mine if r['kind'] in ('pooled_d0', 'arm_mean')), key=lambda r: r['day'])
                x = [r['day'] + (OFFSET[orientation] * .25 if r['day'] == 0 else 0) for r in means]
                ax.plot(x, [scale * r[key] for r in means], color=color, lw=1.8, linestyle=STYLES[orientation], zorder=3)
                ax.plot(x[0], scale * means[0][key], marker='D', color=color, markersize=5, linestyle='none', zorder=5)
                for r in means[1:]:
                    ax.plot(r['day'], scale * r[key], marker='D', color=color, markerfacecolor='white', markersize=7,
                            markeredgewidth=1.4, linestyle='none', zorder=6)
                for r in mine:
                    if r['kind'] not in ('well', 'biological_replicate') or r[key] is None:
                        continue
                    big = r['kind'] == 'biological_replicate'
                    ax.plot(r['day'] + OFFSET[orientation] + (r['replicate'] - 2) * .17, scale * r[key], marker='o',
                            markersize=5.5 if big else 3.5, color=color, markerfacecolor='white' if r['flagged'] else color,
                            markeredgewidth=1.1 if r['flagged'] else .5, markeredgecolor=color if r['flagged'] else 'white',
                            linestyle='none', zorder=5, alpha=1 if big else .75)
        for ax in axes[:, 0]:
            ax.set_ylabel('Change in LOY log-odds from D0' if log_odds else 'LOY among retained cells (%)', color=INK)
        footnote(fig, .085, .075, [
            f'D0: pooled counts of the original mixing batch. {prev}: individual wells and their mean. {tp}: two technical '
            'wells averaged per biological replicate; arm mean across the three replicates.',
            f'{tp} was acquired {current["elapsed_days_from_d0"]} days after the D0 plate and is drawn at its nominal day. '
            'Replicate k is linked to earlier replicate k by number only. Open circles: a well of that replicate carries '
            'a QC flag (see report).',
            'No arm received drug: the C1-labelled arms are second replicate sets. WT C1 was a separate D0 mix; DDX3X C1 '
            'was dispensed from the DDX3X mix. TADA2B uses the day-specific unlabelled-cell gate; other arms the standard gate.'],
            width=150)
        save(fig, figures / f'trajectories_{mode}', paths)


def controls(current: dict, previous: dict, rows: list[dict], figures: Path, paths: list[Path]) -> None:
    tp, prev = current['timepoint'], previous['timepoint']
    days = [0, day(prev), day(tp)]
    fig, axes = plt.subplots(1, 2, figsize=(10, 5.2), sharey=True)
    fig.subplots_adjust(left=.085, right=.97, top=.78, bottom=.2, wspace=.12)
    fig.suptitle(f'Same-genotype mixtures · GFP fraction, D0 to {tp}', x=.085, y=.97, ha='left', fontsize=17,
                 fontweight='bold', color=INK)
    fig.text(.085, .885, 'WT background · equal GFP and mCherry sublines of one genotype · 50% = no subline difference',
             fontsize=10.5, color=MUTED)
    for ax, genotype in zip(axes, ('LOY', 'CTRL')):
        ax.set_title(f'{genotype} / {genotype}', loc='left', color=INK, fontweight='bold', fontsize=13)
        style_axis(ax, days, (20, 80), 10)
        ax.axhline(50, color='#a9b9af', linewidth=.9)
        mine = [r for r in rows if r['genotype'] == genotype]
        means = sorted((r for r in mine if r['kind'] in ('pooled_d0', 'arm_mean')), key=lambda r: r['day'])
        ax.plot([r['day'] for r in means], [100 * r['fraction'] for r in means], color=COLORS['same_genotype'], lw=1.8)
        ax.plot(0, 100 * means[0]['fraction'], marker='D', color=COLORS['same_genotype'], markersize=5, zorder=5)
        for r in means[1:]:
            ax.plot(r['day'], 100 * r['fraction'], marker='D', color=COLORS['same_genotype'], markerfacecolor='white',
                    markersize=7, markeredgewidth=1.4, zorder=6)
        for r in mine:
            if r['kind'] in ('well', 'biological_replicate'):
                big = r['kind'] == 'biological_replicate'
                ax.plot(r['day'] + (r['replicate'] - 2) * .25, 100 * r['fraction'], marker='o',
                        markersize=5.5 if big else 3.5, color=COLORS['same_genotype'], markeredgecolor='white',
                        markeredgewidth=.5, linestyle='none', zorder=5)
        end = means[-1]
        ax.annotate(f"{100 * end['fraction']:.1f}%", (end['day'], 100 * end['fraction']), xytext=(9, 0),
                    textcoords='offset points', va='center', fontsize=9, color=INK)
    axes[0].set_ylabel('GFP among retained cells (%)', color=INK)
    fig.text(.085, .075, f'Diamonds: pooled D0 and arm means. Small dots: {prev} wells. Larger dots: {tp} biological '
                         'replicates (mean of two technical wells).', fontsize=8.8, color=MUTED)
    fig.text(.085, .045, 'A shift away from 50% means one reporter subline outgrew the other within the same genotype.',
             fontsize=8.8, color=MUTED)
    save(fig, figures / 'same_genotype_gfp', paths)


def technical(current: dict, figures: Path, paths: list[Path]) -> list[dict]:
    pairs = []
    for b in current['biological_replicates']:
        wells = {w['well']: w for w in current['wells']}
        first, second = (wells[x] for x in b['wells'].split(';'))
        readout = 'gfp_fraction' if b['group'] == 'same_genotype' else 'loy_fraction'
        pairs.append(dict(group=b['group'], orientation=b['orientation'], replicate=b['replicate'], wells=b['wells'],
                          readout=readout, technical_1=first[readout], technical_2=second[readout],
                          gap_pp=b['technical_gap_pp']))
    limit = current['flag_thresholds']['technical_gap_pp']
    fig, ax = plt.subplots(figsize=(7.2, 7.4))
    fig.subplots_adjust(left=.13, right=.96, top=.84, bottom=.13)
    fig.suptitle(f'{current["timepoint"]} technical replicates', x=.13, y=.97, ha='left', fontsize=17,
                 fontweight='bold', color=INK)
    fig.text(.13, .915, 'Each point is one biological replicate: its two wells, in column order', fontsize=10.5, color=MUTED)
    lo, hi = 10, 70
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    ax.set_aspect('equal')
    ax.grid(color=GRID, linewidth=.7)
    ax.set_axisbelow(True)
    ax.plot([lo, hi], [lo, hi], color='#a9b9af', lw=1)
    for d in (-limit, limit):
        ax.plot([lo, hi], [lo + d, hi + d], color='#a9b9af', lw=.8, linestyle=(0, (4, 3)))
    for p in pairs:
        kind = 'same_genotype' if p['group'] == 'same_genotype' else p['orientation']
        ax.plot(100 * p['technical_1'], 100 * p['technical_2'], marker='s' if kind == 'same_genotype' else 'o',
                markersize=6.5, color=COLORS[kind], markeredgecolor='white', markeredgewidth=.6, linestyle='none')
        if p['gap_pp'] > limit:
            ax.annotate(p['wells'].replace(';', '/'), (100 * p['technical_1'], 100 * p['technical_2']), xytext=(-8, 4),
                        textcoords='offset points', ha='right', fontsize=8.5, color=INK)
    ax.set_xlabel('First well (%)', color=INK)
    ax.set_ylabel('Second well (%)', color=INK)
    ax.tick_params(colors=MUTED)
    ax.spines[['left', 'bottom']].set_color('#a9b9af')
    ax.legend(handles=[Line2D([], [], color=COLORS['LOY_GFP'], marker='o', linestyle='none', label='LOY % · LOY = GFP'),
                       Line2D([], [], color=COLORS['LOY_mCherry'], marker='o', linestyle='none', label='LOY % · LOY = mCherry'),
                       Line2D([], [], color=COLORS['same_genotype'], marker='s', linestyle='none',
                              label='GFP % · same-genotype'),
                       Line2D([], [], color='#a9b9af', linestyle=(0, (4, 3)), label=f'±{limit:g} pp')],
              loc='upper left', frameon=False, fontsize=9)
    fig.text(.13, .035, f'Labelled pairs differ by more than {limit:g} percentage points.', fontsize=8.8, color=MUTED)
    save(fig, figures / 'technical_replicates', paths)
    return pairs


def identity(current: dict, figures: Path, paths: list[Path]) -> None:
    checks = current['identity_checks']
    fig, (left, right) = plt.subplots(1, 2, figsize=(12, 5.4), gridspec_kw=dict(width_ratios=[1.35, 1]))
    fig.subplots_adjust(left=.07, right=.98, top=.78, bottom=.17, wspace=.22)
    fig.suptitle(f'{current["timepoint"]} identity checks · WT background', x=.07, y=.97, ha='left', fontsize=17,
                 fontweight='bold', color=INK)
    fig.text(.07, .885, 'LOY cells are larger on LightLoss (Violet)-A than CTRL cells in WT backgrounds', fontsize=10.5,
             color=MUTED)
    size = [r for r in checks['size_fingerprint'] if r['role'] == 'mixture']
    for i, r in enumerate(size):
        left.plot(i, r['gfp_to_mcherry_ratio'], marker='o', markersize=6.5, color=COLORS[r['orientation']],
                  markeredgecolor='white', markeredgewidth=.6, linestyle='none')
    left.axhline(1, color='#a9b9af', lw=1)
    left.set_xticks(range(len(size)), [r['well'] for r in size], rotation=90, fontsize=8)
    left.set_ylabel('GFP-class / mCherry-class median LightLoss', color=INK)
    left.set_title('Competition wells: recorded LOY reporter on the larger cells', loc='left', fontsize=11, color=INK,
                   fontweight='bold', pad=20)
    groups = [r['group'] for r in size]
    for group in dict.fromkeys(groups):
        first, last = groups.index(group), len(groups) - 1 - groups[::-1].index(group)
        left.text((first + last) / 2, 1.015, group.replace('_', ' '), transform=left.get_xaxis_transform(), ha='center',
                  va='bottom', fontsize=9, color=MUTED)
        if first:
            left.axvline(first - .5, color=GRID, lw=1)
    left.legend(handles=[Line2D([], [], color=COLORS[o], marker='o', linestyle='none',
                                label='recorded LOY = GFP (expect > 1)' if o == 'LOY_GFP' else 'recorded LOY = mCherry (expect < 1)')
                         for o in ('LOY_GFP', 'LOY_mCherry')], loc='center left', bbox_to_anchor=(0, .58),
                frameon=False, fontsize=9)
    genotype = checks['same_genotype']
    for i, r in enumerate(genotype):
        values = [100 * v for k, v in r.items() if k.startswith('loy_content_')]
        right.plot([i, i], [min(values), max(values)], color=INK, lw=1.4)
        right.plot([i] * len(values), values, marker='o', markersize=4, color=INK, linestyle='none')
    for y in (0, 50, 100):
        right.axhline(y, color='#a9b9af', lw=1 if y == 50 else .7)
    right.set_xticks(range(len(genotype)), [f"{r['well']}\n{r['analysis_genotype']}" for r in genotype], fontsize=8)
    right.set_ylabel('Estimated LOY content (%)', color=INK)
    right.set_title('Same-genotype wells: classifier estimate', loc='left', fontsize=11, color=INK, fontweight='bold',
                    pad=20)
    for ax in (left, right):
        ax.grid(axis='y', color=GRID, linewidth=.7)
        ax.set_axisbelow(True)
        ax.tick_params(colors=MUTED)
        ax.spines[['left', 'bottom']].set_color('#a9b9af')
    fig.text(.07, .045, 'Right: four estimates per well (two training arms × two reporter classes), each corrected '
                        'with held-out TPR/FPR; labels under the wells are the analysis genotypes.', fontsize=8.8, color=MUTED)
    save(fig, figures / 'identity_checks', paths)


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(dict.fromkeys(k for r in rows for k in r)))
        writer.writeheader()
        writer.writerows(rows)


def render(report: Path) -> None:
    report = report.resolve()
    results_path = report / 'results.json'
    results_hash = digest(results_path)
    current = json.loads(results_path.read_text(encoding='utf-8'))
    previous_path = REPO / current['previous_report']
    if digest(previous_path) != current['previous_sha256']:
        raise ValueError('The previous report changed after this report was built')
    previous = json.loads(previous_path.read_text(encoding='utf-8'))
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'axes.titlesize': 13, 'axes.labelsize': 10,
                         'svg.fonttype': 'none', 'pdf.fonttype': 42, 'axes.spines.top': False,
                         'axes.spines.right': False})
    figures = report / 'figures'
    figures.mkdir(exist_ok=True)
    paths: list[Path] = []
    rows, control_rows = series(current, previous), control_series(current, previous)
    trajectories(current, previous, rows, figures, paths)
    controls(current, previous, control_rows, figures, paths)
    pairs = technical(current, figures, paths)
    identity(current, figures, paths)
    for name, data in (('trajectories_data.csv', rows), ('same_genotype_data.csv', control_rows),
                       ('technical_replicates_data.csv', pairs)):
        write_csv(figures / name, data)
        paths.append(figures / name)
    if digest(results_path) != results_hash:
        raise ValueError('Analysis changed during plotting')
    manifest = dict(generated_utc=datetime.now(timezone.utc).isoformat(), results_sha256=results_hash,
                    previous_sha256=current['previous_sha256'], matplotlib_version=matplotlib.__version__,
                    script_sha256=digest(Path(__file__)), files={str(p.relative_to(report)): digest(p) for p in paths})
    (figures / 'plots_manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    print(json.dumps(dict(figures=[str(p) for p in paths], measurements_unchanged=True), indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--report-dir', type=Path, required=True)
    render(parser.parse_args().report_dir)
