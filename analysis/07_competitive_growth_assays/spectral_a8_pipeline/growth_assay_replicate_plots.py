"""Export individual D4 trajectories against their shared D0 references."""
from __future__ import annotations

import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator

REPLICATES = {1: ('#2563ab', 'o', '-'), 2: ('#ad5b14', 's', '--'), 3: ('#7c3ca0', '^', ':')}
ORDER = ('WT', 'WT_C1', 'same_genotype', 'DDX3X', 'DDX3X_C1', 'SREBF1', 'TADA2B')


def panels(data: dict) -> list[dict]:
    """Keep approved replicate identities and distinguish LOY from GFP readouts."""
    result = []
    for group in ORDER:
        controls = group == 'same_genotype'
        for orientation in (('CTRL', 'LOY') if controls else ('LOY_GFP', 'LOY_mCherry')):
            if controls:
                reference = next(s for s in data['dye_controls'] if s['genotype'] == orientation)
                baseline = reference['d0_gfp_fraction_pooled']
                source_wells = reference['d0_wells']
                wells = [w for w in data['wells'] if w['role'] == 'dye_neutral' and w['genotype'] == orientation]
                title = f'{orientation} / {orientation}'
                subtitle = 'Same genotype · GFP + mCherry · WT background'
            else:
                wells = [w for w in data['wells'] if w['role'] == 'mixture' and w['group'] == group
                         and w['orientation'] == orientation]
                baseline = wells[0]['d0_loy_fraction']
                source_wells = wells[0]['d0_wells']
                title = group.replace('_', ' ') + ' · LOY = ' + ('GFP' if orientation == 'LOY_GFP' else 'mCherry')
                subtitle = 'LOY / CTRL competition'
            wells.sort(key=lambda w: w['replicate'])
            if [w['replicate'] for w in wells] != [1, 2, 3]:
                raise ValueError(f'Expected exactly three approved replicates: {title}')
            rows = []
            for w in wells:
                fraction = w['gfp_fraction'] if controls else w['loy_fraction']
                delta = (math.log(fraction / (1 - fraction)) - math.log(baseline / (1 - baseline))
                         if fraction is not None and 0 < fraction < 1 and 0 < baseline < 1 else None)
                if not controls and delta is not None and not math.isclose(delta, w['delta_log_odds'], abs_tol=1e-12):
                    raise ValueError('Replicate comparison differs from reviewed analysis')
                rows.append(dict(sample=w['analysis_sample'], original_filename=w['original_filename'],
                                 well=w['well'], replicate=w['replicate'], fraction=fraction, delta_log_odds=delta))
            result.append(dict(group=group, orientation=orientation, title=title, subtitle=subtitle,
                               readout='GFP' if controls else 'LOY', baseline=baseline,
                               d0_wells=source_wells, rows=rows))
    return result


def spread_labels(values: list[float], gap: float, lo: float, hi: float) -> list[float]:
    order = sorted(range(len(values)), key=lambda i: values[i])
    labels = [max(lo, min(hi, values[i])) for i in order]
    for i in range(1, len(labels)):
        labels[i] = max(labels[i], labels[i - 1] + gap)
    if labels and labels[-1] > hi:
        labels[-1] = hi
        for i in range(len(labels) - 2, -1, -1):
            labels[i] = min(labels[i], labels[i + 1] - gap)
    result = values.copy()
    for i, value in zip(order, labels):
        result[i] = value
    return result


def render_replicates(data: dict, figures: Path) -> list[Path]:
    groups = panels(data)
    outputs = []
    deltas = [r['delta_log_odds'] for p in groups for r in p['rows'] if r['delta_log_odds'] is not None]
    limits = (min(-.5, math.floor(min(deltas)*2)/2), max(.5, math.ceil(max(deltas)*2)/2))
    for normalized in (False, True):
        # Paired columns keep reciprocal dye orientations side by side in plate-row order.
        fig, axes = plt.subplots(7, 2, figsize=(11, 20), sharex=True, sharey=True)
        fig.subplots_adjust(left=.09, right=.96, top=.91, bottom=.055, hspace=.55, wspace=.27)
        fig.suptitle('Each replicate · D0 to D4', x=.09, y=.982, ha='left',
                     fontsize=23, fontweight='bold', color='#172d24')
        fig.text(.09, .959, ('Change in log-odds from D0' if normalized else 'Fraction among retained cells (%)')
                 + ' · approved replicate labels · broad-negative exclusion', fontsize=11, color='#52635b')
        legend = [Line2D([], [], color=c, marker=m, linestyle=s, lw=1.6, markersize=5, label=f'Rep {r}')
                  for r, (c, m, s) in REPLICATES.items()]
        legend.append(Line2D([], [], color='#52635b', marker='D', linestyle='none', markersize=5,
                             label='Shared pooled D0'))
        fig.legend(handles=legend, loc='upper left', bbox_to_anchor=(.08, .946), ncol=4, frameon=False)
        lo, hi = limits if normalized else (0, 100)
        for ax, panel in zip(axes.flat, groups):
            ax.set_title(panel['title'], loc='left', fontsize=13, color='#172d24', fontweight='bold', pad=21)
            ax.text(0, 1.04, panel['subtitle'], transform=ax.transAxes, fontsize=9, color='#52635b')
            ax.set_xlim(-.4, 5.05)
            ax.set_ylim(lo, hi)
            ax.set_xticks([0, 4], ['D0', 'D4'])
            ax.tick_params(axis='x', labelbottom=True)
            ax.tick_params(axis='y', labelleft=True)
            ax.yaxis.set_major_locator(MultipleLocator(.5 if normalized else 20))
            ax.set_ylabel(('Δ ' + panel['readout'] + ' log-odds') if normalized else panel['readout'] + ' (%)')
            ax.set_axisbelow(True)
            ax.grid(axis='y', color='#dfe7e1', lw=.7)
            ax.spines[['left', 'bottom']].set_color('#a9b9af')
            ax.tick_params(colors='#52635b', length=3)
            if normalized:
                ax.axhline(0, color='#84968b', lw=.8, linestyle=':')
            start = 0 if normalized else 100 * panel['baseline']
            valid = [r for r in panel['rows'] if r['delta_log_odds' if normalized else 'fraction'] is not None]
            ends = [r['delta_log_odds'] if normalized else 100*r['fraction'] for r in valid]
            label_y = spread_labels(ends, (hi-lo)*.085, lo+(hi-lo)*.04, hi-(hi-lo)*.04)
            for row, end, label in zip(valid, ends, label_y):
                color, marker, style = REPLICATES[row['replicate']]
                ax.plot([0, 4], [start, end], color=color, lw=1.5, linestyle=style, alpha=.9)
                ax.plot(4, end, color=color, marker=marker, markersize=5, markeredgecolor='white',
                        markeredgewidth=.5, linestyle='none', zorder=4)
                ax.plot([4.08, 4.35], [end, label], color=color, lw=.6)
                ax.text(4.4, label, f'R{row["replicate"]} · {row["well"]}', va='center', fontsize=8, color=color)
            ax.plot(0, start, color='#52635b', marker='D', markersize=5, linestyle='none', zorder=5)
        fig.text(.09, .029, 'Each line follows an approved D4 replicate from its shared pooled D0 reference. D0 samples are not paired to individual D4 wells.',
                 fontsize=8.5, color='#52635b')
        fig.text(.09, .017, 'Competition panels track LOY; same-genotype panels track GFP. Only D0 and D4 were measured. DDX3X arms share their D0 mix.',
                 fontsize=8.5, color='#52635b')
        stem = 'replicate_trajectories_log_odds' if normalized else 'replicate_trajectories_fraction'
        for ext in ('png', 'svg', 'pdf'):
            path = figures / f'{stem}.{ext}'
            fig.savefig(path, dpi=180, facecolor='white')
            outputs.append(path)
        plt.close(fig)
    csv_path = figures / 'replicate_trajectories_data.csv'
    rows = [dict(group=p['group'], orientation=p['orientation'], readout=p['readout'],
                 d0_fraction=p['baseline'], d0_source_wells=p['d0_wells'],
                 baseline_type='shared_pooled_reference', **r) for p in groups for r in p['rows']]
    with csv_path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    outputs.append(csv_path)
    return outputs
