"""Main draft: do the CRISPR screen hits replicate in the competition assay?

Time-course layout (a: per-arm trajectories with a WT reference, b: knockout - WT,
c: same-genotype controls). Each follow-up report is listed in TIMEPOINTS; today only
D4 exists, and later days are added as one entry each. Values are read from the approved
reports (adopted gating: TADA2B unlabelled-cell gate, standard gate elsewhere); no
reprocessing or inferential test.
"""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import json
import os
import sys

THESIS = Path(__file__).resolve().parents[3]
COMP = THESIS.parent / 'comp_assay'
sys.path.insert(0, str(THESIS.parent / '.figure_python_deps'))
sys.path.insert(0, str(THESIS / 'skills/thesis-figure-style/scripts'))
OUT = THESIS / 'figures/chapter4'
DATA = OUT / 'source_data/competition_hits_draft'
os.environ.setdefault('MPLCONFIGDIR', str(DATA / 'mpl_cache'))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from thesis_style import apply_style

# (day, approved follow-up report directory). Add later timepoints here.
TIMEPOINTS = [(4, COMP / 'analysis/spectral_growth_assay_D4_20261020')]
ARMS = ('WT', 'DDX3X', 'SREBF1', 'TADA2B')
HITS = ('DDX3X', 'SREBF1', 'TADA2B')
ROLE = {'WT': 'no knockout', 'DDX3X': 'positive control', 'SREBF1': 'screen hit', 'TADA2B': 'screen hit'}
GREEN, CHERRY, GREY, INK, WT_GHOST = '#32965B', '#D34F4F', '#777E86', '#222222', '#B9BEC4'
COLOURS = json.loads((THESIS / 'skills/thesis-figure-style/assets/condition_colours.json').read_text())['palette']
CTRL_BLUE, LOY_RED = COLOURS['Y-retaining'], COLOURS['Y-depleted']
ORIENT = {'LOY_GFP': (GREEN, 'o', '-', 'LOY–GFP'), 'LOY_mCherry': (CHERRY, 's', (0, (3, 1.6)), 'LOY–mCherry')}
WIDTH, HEIGHT = 183, 84


def digest(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def signed(v):
    v = round(float(v), 2)
    return '0.00' if v == 0 else f'{v:+.2f}'.replace('-', '−')


def load():
    wells, controls, hashes, gating = [], [], {}, {}
    for day, report in TIMEPOINTS:
        results = json.loads((report / 'results.json').read_text(encoding='utf-8'))
        if results['identity_status'] != 'user_confirmed':
            raise ValueError(f'{report.name}: labels not approved')
        gating[day] = results.get('gating', {}).get('unlabelled_gate_primary_backgrounds', [])
        for name in ('results.json', 'wells.csv', 'same_genotype/wells.csv', 'dye_controls.csv'):
            hashes[f'{report.name}/{name}'] = digest(report / name)
        w = pd.read_csv(report / 'wells.csv')
        w = w[(w.role == 'mixture') & w.group.isin(ARMS)]
        wells.append(w.assign(day=day)[['day', 'well', 'group', 'orientation', 'replicate', 'd0_loy_fraction',
                                        'loy_fraction', 'delta_log_odds', 'primary_gate']])
        s = pd.read_csv(report / 'same_genotype/wells.csv')
        p0 = s[s.day == 0].groupby('genotype').apply(lambda g: g.gfp_count.sum() / g.reporter_positive_singlets.sum(),
                                                     include_groups=False)
        d = s[s.day == 4].copy()
        d['delta_log_odds'] = np.log(d.gfp_fraction / (1 - d.gfp_fraction)) - d.genotype.map(np.log(p0 / (1 - p0)))
        controls.append(d.assign(day=day)[['day', 'well', 'genotype', 'replicate', 'gfp_fraction', 'delta_log_odds']])
    wells, controls = pd.concat(wells), pd.concat(controls)
    days = [0] + [d for d, _ in TIMEPOINTS]
    # Orientation means per day, the mean of orientations, and knockout - WT at each day.
    means = wells.groupby(['day', 'group', 'orientation']).delta_log_odds.agg(['mean', 'min', 'max']).reset_index()
    arm = means.groupby(['day', 'group'])['mean'].mean().rename('mean_of_orientations').reset_index()
    contrast = []
    for day in days[1:]:
        for hit in HITS:
            for orientation in ORIENT:
                ko = wells[(wells.day == day) & (wells.group == hit) & (wells.orientation == orientation)].delta_log_odds
                wt = means[(means.day == day) & (means.group == 'WT') & (means.orientation == orientation)]['mean'].iloc[0]
                contrast += [dict(day=day, hit=hit, orientation=orientation, well_value=v - wt) for v in ko]
    contrast = pd.DataFrame(contrast)
    DATA.mkdir(parents=True, exist_ok=True)
    wells.to_csv(DATA / 'competition_wells.csv', index=False)
    means.to_csv(DATA / 'orientation_means.csv', index=False)
    arm.to_csv(DATA / 'arm_means.csv', index=False)
    contrast.to_csv(DATA / 'knockout_minus_wt.csv', index=False)
    controls.to_csv(DATA / 'same_genotype_wells.csv', index=False)
    return dict(wells=wells, means=means, arm=arm, contrast=contrast, controls=controls, days=days,
                hashes=hashes, gating=gating)


def render(d):
    apply_style(7.5)
    plt.rcParams['axes.unicode_minus'] = True
    fig = plt.figure(figsize=(WIDTH / 25.4, HEIGHT / 25.4))
    label = lambda x, y, s, **kw: fig.text(x / WIDTH, y / HEIGHT, s, **kw)
    axes = lambda x, y, w, h: fig.add_axes([x / WIDTH, y / HEIGHT, w / WIDTH, h / HEIGHT])
    days, last = d['days'], d['days'][-1]
    xpad = max(1.2, .08 * last)

    def tidy(ax):
        ax.tick_params(labelsize=6.3, length=1.6, pad=1.2, width=.55)
        ax.xaxis.labelpad = 1.5
        ax.yaxis.labelpad = 1.5

    def trajectory(arm):
        series = {}
        for orientation in ORIENT:
            m = d['means'][(d['means'].group == arm) & (d['means'].orientation == orientation)].sort_values('day')
            series[orientation] = (days, [0] + m['mean'].tolist())
        a = d['arm'][d['arm'].group == arm].sort_values('day')
        series['mean'] = (days, [0] + a.mean_of_orientations.tolist())
        return series

    # a: per-arm trajectories from each arm's own D0, with the WT mean as a grey reference.
    label(1.5, 81.5, 'a', fontsize=9, fontweight='semibold', va='center')
    wt = trajectory('WT')['mean']
    lo = min(-.45, d['wells'].delta_log_odds.min() - .1)
    hi = max(1.25, d['wells'].delta_log_odds.max() + .1)
    for i, arm in enumerate(ARMS):
        ax = axes(14 + 42.5 * i, 46, 37, 28)
        ax.axhline(0, color='#A7ADB5', lw=.6, ls=(0, (3, 3)), zorder=0)
        if arm != 'WT':
            ax.plot(*wt, color=WT_GHOST, lw=2.2, solid_capstyle='round', zorder=1)
        series = trajectory(arm)
        for orientation, (colour, marker, style, _) in ORIENT.items():
            ax.plot(*series[orientation], color=colour, lw=.9, ls=style, zorder=2)
            w = d['wells'][(d['wells'].group == arm) & (d['wells'].orientation == orientation)]
            offset = (-1 if orientation == 'LOY_GFP' else 1) * .025 * last
            ax.plot(w.day + offset, w.delta_log_odds, marker, color=colour, ms=2.6, mfc='white', mew=.6, zorder=3)
        ax.plot(*series['mean'], color=INK, lw=1.3, zorder=4)
        ax.plot(series['mean'][0][1:], series['mean'][1][1:], 'D', color=INK, ms=3, zorder=5)
        ax.text(last + xpad * .55, series['mean'][1][-1], signed(series['mean'][1][-1]), fontsize=6.2, va='center',
                ha='left', clip_on=False)
        ax.set(xlim=(-xpad, last + xpad), ylim=(lo, hi), xticks=days, xticklabels=[f'D{x}' for x in days])
        tidy(ax)
        if i:
            ax.set_yticklabels([])
        else:
            ax.set_ylabel('Δ ln(LOY:CTRL) from D0', fontsize=6.8)
        ax.text(.03, 1.04, arm, transform=ax.transAxes, fontsize=7.2, fontweight='semibold', va='bottom')
        note = ROLE[arm] + (' · unlabelled-cell gate' if arm in d['gating'][last] else '')
        ax.text(.03, .97, note, transform=ax.transAxes, fontsize=5.9, color=GREY, va='top', ha='left')
    handles = [plt.Line2D([], [], color=c, marker=m, ls=s, lw=.9, ms=3, mfc='white', mew=.6, label=n)
               for c, m, s, n in ORIENT.values()]
    handles += [plt.Line2D([], [], color=INK, marker='D', lw=1.3, ms=3, label='Mean of orientations'),
                plt.Line2D([], [], color=WT_GHOST, lw=2.2, label='WT mean')]
    fig.legend(handles=handles, loc='center right', bbox_to_anchor=(181 / WIDTH, 81.5 / HEIGHT), ncol=4, fontsize=6.3,
               handlelength=2, columnspacing=1.2, handletextpad=.4, frameon=False)

    # b: knockout - WT for each hit at the latest day.
    label(1.5, 37, 'b', fontsize=9, fontweight='semibold', va='center')
    ax = axes(14, 8, 88, 25)
    ax.axhline(0, color='#A7ADB5', lw=.6, ls=(0, (3, 3)), zorder=0)
    c = d['contrast'][d['contrast'].day == last]
    jitter = np.array([-.035, 0, .035])
    for i, hit in enumerate(HITS):
        means = []
        for offset, orientation in ((-.16, 'LOY_GFP'), (.16, 'LOY_mCherry')):
            colour, marker, _, _ = ORIENT[orientation]
            v = c[(c.hit == hit) & (c.orientation == orientation)].well_value.to_numpy()
            ax.plot(i + offset + jitter[:len(v)], v, marker, color=colour, ms=2.8, mfc='white', mew=.6)
            ax.plot([i + offset - .07, i + offset + .07], [v.mean()] * 2, color=colour, lw=1.1)
            means.append(v.mean())
        ax.plot([i], [np.mean(means)], 'D', color=INK, ms=3.4)
        ax.text(i, .32, signed(np.mean(means)), ha='center', fontsize=6.3)
    ax.set(xlim=(-.5, 2.5), ylim=(-1.4, .45), xticks=range(len(HITS)),
           xticklabels=[f'{h} − WT' for h in HITS], ylabel=f'Difference at D{last}')
    tidy(ax)
    ax.tick_params(axis='x', labelsize=6.4)

    # c: same-genotype controls over time.
    label(112, 37, 'c', fontsize=9, fontweight='semibold', va='center')
    ax = axes(124, 8, 40, 25)
    ax.axhline(0, color='#A7ADB5', lw=.6, ls=(0, (3, 3)), zorder=0)
    for genotype, colour, marker, offset in (('CTRL', CTRL_BLUE, 'o', -.025), ('LOY', LOY_RED, 's', .025)):
        w = d['controls'][d['controls'].genotype == genotype]
        m = [0] + [w[w.day == x].delta_log_odds.mean() for x in days[1:]]
        ax.plot(days, m, color=colour, lw=1.1)
        ax.plot(w.day + offset * last, w.delta_log_odds, marker, color=colour, ms=2.6, mfc='white', mew=.6)
        ax.text(last + xpad * .55, m[-1], f'{genotype}/{genotype}\n{signed(m[-1])}', color=colour, fontsize=5.9,
                va='center', ha='left', linespacing=1.05, clip_on=False)
    ax.set(xlim=(-xpad, last + xpad), ylim=(-.3, .7), xticks=days, xticklabels=[f'D{x}' for x in days],
           ylabel='Δ ln(GFP:mCherry)')
    tidy(ax)
    ax.text(.03, 1.04, 'Same-genotype mixtures', transform=ax.transAxes, fontsize=6.4, color=GREY, va='bottom')

    stem = OUT / 'ch04_main_competition_hits_draft'
    fig.savefig(stem.with_suffix('.pdf'))
    fig.savefig(stem.with_suffix('.svg'))
    fig.savefig(stem.with_suffix('.png'), dpi=600)
    plt.close(fig)
    return stem


def main():
    d = load()
    stem = render(d)
    for (day, report) in TIMEPOINTS:
        for name in ('results.json', 'wells.csv', 'same_genotype/wells.csv', 'dye_controls.csv'):
            if digest(report / name) != d['hashes'][f'{report.name}/{name}']:
                raise ValueError('Source report changed while rendering')
    manifest = dict(generated_utc=datetime.now(timezone.utc).isoformat(), timepoints=[[day, str(r)] for day, r in TIMEPOINTS],
                    input_sha256=d['hashes'], gating={f'D{k}': v for k, v in d['gating'].items()},
                    script_sha256=digest(__file__), matplotlib=matplotlib.__version__, figure_mm=[WIDTH, HEIGHT],
                    exports={p.name: digest(p) for p in (stem.with_suffix(s) for s in ('.pdf', '.svg', '.png'))},
                    note='Values copied from approved reports; knockout - WT uses each knockout well minus the WT '
                         'orientation mean at the same day. No inferential test.')
    (DATA / 'manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    print(stem.with_suffix('.png'))


if __name__ == '__main__':
    main()
