"""Main and supplementary drafts: knockout competition assay from D0 to D14.

Main: per-arm trajectories with both replicate sets (a), D14 effects against the range of the
two WT sets (b), same-genotype mixtures (c), and the WT orientation gap against the sum of the
same-genotype changes (d). Supplement: technical duplicates (a), reporter-negative fractions (b),
bounds for reporter-negative events (c) and identity checks (d, e).

Values are read from the D14 report; its D4 points are checked against the approved D4 report.
The only new arithmetic is means and sums of saved values. No reprocessing or inferential test.
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
DATA = OUT / 'source_data/competition_d14_draft'
REPORT = COMP / 'analysis/spectral_growth_assay_D14_20261029'
D4_REPORT = COMP / 'analysis/spectral_growth_assay_D4_20261020'
os.environ.setdefault('MPLCONFIGDIR', str(DATA / 'mpl_cache'))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from thesis_style import apply_style

# Reporter colours follow the Chapter 4 competition figures; genotype colours follow the thesis mapping.
GREEN, CHERRY, GREY, INK, BAND, ZERO = '#32965B', '#D34F4F', '#777E86', '#222222', '#E3E6EA', '#A7ADB5'
PALE = {'LOY_GFP': '#A6D6B7', 'LOY_mCherry': '#EFAAAA'}
COLOURS = json.loads((THESIS / 'skills/thesis-figure-style/assets/condition_colours.json').read_text())['palette']
CTRL_BLUE, LOY_RED = COLOURS['Y-retaining'], COLOURS['Y-depleted']
ORIENT = {'LOY_GFP': (GREEN, 'o', '-', 'LOY–GFP'), 'LOY_mCherry': (CHERRY, 's', (0, (3, 1.6)), 'LOY–mCherry')}
# The second sets carry the C1 label but received no drug (experimenter, 29 October 2026).
ARMS = [('WT', 'WT', 'set 1'), ('WT_C1', 'WT', 'set 2'), ('DDX3X', 'DDX3X', 'set 1'),
        ('DDX3X_C1', 'DDX3X', 'set 2'), ('SREBF1', 'SREBF1', ''), ('TADA2B', 'TADA2B', '')]
ROLE = {'WT': 'no knockout', 'DDX3X': 'positive control', 'SREBF1': 'screen hit', 'TADA2B': 'screen hit'}
DAYS = (0, 4, 14)
MAIN_MM, SUPP_MM = (183, 112), (183, 100)
INPUTS = {REPORT: ['results.json', 'wells.csv', 'biological_replicates.csv', 'figures/trajectories_data.csv',
                   'figures/same_genotype_data.csv', 'figures/technical_replicates_data.csv',
                   'identity_size_fingerprint.csv', 'identity_same_genotype.csv'],
          D4_REPORT: ['results.json', 'wells.csv']}


def digest(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def signed(v):
    v = round(float(v), 2)
    return '0.00' if v == 0 else f'{v:+.2f}'.replace('-', '−')


def load():
    hashes = {f'{r.name}/{n}': digest(r / n) for r, names in INPUTS.items() for n in names}
    results = json.loads((REPORT / 'results.json').read_text(encoding='utf-8'))
    if json.loads((D4_REPORT / 'results.json').read_text(encoding='utf-8'))['identity_status'] != 'user_confirmed':
        raise ValueError('D4 labels are not approved')
    if results['identity_status'] not in ('user_confirmed', 'recorded_labels_checked'):
        raise ValueError(f"Unexpected D14 identity status {results['identity_status']}")
    if Path(results['previous_report']).parent.name != D4_REPORT.name:
        raise ValueError('D14 report was not built on the approved D4 report')
    traj = pd.read_csv(REPORT / 'figures/trajectories_data.csv')
    bio = pd.read_csv(REPORT / 'biological_replicates.csv')
    wells = pd.read_csv(REPORT / 'wells.csv')
    d4 = pd.read_csv(D4_REPORT / 'wells.csv')
    d4 = d4[d4.role == 'mixture']

    # The D14 report's D4 points must equal the approved D4 report, and its D14 points the replicate table.
    check4 = traj[traj.kind == 'well'].merge(d4, left_on=['group', 'orientation', 'wells'],
                                             right_on=['group', 'orientation', 'well'], suffixes=('', '_d4'))
    reps = traj[traj.kind == 'biological_replicate'].astype({'replicate': int})
    check14 = reps.merge(bio, on=['group', 'orientation', 'replicate'], suffixes=('', '_bio'))
    if not (len(check4) == 36 and np.allclose(check4.delta_log_odds, check4.delta_log_odds_d4, rtol=0, atol=1e-9)
            and len(check14) == 36 and np.allclose(check14.delta_log_odds, check14.delta_log_odds_bio, rtol=0, atol=1e-9)):
        raise ValueError('Trajectory points do not match the D4 report or the D14 replicate table')

    # Mean of orientations per arm and day (D0 is zero by definition).
    means = traj[traj.kind.isin(['pooled_d0', 'arm_mean'])][['group', 'orientation', 'day', 'delta_log_odds']]
    arm = means.groupby(['group', 'day']).delta_log_odds.mean().rename('mean_of_orientations').reset_index()

    # Biological replicates whose wells carry the relative reporter-negative flag (>20% and >2x the D0 mix).
    flagged = wells[wells.qc_flags.fillna('').str.contains('reporter-negative')].well
    bio['reporter_negative_flag'] = bio.wells.str.split(';').apply(lambda ws: bool(set(ws) & set(flagged)))

    same = pd.read_csv(REPORT / 'figures/same_genotype_data.csv')
    sg = same[same.kind.isin(['pooled_d0', 'arm_mean'])].groupby(['genotype', 'day']).gfp_delta_log_odds.first().unstack(0)
    gaps = []
    for day in DAYS:
        row = dict(day=day, predicted_from_same_genotype=sg.loc[day, 'LOY'] + sg.loc[day, 'CTRL'],
                   loy_loy=sg.loc[day, 'LOY'], ctrl_ctrl=sg.loc[day, 'CTRL'])
        for group, name in (('WT', 'wt_set_1'), ('WT_C1', 'wt_set_2')):
            m = means[(means.group == group) & (means.day == day)].set_index('orientation').delta_log_odds
            row[name] = m['LOY_GFP'] - m['LOY_mCherry']
        gaps.append(row)
    gaps = pd.DataFrame(gaps)

    excluded = pd.concat([
        wells[wells.role == 'mixture'].groupby(['group', 'orientation']).d0_excluded_fraction.first()
        .rename('excluded_fraction').reset_index().assign(day=0, unit='pooled D0 mix'),
        d4[['group', 'orientation', 'well', 'excluded_broad_fraction']].rename(
            columns={'excluded_broad_fraction': 'excluded_fraction'}).assign(day=4, unit='well'),
        bio[bio.group != 'same_genotype'][['group', 'orientation', 'replicate', 'excluded_broad_fraction']].rename(
            columns={'excluded_broad_fraction': 'excluded_fraction'}).assign(day=14, unit='culture (mean of 2 wells)')])

    size = pd.read_csv(REPORT / 'identity_size_fingerprint.csv')
    size = size[size.role == 'mixture']
    genotype = pd.read_csv(REPORT / 'identity_same_genotype.csv')
    estimates = [c for c in genotype.columns if c.startswith('loy_content_') and '_to_' in c]
    genotype['loy_content_mean'] = genotype[estimates].mean(axis=1)
    tech = pd.read_csv(REPORT / 'figures/technical_replicates_data.csv')

    DATA.mkdir(parents=True, exist_ok=True)
    for frame, name in [(traj, 'trajectory_points.csv'), (arm, 'mean_of_orientations.csv'), (bio, 'd14_cultures.csv'),
                        (same, 'same_genotype_points.csv'), (gaps, 'wt_orientation_gap.csv'),
                        (excluded, 'reporter_negative_fractions.csv'), (tech, 'technical_duplicates.csv'),
                        (size, 'size_fingerprint.csv'), (genotype, 'same_genotype_identity.csv')]:
        frame.to_csv(DATA / name, index=False)
    return dict(traj=traj, means=means, arm=arm, bio=bio, same=same, gaps=gaps, excluded=excluded, size=size,
                genotype=genotype, tech=tech, hashes=hashes, status=results['identity_status'],
                gating=results['gating']['unlabelled_gate_primary_backgrounds'], elapsed=results['elapsed_days_from_d0'])


def canvas(width, height):
    apply_style(7.5)
    plt.rcParams['axes.unicode_minus'] = True
    fig = plt.figure(figsize=(width / 25.4, height / 25.4))
    label = lambda x, y, s, **kw: fig.text(x / width, y / height, s, **kw)
    axes = lambda x, y, w, h: fig.add_axes([x / width, y / height, w / width, h / height])
    return fig, label, axes


def letter(label, x, y, s):
    label(x, y, s, fontsize=9, fontweight='semibold', va='center')


def tidy(ax):
    ax.tick_params(labelsize=6.3, length=1.6, pad=1.2, width=.55)
    ax.xaxis.labelpad = 1.5
    ax.yaxis.labelpad = 1.5


def zero(ax):
    ax.axhline(0, color=ZERO, lw=.6, ls=(0, (3, 3)), zorder=0)


def day_axis(ax):
    ax.set(xlim=(-1.6, 15.6), xticks=DAYS, xticklabels=[f'D{x}' for x in DAYS])


def render_main(d):
    width, height = MAIN_MM
    fig, label, axes = canvas(width, height)
    traj, means, arm, bio = d['traj'], d['means'], d['arm'], d['bio']
    series = lambda group: [arm[(arm.group == group) & (arm.day == x)].mean_of_orientations.iloc[0] for x in DAYS]
    orient = lambda group, o: [means[(means.group == group) & (means.orientation == o) & (means.day == x)]
                               .delta_log_odds.iloc[0] for x in DAYS]
    wt = np.array([series('WT'), series('WT_C1')])
    wt_lo, wt_hi = wt.min(axis=0), wt.max(axis=0)
    negative = {(r.group, r.orientation, r.replicate) for r in bio.itertuples() if r.reporter_negative_flag}

    # a: trajectories per arm and set, with the range of the two WT sets as a grey band.
    letter(label, 1.5, 108.5, 'a')
    points = traj[traj.kind.isin(['well', 'biological_replicate'])].astype({'replicate': int})
    lo, hi = min(-.75, points.delta_log_odds.min() - .12), max(1.2, points.delta_log_odds.max() + .12)
    x0, w, gap = 13, 25.6, 2.9
    for i, (group, background, set_label) in enumerate(ARMS):
        ax = axes(x0 + (w + gap) * i, 61, w, 34)
        zero(ax)
        ax.fill_between(DAYS, wt_lo, wt_hi, color=BAND, lw=0, zorder=.5)
        for orientation, (colour, marker, style, _) in ORIENT.items():
            ax.plot(DAYS, orient(group, orientation), color=colour, lw=.9, ls=style, zorder=2)
            offset = -.5 if orientation == 'LOY_GFP' else .5
            for p in points[(points.group == group) & (points.orientation == orientation)].itertuples():
                fill = PALE[orientation] if p.day == 14 and (group, orientation, p.replicate) in negative else 'white'
                ax.plot(p.day + offset, p.delta_log_odds, marker, color=colour, ms=2.6, mfc=fill, mew=.6, zorder=3)
        ax.plot(DAYS, series(group), color=INK, lw=1.3, zorder=4)
        ax.plot(DAYS[1:], series(group)[1:], 'D', color=INK, ms=3, zorder=5)
        ax.set(ylim=(lo, hi))
        day_axis(ax)
        tidy(ax)
        title = background + (f' · {set_label}' if set_label else '')
        ax.text(0, 1.035, title, transform=ax.transAxes, fontsize=6.9, fontweight='semibold', va='bottom', ha='left')
        ax.text(1, 1.035, signed(series(group)[-1]), transform=ax.transAxes, fontsize=6.3, va='bottom', ha='right')
        note = ROLE[background] + ('\nunlabelled-cell gate' if group in d['gating'] else '')
        ax.text(.04, .03, note, transform=ax.transAxes, fontsize=5.7, color=GREY, va='bottom', ha='left', linespacing=1.1)
        if i:
            ax.set_yticklabels([])
        else:
            ax.set_ylabel('Δ ln(LOY:CTRL) from D0', fontsize=6.8)
    # Brackets mark how the paired sets were mixed.
    for first, text in ((0, 'separate D0 mixes'), (2, 'one D0 mix per orientation, split')):
        left, right = x0 + (w + gap) * first, x0 + (w + gap) * (first + 1) + w
        fig.add_artist(plt.Line2D([left / width, left / width, right / width, right / width],
                                  [99.2 / height, 100.2 / height, 100.2 / height, 99.2 / height], color=GREY, lw=.5))
        label((left + right) / 2, 100.8, text, fontsize=5.9, color=GREY, ha='center', va='bottom')
    handles = [plt.Line2D([], [], color=c, marker=m, ls=s, lw=.9, ms=3, mfc='white', mew=.6, label=n)
               for c, m, s, n in ORIENT.values()]
    handles += [plt.Line2D([], [], color=INK, marker='D', lw=1.3, ms=3, label='Mean of orientations'),
                Patch(facecolor=BAND, edgecolor='none', label='Range of the two WT sets'),
                plt.Line2D([], [], color=GREY, marker='o', ls='none', ms=3, mfc='#DADDE1', mew=.6,
                           label='Shaded: >20% reporter-negative')]
    fig.legend(handles=handles, loc='center right', bbox_to_anchor=(181 / width, 108.5 / height), ncol=5, fontsize=6.2,
               handlelength=2, columnspacing=1.1, handletextpad=.4, frameon=False)

    # b: D14 effect per culture, orientation means and mean of orientations, against the WT-set range.
    letter(label, 1.5, 51, 'b')
    ax = axes(13, 10.5, 74, 34.5)
    zero(ax)
    ax.axhspan(wt_lo[-1], wt_hi[-1], color=BAND, lw=0, zorder=.5)
    ax.text(-.5, wt_lo[-1] - .025, 'WT sets', fontsize=5.7, color=GREY, va='top', ha='left')
    jitter = np.array([-.05, 0, .05])
    for i, (group, background, set_label) in enumerate(ARMS):
        orientation_means = []
        for offset, orientation in ((-.19, 'LOY_GFP'), (.19, 'LOY_mCherry')):
            colour, marker, _, _ = ORIENT[orientation]
            b = bio[(bio.group == group) & (bio.orientation == orientation)].sort_values('replicate')
            for j, r in enumerate(b.itertuples()):
                fill = PALE[orientation] if r.reporter_negative_flag else 'white'
                ax.plot(i + offset + jitter[j], r.delta_log_odds, marker, color=colour, ms=2.8, mfc=fill, mew=.6, zorder=3)
            m = b.delta_log_odds.mean()
            ax.plot([i + offset - .09, i + offset + .09], [m, m], color=colour, lw=1.1, solid_capstyle='butt', zorder=2)
            orientation_means.append(m)
        mean = np.mean(orientation_means)
        ax.plot([i], [mean], 'D', color=INK, ms=3.4, zorder=4)
        top = max(bio[bio.group == group].delta_log_odds.max(), mean)
        ax.text(i, top + .06, signed(mean), ha='center', va='bottom', fontsize=6.2)
    blo = min(-.8, bio.delta_log_odds.min() - .1)
    ax.set(xlim=(-.55, 5.55), ylim=(blo, .75), xticks=range(len(ARMS)),
           xticklabels=[f'{b}\n{s}' if s else b for _, b, s in ARMS], ylabel='Δ ln(LOY:CTRL), D0 to D14')
    tidy(ax)
    ax.tick_params(axis='x', labelsize=6.2)

    # c: same-genotype WT mixtures over time.
    letter(label, 93, 51, 'c')
    ax = axes(103, 10.5, 28, 34.5)
    zero(ax)
    same = d['same']
    for genotype, colour, marker, offset in (('CTRL', CTRL_BLUE, 'o', -.5), ('LOY', LOY_RED, 's', .5)):
        g = same[same.genotype == genotype]
        m = [0] + [g[(g.day == x) & (g.kind == 'arm_mean')].gfp_delta_log_odds.iloc[0] for x in DAYS[1:]]
        ax.plot(DAYS, m, color=colour, lw=1.1, zorder=2)
        p = g[g.kind.isin(['well', 'biological_replicate'])]
        ax.plot(p.day + offset, p.gfp_delta_log_odds, marker, color=colour, ms=2.6, mfc='white', mew=.6, zorder=3)
        ax.text(16.2, m[-1], f'{genotype}/{genotype}\n{signed(m[-1])}', color=colour, fontsize=5.9, va='center',
                ha='left', linespacing=1.05, clip_on=False)
    ax.set(ylim=(-.65, .55), ylabel='Δ ln(GFP:mCherry) from D0')
    day_axis(ax)
    tidy(ax)
    ax.text(0, 1.035, 'Same-genotype WT mixtures', transform=ax.transAxes, fontsize=6.3, color=GREY, va='bottom')

    # d: WT orientation gap against the sum of the two same-genotype changes.
    letter(label, 145, 51, 'd')
    ax = axes(155.5, 10.5, 25.5, 34.5)
    zero(ax)
    g = d['gaps']
    ax.plot(g.day, g.predicted_from_same_genotype, color=INK, lw=1.2, marker='D', ms=3, zorder=4,
            label='Sum of c')
    ax.plot(g.day, g.wt_set_1, color=GREY, lw=.9, marker='o', ms=3, mfc=GREY, mec='white', mew=.4, zorder=3,
            label='WT set 1')
    ax.plot(g.day, g.wt_set_2, color=GREY, lw=.9, ls=(0, (2, 1.4)), marker='^', ms=3.2, mfc='white', mew=.6, zorder=3,
            label='WT set 2')
    ax.set(ylim=(-.65, .7), ylabel='LOY–GFP minus LOY–mCherry')
    day_axis(ax)
    tidy(ax)
    ax.legend(loc='upper right', fontsize=5.5, handlelength=1.8, handletextpad=.35, borderaxespad=.3, frameon=False,
              labelspacing=.25)
    ax.text(0, 1.035, 'WT orientation gap', transform=ax.transAxes, fontsize=6.3, color=GREY, va='bottom')

    return save(fig, 'ch04_main_competition_d14_draft')


def render_supp(d):
    width, height = SUPP_MM
    fig, label, axes = canvas(width, height)
    bio = d['bio']
    mix = bio[bio.group != 'same_genotype']
    slots = lambda ax: (ax.set(xlim=(-.55, 5.55), xticks=range(len(ARMS))),
                        ax.set_xticklabels([f'{b}\n{s}' if s else b for _, b, s in ARMS], fontsize=6.2))

    # a: the two technical wells of every culture, with the identity line and a 5-point band.
    letter(label, 1.5, 97, 'a')
    ax = axes(13, 57, 33, 33)
    t = d['tech']
    ax.fill_between([0, 100], [-5, 95], [5, 105], color=BAND, lw=0, zorder=0)
    ax.plot([0, 100], [0, 100], color=ZERO, lw=.6, ls=(0, (3, 3)), zorder=1)
    style = {'LOY_GFP': (GREEN, 'o'), 'LOY_mCherry': (CHERRY, 's'), 'CTRL': (CTRL_BLUE, 'o'), 'LOY': (LOY_RED, 's')}
    for r in t.itertuples():
        colour, marker = style[r.orientation]
        ax.plot(100 * r.technical_1, 100 * r.technical_2, marker, color=colour, ms=2.6, mfc='white', mew=.6, zorder=3)
        if r.gap_pp > 5:
            ax.text(100 * r.technical_1 - 1.2, 100 * r.technical_2 + 1.2, r.wells.replace(';', '/'), fontsize=5.5,
                    color=GREY, va='bottom', ha='right')
    ax.set(xlim=(15, 65), ylim=(15, 65), xticks=[20, 40, 60], yticks=[20, 40, 60], aspect='equal',
           xlabel='First well (%)', ylabel='Second well (%)')
    tidy(ax)
    big = int((t.gap_pp > 5).sum())
    ax.text(.96, .04, f'median gap {np.median(t.gap_pp):.1f} pp\n{big} of {len(t)} > 5 pp', transform=ax.transAxes,
            fontsize=5.7, color=GREY, va='bottom', ha='right', linespacing=1.1)
    ax.text(0, 1.035, 'Technical duplicates, D14', transform=ax.transAxes, fontsize=6.3, color=GREY, va='bottom')

    # b: reporter-negative events by day: D0 mix (bar), D4 wells (open), D14 cultures (filled).
    letter(label, 54, 97, 'b')
    ax = axes(64, 57, 117, 33)
    ex = d['excluded']
    jitter = np.array([-.03, 0, .03])
    for i, (group, _, _) in enumerate(ARMS):
        for offset, orientation in ((-.21, 'LOY_GFP'), (.21, 'LOY_mCherry')):
            colour, marker, _, _ = ORIENT[orientation]
            e = ex[(ex.group == group) & (ex.orientation == orientation)]
            x = i + offset
            d0 = 100 * e[e.day == 0].excluded_fraction.iloc[0]
            ax.plot([x - .15, x + .15], [d0, d0], color=INK, lw=1, solid_capstyle='butt', zorder=2)
            ax.plot(x - .07 + jitter, 100 * e[e.day == 4].excluded_fraction, marker, color=colour, ms=2.5,
                    mfc='white', mew=.6, zorder=3)
            ax.plot(x + .07 + jitter, 100 * e[e.day == 14].excluded_fraction, marker, color=colour, ms=2.7,
                    mfc=colour, mec='white', mew=.35, zorder=3)
    slots(ax)
    ax.set(ylim=(0, 75), ylabel='Reporter-negative (% of singlets)')
    tidy(ax)
    handles = [plt.Line2D([], [], color=INK, lw=1, label='D0 mix'),
               plt.Line2D([], [], color=GREY, marker='o', ls='none', ms=2.8, mfc='white', mew=.6, label='D4 wells'),
               plt.Line2D([], [], color=GREY, marker='o', ls='none', ms=2.9, mfc=GREY, mec='white', mew=.35,
                          label='D14 cultures')]
    ax.legend(handles=handles, loc='upper left', fontsize=5.8, frameon=False, ncol=3, handlelength=1.6,
              columnspacing=1, handletextpad=.4, borderaxespad=.3)
    ax.text(5, 100 * ex[(ex.group == 'TADA2B') & (ex.day == 14)].excluded_fraction.max() + 3,
            'unlabelled-cell gate', fontsize=5.6, color=GREY, ha='center', va='bottom')

    # c: D14 change with every reporter-negative event counted as CTRL or as LOY (line ends).
    letter(label, 1.5, 47, 'c')
    ax = axes(13, 10.5, 84, 32)
    zero(ax)
    for i, (group, _, _) in enumerate(ARMS):
        for offset, orientation in ((-.19, 'LOY_GFP'), (.19, 'LOY_mCherry')):
            colour, marker, _, _ = ORIENT[orientation]
            b = mix[(mix.group == group) & (mix.orientation == orientation)].sort_values('replicate')
            for j, r in enumerate(b.itertuples()):
                x = i + offset + 1.6 * jitter[j]
                ax.plot([x, x], [r.delta_log_odds_excluded_as_ctrl, r.delta_log_odds_excluded_as_loy], color=colour,
                        lw=.8, alpha=.75, solid_capstyle='butt', zorder=2)
                ax.plot(x, r.delta_log_odds, marker, color=colour, ms=2.6, mfc='white', mew=.6, zorder=3)
    slots(ax)
    ax.set(ylim=(-1.4, 1.8), ylabel='Δ ln(LOY:CTRL), D0 to D14')
    tidy(ax)

    # d: WT size fingerprint: LOY cells are larger, so the GFP:mCherry size ratio follows the orientation.
    letter(label, 101, 47, 'd')
    ax = axes(112, 10.5, 25, 32)
    s = d['size']
    ax.axhline(1, color=ZERO, lw=.6, ls=(0, (3, 3)), zorder=0)
    for x, orientation in enumerate(ORIENT):
        colour, marker, _, _ = ORIENT[orientation]
        for k, (group, fill) in enumerate((('WT', 'white'), ('WT_C1', colour))):
            v = s[(s.group == group) & (s.orientation == orientation)].gfp_to_mcherry_ratio.to_numpy()
            xs = x + (-.14 if k == 0 else .14) + np.linspace(-.07, .07, len(v))
            ax.plot(xs, v, marker, color=colour, ms=2.6, mfc=fill, mec=colour if fill == 'white' else 'white',
                    mew=.6 if fill == 'white' else .35, zorder=3)
    ax.set_yscale('log', base=2)
    ax.set(xlim=(-.5, 1.5), ylim=(.7, 1.65), xticks=[0, 1], yticks=[.75, 1, 1.25, 1.5])
    ax.set_xticklabels(['LOY–GFP', 'LOY–mCherry'], fontsize=6.1)
    ax.set_yticklabels(['0.75', '1.00', '1.25', '1.50'])
    ax.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.set_ylabel('Size ratio, GFP:mCherry cells', fontsize=6.5)
    tidy(ax)
    ok = int(s.consistent.sum())
    ax.text(.5, .53, f'{ok}/{len(s)} as expected\nopen set 1, filled set 2', transform=ax.transAxes, fontsize=5.5,
            color=GREY, va='center', ha='center', linespacing=1.1)
    ax.text(0, 1.035, 'WT size check', transform=ax.transAxes, fontsize=6.3, color=GREY, va='bottom')

    # e: estimated LOY content of the same-genotype wells (mean of four estimates; lines span them).
    letter(label, 143, 47, 'e')
    ax = axes(154, 10.5, 27, 32)
    g = d['genotype'].sort_values('well')
    for ref in (0, 100):
        ax.axhline(ref, color=ZERO, lw=.6, ls=(0, (3, 3)), zorder=0)
    for k, r in enumerate(g.itertuples()):
        x = k + (.8 if r.recorded_genotype == 'CTRL' else 0)
        colour, marker = (LOY_RED, 's') if r.recorded_genotype == 'LOY' else (CTRL_BLUE, 'o')
        ax.plot([x, x], [100 * r.loy_content_minimum, 100 * r.loy_content_maximum], color=colour, lw=.8, zorder=2)
        ax.plot(x, 100 * r.loy_content_mean, marker, color=colour, ms=2.6, mfc='white', mew=.6, zorder=3)
    ax.set(xlim=(-.8, 12.6), ylim=(-30, 130), xticks=[2.5, 9.3], yticks=[0, 50, 100])
    ax.set_xticklabels(['LOY/LOY', 'CTRL/CTRL'], fontsize=6.1)
    ax.set_ylabel('Estimated LOY content (%)', fontsize=6.5)
    tidy(ax)
    ax.text(0, 1.035, 'Same-genotype wells', transform=ax.transAxes, fontsize=6.3, color=GREY, va='bottom')

    return save(fig, 'ch04_supp_competition_d14_checks_draft')


def save(fig, name):
    stem = OUT / name
    fig.savefig(stem.with_suffix('.pdf'))
    fig.savefig(stem.with_suffix('.svg'))
    fig.savefig(stem.with_suffix('.png'), dpi=600)
    plt.close(fig)
    return stem


def main():
    d = load()
    stems = [render_main(d), render_supp(d)]
    roots = {r.name: r for r in INPUTS}
    if any(digest(roots[k.split('/', 1)[0]] / k.split('/', 1)[1]) != h for k, h in d['hashes'].items()):
        raise ValueError('Source report changed while rendering')
    manifest = dict(generated_utc=datetime.now(timezone.utc).isoformat(), report=str(REPORT), d4_report=str(D4_REPORT),
                    d14_identity_status=d['status'], elapsed_days_from_d0_read=d['elapsed'],
                    unlabelled_gate_primary_backgrounds=d['gating'], input_sha256=d['hashes'],
                    script_sha256=digest(__file__), matplotlib=matplotlib.__version__,
                    figure_mm={'main': MAIN_MM, 'supp': SUPP_MM},
                    exports={p.name: digest(p) for s in stems for p in (s.with_suffix(x) for x in ('.pdf', '.svg', '.png'))},
                    note='Draft. D14 labels are recorded labels checked by the size and genotype classifiers, not yet '
                         'confirmed by the experimenter. Values are copied from the reports; means of orientations and '
                         'the predicted orientation gap are the only arithmetic. No inferential test.')
    (DATA / 'manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    for s in stems:
        print(s.with_suffix('.png'))


if __name__ == '__main__':
    main()
