"""Placeholder D0-to-D4 knockout competition figure from the approved D4 report tables.

Reads the saved comp_assay outputs (row-C genotype swap applied), writes tidy source
tables, and renders a-d. No FCS processing, refitting or new statistical test.
"""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import json
import os
import sys

THESIS = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(THESIS.parent / '.figure_python_deps'))
sys.path.insert(0, str(THESIS / 'skills/thesis-figure-style/scripts'))
OUT = THESIS / 'figures/chapter4'
DATA = OUT / 'source_data/competition_d0_d4_draft'
REPORT = THESIS.parent / 'comp_assay/analysis/spectral_growth_assay_D4_20261020'
os.environ.setdefault('MPLCONFIGDIR', str(DATA / 'mpl_cache'))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from thesis_style import apply_style

# Reporter colours follow the Chapter 4 competition QC figure; genotype colours follow the thesis mapping.
GREEN, CHERRY, GREY, INK = '#32965B', '#D34F4F', '#777E86', '#222222'
COLOURS = json.loads((THESIS / 'skills/thesis-figure-style/assets/condition_colours.json').read_text())['palette']
CTRL_BLUE, LOY_RED = COLOURS['Y-retaining'], COLOURS['Y-depleted']
ARMS = ('WT', 'WT_C1', 'DDX3X', 'DDX3X_C1', 'SREBF1', 'TADA2B')
NAMES = {'WT': 'WT', 'WT_C1': 'WT C1', 'DDX3X': 'DDX3X', 'DDX3X_C1': 'DDX3X C1', 'SREBF1': 'SREBF1', 'TADA2B': 'TADA2B'}
ORIENT = {'LOY_GFP': (GREEN, 'o', '-', 'LOY–GFP'), 'LOY_mCherry': (CHERRY, 's', (0, (3, 1.6)), 'LOY–mCherry')}
WIDTH, HEIGHT = 183, 101
INPUTS = ['wells.csv', 'orientation_summaries.csv', 'dye_agreement.csv', 'dye_controls.csv', 'd0_baselines.csv',
          'identity_mapping.json', 'results.json', 'same_genotype/wells.csv']


def digest(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def signed(v):
    return f'{v:+.2f}'.replace('-', '−')


def logit(p):
    return np.log(p / (1 - p))


def load():
    hashes = {name: digest(REPORT / name) for name in INPUTS}
    results = json.loads((REPORT / 'results.json').read_text(encoding='utf-8'))
    if results['identity_status'] != 'user_confirmed' or results['qc'].get('corrected_genotype_count') != 6:
        raise ValueError('Expected the approved D4 labels with the row-C genotype swap')
    wells = pd.read_csv(REPORT / 'wells.csv')
    mix = wells.loc[wells.role.eq('mixture'), ['well', 'original_filename', 'group', 'orientation', 'replicate',
                                               'd0_loy_fraction', 'd0_wells', 'loy_fraction', 'delta_log_odds',
                                               'excluded_broad_fraction']].copy()
    summary = pd.read_csv(REPORT / 'orientation_summaries.csv')
    agree = pd.read_csv(REPORT / 'dye_agreement.csv').set_index('group')
    # Balanced effect: mean of the two orientation means; contrasts are differences of these saved means.
    agree['balanced'] = (agree.gfp_delta_log_odds + agree.mcherry_delta_log_odds) / 2
    contrasts = []
    for arm, ref, family in [('DDX3X', 'WT', 'knockout vs WT'), ('SREBF1', 'WT', 'knockout vs WT'),
                             ('TADA2B', 'WT', 'knockout vs WT'), ('WT_C1', 'WT', 'C1 arm vs untreated'),
                             ('DDX3X_C1', 'DDX3X', 'C1 arm vs untreated')]:
        contrasts.append(dict(arm=arm, reference=ref, family=family,
                              **{k: agree.loc[arm, c] - agree.loc[ref, c] for k, c in
                                 [('LOY_GFP', 'gfp_delta_log_odds'), ('LOY_mCherry', 'mcherry_delta_log_odds'),
                                  ('balanced', 'balanced')]}))
    contrasts = pd.DataFrame(contrasts)
    same = pd.read_csv(REPORT / 'same_genotype/wells.csv')
    controls = pd.read_csv(REPORT / 'dye_controls.csv').set_index('genotype')
    for frame, name in [(mix, 'd4_competition_wells.csv'), (summary, 'orientation_summaries.csv'),
                        (agree.reset_index(), 'arm_effects.csv'), (contrasts, 'arm_contrasts.csv'),
                        (same, 'same_genotype_wells.csv'), (controls.reset_index(), 'same_genotype_summary.csv')]:
        frame.to_csv(DATA / name, index=False)
    return mix, summary, agree, contrasts, same, controls, hashes


def mm_axes(fig, x, y, w, h):
    return fig.add_axes([x / WIDTH, y / HEIGHT, w / WIDTH, h / HEIGHT])


def style_axis(ax):
    ax.tick_params(labelsize=6.4, length=1.6, pad=1.3, width=.55)
    ax.xaxis.labelpad = 2
    ax.yaxis.labelpad = 2


def render(mix, summary, agree, contrasts, same, controls):
    apply_style(7.5)
    plt.rcParams['axes.unicode_minus'] = True
    fig = plt.figure(figsize=(WIDTH / 25.4, HEIGHT / 25.4))
    label = lambda x, y, s, **kw: fig.text(x / WIDTH, y / HEIGHT, s, **kw)
    jitter = np.array([-.09, 0, .09])

    # a: LOY share in each arm, pooled D0 mixture to the three D4 wells.
    label(1.5, 98.6, 'a', fontsize=9, fontweight='semibold', va='center')
    for i, arm in enumerate(ARMS):
        ax = mm_axes(fig, 13 + 28.6 * i, 64, 23.5, 28)
        for orientation, (colour, marker, style, _) in ORIENT.items():
            d = mix[(mix.group == arm) & (mix.orientation == orientation)].sort_values('replicate')
            d0 = 100 * d.d0_loy_fraction.iloc[0]
            for x4, y4 in zip(1 + jitter, 100 * d.loy_fraction):
                ax.plot([0, x4], [d0, y4], color=colour, lw=.45, alpha=.55, ls=style, zorder=1)
            ax.plot([0], [d0], marker=marker, color=colour, ms=3.4, mfc=colour, mec='white', mew=.4, zorder=3)
            ax.plot(1 + jitter, 100 * d.loy_fraction, marker, color=colour, ms=2.9, mfc='white', mew=.65, zorder=3)
        ax.set(xlim=(-.35, 1.35), ylim=(20, 72), xticks=[0, 1], xticklabels=['D0', 'D4'], yticks=[20, 40, 60])
        style_axis(ax)
        if i:
            ax.set_yticklabels([])
        else:
            ax.set_ylabel('LOY of retained cells (%)', fontsize=6.8)
        ax.text(.5, 1.03, NAMES[arm], transform=ax.transAxes, ha='center', va='bottom', fontsize=7, fontweight='semibold')
        if arm in ('DDX3X', 'DDX3X_C1'):
            ax.text(.5, .955, 'shared D0 mix', transform=ax.transAxes, ha='center', va='top', fontsize=5.8, color=GREY)
        if arm == 'TADA2B':
            ax.text(.5, .955, 'unlabelled-cell gate', transform=ax.transAxes, ha='center', va='top', fontsize=5.8, color=GREY)
    handles = [plt.Line2D([], [], color=c, marker=m, ls=s, lw=.9, ms=3, mfc='white', mew=.65, label=n)
               for c, m, s, n in ORIENT.values()]
    handles.append(plt.Line2D([], [], color=INK, marker='D', ls='none', ms=3.2, label='Mean of orientations (b, c)'))
    fig.legend(handles=handles, loc='center left', bbox_to_anchor=(92 / WIDTH, 98.5 / HEIGHT), ncol=3,
               fontsize=6.5, handlelength=2.2, columnspacing=1.2, handletextpad=.45, frameon=False)
    label(13, 98.5, 'Pooled D0 mixture (filled) to three D4 wells (open)', fontsize=6.4, color=GREY, va='center')

    # b: change in log-odds per well, orientation means and the dye-balanced mean.
    label(1.5, 54, 'b', fontsize=9, fontweight='semibold', va='center')
    ax = mm_axes(fig, 13, 8.5, 63, 42)
    ax.axhline(0, color='#A7ADB5', lw=.6, ls=(0, (3, 3)), zorder=0)
    for i, arm in enumerate(ARMS):
        for offset, orientation in [(-.2, 'LOY_GFP'), (.2, 'LOY_mCherry')]:
            colour, marker, _, _ = ORIENT[orientation]
            d = mix[(mix.group == arm) & (mix.orientation == orientation)]
            ax.plot(i + offset + jitter * .45, d.delta_log_odds, marker, color=colour, ms=2.8, mfc='white', mew=.6, zorder=3)
            m = d.delta_log_odds.mean()
            ax.plot([i + offset - .1, i + offset + .1], [m, m], color=colour, lw=1.1, solid_capstyle='butt', zorder=2)
        b = agree.loc[arm, 'balanced']
        ax.plot([i], [b], 'D', color=INK, ms=3.4, zorder=4)
        top = max(mix[(mix.group == arm)].delta_log_odds.max(), b)
        ax.text(i, top + .07, signed(b), ha='center', va='bottom', fontsize=6.3)
    ax.set(xlim=(-.6, 5.6), ylim=(-.45, 1.35), xticks=range(6), xticklabels=[NAMES[a] for a in ARMS],
           ylabel='Δ ln(LOY:CTRL), D4 − D0')
    style_axis(ax)
    ax.tick_params(axis='x', labelsize=6.3)

    # c: differences between arms in the same saved effects.
    label(79, 54, 'c', fontsize=9, fontweight='semibold', va='center')
    ax = mm_axes(fig, 89.5, 8.5, 45, 42)
    ax.axhline(0, color='#A7ADB5', lw=.6, ls=(0, (3, 3)), zorder=0)
    xs = [0, 1.15, 2.3, 3.9, 5.05]
    for x, (_, r) in zip(xs, contrasts.iterrows()):
        for offset, orientation in [(-.18, 'LOY_GFP'), (.18, 'LOY_mCherry')]:
            colour, marker, _, _ = ORIENT[orientation]
            ax.plot([x + offset], [r[orientation]], marker, color=colour, ms=3, mfc='white', mew=.65, zorder=3)
        ax.plot([x], [r.balanced], 'D', color=INK, ms=3.4, zorder=4)
        ax.text(x, .2, signed(r.balanced), ha='center', va='center', fontsize=6.3)
    ax.set(xlim=(-.5, 5.55), ylim=(-1.4, .55), xticks=xs, yticks=[-1.25, -1, -.75, -.5, -.25, 0, .25],
           xticklabels=['DDX3X', 'SREBF1', 'TADA2B', 'WT C1', 'DDX3X\nC1'], ylabel='Difference in Δ ln(LOY:CTRL)')
    style_axis(ax)
    ax.tick_params(axis='x', labelsize=6.3)
    ax.axvline(3.1, color='#D5D8DC', lw=.5, zorder=0)
    ax.text(1.15, .5, 'Knockout − WT', ha='center', va='top', fontsize=6.3, color=GREY)
    ax.text(4.475, .5, 'C1 − untreated', ha='center', va='top', fontsize=6.3, color=GREY)

    # d: same-genotype reporter-balance controls (corrected row-C labels).
    label(139, 54, 'd', fontsize=9, fontweight='semibold', va='center')
    ax = mm_axes(fig, 149.5, 8.5, 23, 42)
    for genotype, colour, marker, xoff in [('CTRL', CTRL_BLUE, 'o', -.08), ('LOY', LOY_RED, 's', .08)]:
        d0 = 100 * same[(same.genotype == genotype) & (same.day == 0)].gfp_fraction.to_numpy()
        d4 = 100 * same[(same.genotype == genotype) & (same.day == 4)].gfp_fraction.to_numpy()
        p0, p4 = 100 * controls.loc[genotype, 'd0_gfp_fraction_pooled'], 100 * controls.loc[genotype, 'd4_gfp_fraction_mean']
        ax.plot([0 + xoff, 1 + xoff], [p0, p4], color=colour, lw=1, zorder=2)
        ax.plot(0 + xoff + jitter * .4, d0, marker, color=colour, ms=2.7, mfc='white', mew=.6, zorder=3)
        ax.plot(1 + xoff + jitter * .4, d4, marker, color=colour, ms=2.7, mfc='white', mew=.6, zorder=3)
        ax.text(1.42, p4, f"{genotype}\n{signed(controls.loc[genotype, 'mean_delta_log_odds'])}", color=colour,
                ha='left', va='center', fontsize=6.2, linespacing=1.05, clip_on=False)
    ax.axhline(50, color='#A7ADB5', lw=.6, ls=(0, (3, 3)), zorder=0)
    ax.set(xlim=(-.35, 1.35), ylim=(40, 65), xticks=[0, 1], xticklabels=['D0', 'D4'], yticks=[40, 50, 60],
           ylabel='GFP of retained cells (%)')
    style_axis(ax)
    label(161, 52.6, 'Same-genotype mixtures', fontsize=6.4, color=GREY, ha='center', va='center')

    stem = OUT / 'ch04_main_competition_d0_d4_draft'
    fig.savefig(stem.with_suffix('.pdf'))
    fig.savefig(stem.with_suffix('.svg'))
    fig.savefig(stem.with_suffix('.png'), dpi=600)
    plt.close(fig)
    return stem


def main():
    DATA.mkdir(parents=True, exist_ok=True)
    mix, summary, agree, contrasts, same, controls, hashes = load()
    stem = render(mix, summary, agree, contrasts, same, controls)
    if any(digest(REPORT / name) != h for name, h in hashes.items()):
        raise ValueError('Source report changed while the figure was rendered')
    manifest = dict(generated_utc=datetime.now(timezone.utc).isoformat(), report=str(REPORT), input_sha256=hashes,
                    script_sha256=digest(__file__), matplotlib=matplotlib.__version__, figure_mm=[WIDTH, HEIGHT],
                    exports={p.name: digest(p) for p in (stem.with_suffix(s) for s in ('.pdf', '.svg', '.png'))},
                    note='Placeholder: D0 and D4 only. Values are copied from the approved D4 report; '
                         'contrasts are differences of saved orientation means. No inferential test.')
    (DATA / 'manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    print(stem.with_suffix('.png'))


if __name__ == '__main__':
    main()
