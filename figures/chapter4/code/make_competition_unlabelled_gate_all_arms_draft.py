"""D0-to-D4 change for every arm and same-genotype control under three reporter-negative gates.

Primary: eGFP < 6 and mCherry < 6.5 on both days. Background-specific unlabelled-cell gate: the
saved sensitivity (valley cuts only where an unlabelled cluster is detected on both days, i.e.
TADA2B). Uniform: the TADA2B valley cuts applied to every background, shown for comparison.
Reads hash-verified FCS files through the comp_assay A8 pipeline; no classifier fitting or test.
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
sys.path.insert(0, str(COMP))
OUT = THESIS / 'figures/chapter4'
DATA = OUT / 'source_data/competition_unlabelled_gate_all_arms_draft'
REPORT = COMP / 'analysis/spectral_growth_assay_D4_20261020'
os.environ.setdefault('MPLCONFIGDIR', str(DATA / 'mpl_cache'))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from thesis_style import apply_style
from spectral.a8_pipeline.competition_report import load_cef_gates
from spectral.a8_pipeline.fcs_io import load_fcs
from spectral.a8_pipeline.gating import points_in_polygon
from spectral.a8_pipeline.growth_assay_followup import GATE_NAMES
from spectral.a8_pipeline.growth_assay_unlabelled_gate import gated_counts

GREEN, CHERRY, GREY, INK = '#32965B', '#D34F4F', '#777E86', '#222222'
COLOURS = json.loads((THESIS / 'skills/thesis-figure-style/assets/condition_colours.json').read_text())['palette']
CTRL_BLUE, LOY_RED = COLOURS['Y-retaining'], COLOURS['Y-depleted']
ARMS = ('WT', 'WT_C1', 'DDX3X', 'DDX3X_C1', 'SREBF1', 'TADA2B')
NAMES = {'WT': 'WT', 'WT_C1': 'WT C1', 'DDX3X': 'DDX3X', 'DDX3X_C1': 'DDX3X C1', 'SREBF1': 'SREBF1', 'TADA2B': 'TADA2B'}
ORIENT = {'LOY_GFP': (GREEN, 'o', 'LOY–GFP'), 'LOY_mCherry': (CHERRY, 's', 'LOY–mCherry')}
RULES = ('primary', 'background', 'uniform')
RULE_STYLE = {'primary': dict(mfc='white', label='Standard gate (mCherry < 6.5)'),
              'background': dict(mfc=None, label='Adopted: unlabelled-cell gate for TADA2B only'),
              'uniform': dict(mfc='half', label='TADA2B cuts applied to every arm')}
WIDTH, HEIGHT = 183, 92


def digest(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def signed(v):
    v = round(float(v), 2)
    return '0.00' if v == 0 else f'{v:+.2f}'.replace('-', '−')


def logit(p):
    return float(np.log(p / (1 - p)))


def compute():
    d4 = json.loads((REPORT / 'results.json').read_text(encoding='utf-8'))
    d0 = json.loads((COMP / d4['baseline_path']).read_text(encoding='utf-8'))
    sens_path = REPORT / 'unlabelled_gate_sensitivity.json'
    sens = json.loads(sens_path.read_text(encoding='utf-8'))
    if sens['identity_mapping_sha256'] != d4['identity_mapping_sha256']:
        raise ValueError('Sensitivity file was built from different labels')
    rule = sens['rule']
    primary_cut = rule['primary_mcherry_below']
    tada2b = rule['mcherry_below']['TADA2B']
    cut = {'primary': lambda bg, day: primary_cut,
           'background': lambda bg, day: rule['mcherry_below'][bg][f'D{day}'],
           'uniform': lambda bg, day: tada2b[f'D{day}']}
    weights, intercept = d0['classifier']['weights'], d0['classifier']['intercept']
    gfp_below = d0['negative_rules']['broad']['gfp_below']
    counts, meta, hashes = {}, {}, {}
    for day, data in ((0, d0), (4, d4)):
        expected = {m['path']: m['sha256'] for m in data['source_manifest']}
        source = COMP / data['source']
        gates = load_cef_gates(next(source.glob('*.cef')))
        raw = {p.stem: p for p in source.rglob('*.fcs')}
        for w in data['wells']:
            if w['role'] not in ('mixture', 'dye_neutral'):
                continue
            path = raw[w['sample'] if day == 0 else w['original_sample']]
            rel = path.relative_to(COMP).as_posix()
            if digest(path) != expected[rel]:
                raise ValueError(f'Raw file hash changed: {rel}')
            hashes[rel] = expected[rel]
            events = load_fcs(path).events
            keep = np.ones(len(events), bool)
            for name in GATE_NAMES:
                gate = gates[name]
                keep &= points_in_polygon(events[list(gate.parameters)].to_numpy(float), gate.vertices)
            f = np.arcsinh(events.loc[keep, ['eGFP-A', 'mCherry-A']].to_numpy(float) / 1000)
            meta[(day, w['well'])] = w
            for r in RULES:
                counts[(r, day, w['well'])] = gated_counts(f, weights, intercept, gfp_below, cut[r](w['background'], day))
    rows = []
    for r in RULES:
        base = {}
        for b in d0['d0_baselines']:
            wells = b['wells'].split(';')
            kept = sum(counts[(r, 0, x)]['retained'] for x in wells)
            green = sum(counts[(r, 0, x)]['gfp_count'] for x in wells)
            base[(b['group'], b['orientation'])] = (green if b['orientation'] == 'LOY_GFP' else kept - green) / kept
        neutral0 = {}
        for g in ('CTRL', 'LOY'):
            wells = [w['well'] for (day, _), w in meta.items() if day == 0 and w['role'] == 'dye_neutral' and w['genotype'] == g]
            neutral0[g] = sum(counts[(r, 0, x)]['gfp_count'] for x in wells) / sum(counts[(r, 0, x)]['retained'] for x in wells)
        for (day, well), w in meta.items():
            if day != 4:
                continue
            c = counts[(r, 4, well)]
            gfp = c['gfp_count'] / c['retained']
            if w['role'] == 'mixture':
                loy = gfp if w['orientation'] == 'LOY_GFP' else 1 - gfp
                p0 = base[(w['group'], w['orientation'])]
                rows.append(dict(rule=r, kind='competition', group=w['group'], orientation=w['orientation'], well=well,
                                 replicate=w['replicate'], d0_fraction=p0, d4_fraction=loy, delta=logit(loy) - logit(p0)))
            else:
                p0 = neutral0[w['genotype']]
                rows.append(dict(rule=r, kind='same_genotype', group=w['genotype'], orientation='GFP fraction', well=well,
                                 replicate=w['replicate'], d0_fraction=p0, d4_fraction=gfp, delta=logit(gfp) - logit(p0)))
    df = pd.DataFrame(rows)
    # Reproduce the saved primary report and background-specific sensitivity before plotting.
    primary = pd.read_csv(REPORT / 'wells.csv').set_index('well')
    sens_wells = {w['well']: w for w in sens['wells']}
    for row in df[df.kind == 'competition'].itertuples():
        expected = primary.loc[row.well, 'delta_log_odds_standard_gate'] if row.rule == 'primary' else sens_wells[row.well]['delta_log_odds'] if row.rule == 'background' else None
        if expected is not None and abs(expected - row.delta) > 1e-9:
            raise ValueError(f'{row.rule} result not reproduced for {row.well}')
    DATA.mkdir(parents=True, exist_ok=True)
    df.to_csv(DATA / 'well_deltas_by_rule.csv', index=False)
    summary = (df.groupby(['rule', 'kind', 'group', 'orientation']).delta.agg(['mean', 'min', 'max', 'count']).reset_index())
    summary.to_csv(DATA / 'summary_by_rule.csv', index=False)
    inputs = {str(p): digest(p) for p in (REPORT / 'results.json', sens_path, REPORT / 'wells.csv', COMP / d4['baseline_path'])}
    return df, summary, dict(primary=primary_cut, tada2b=tada2b, background=rule['mcherry_below']), hashes, inputs


def marker_kwargs(rule, colour, marker, ms):
    style = RULE_STYLE[rule]
    kw = dict(marker=marker, color=colour, ms=ms, mew=.65, ls='none')
    if style['mfc'] == 'half':
        kw.update(fillstyle='bottom', mfc=colour, mfcalt='white')
    else:
        kw.update(mfc=style['mfc'] or colour)
    return kw


def render(df, summary, cuts):
    apply_style(7.5)
    plt.rcParams['axes.unicode_minus'] = True
    fig = plt.figure(figsize=(WIDTH / 25.4, HEIGHT / 25.4))
    label = lambda x, y, s, **kw: fig.text(x / WIDTH, y / HEIGHT, s, **kw)
    axes = lambda x, y, w, h: fig.add_axes([x / WIDTH, y / HEIGHT, w / WIDTH, h / HEIGHT])

    def tidy(ax):
        ax.tick_params(labelsize=6.3, length=1.6, pad=1.2, width=.55)
        ax.yaxis.labelpad = 1.5
    offsets = dict(zip(RULES, (-.075, 0, .075)))
    mean = lambda rule, kind, group, orientation: summary[(summary.rule == rule) & (summary.kind == kind) & (summary.group == group)
                                                         & (summary.orientation == orientation)].iloc[0]

    # a: every arm, both orientations, three gates (mean and well range).
    label(1.5, 89, 'a', fontsize=9, fontweight='semibold', va='center')
    ax = axes(14, 51, 167, 33)
    ax.axhline(0, color='#A7ADB5', lw=.6, ls=(0, (3, 3)), zorder=0)
    for i, arm in enumerate(ARMS):
        if i:
            ax.axvline(i - .5, color='#E3E5E8', lw=.5, zorder=0)
        for j, (orientation, (colour, marker, _)) in enumerate(ORIENT.items()):
            x0 = i + (-.2 if j == 0 else .2)
            for rule in RULES:
                s = mean(rule, 'competition', arm, orientation)
                x = x0 + offsets[rule]
                ax.plot([x, x], [s['min'], s['max']], color=colour, lw=.6, alpha=.7)
                ax.plot([x], [s['mean']], **marker_kwargs(rule, colour, marker, 3.1))
    ax.set(xlim=(-.55, 5.55), ylim=(-.45, 1.35), xticks=range(6), xticklabels=[NAMES[a] for a in ARMS],
           ylabel='Δ ln(LOY:CTRL), D4 − D0')
    tidy(ax)
    handles = [plt.Line2D([], [], **marker_kwargs(r, GREY, 'o', 3.4), label=RULE_STYLE[r]['label']) for r in RULES]
    fig.legend(handles=handles, loc='center left', bbox_to_anchor=(12 / WIDTH, 89 / HEIGHT), ncol=3, fontsize=6.3,
               handletextpad=.3, columnspacing=1.4, frameon=False)
    label(181, 85.6, 'Markers: orientation mean; lines: range of three D4 wells. Green: LOY–GFP; red: LOY–mCherry.',
          fontsize=5.9, color=GREY, ha='right', va='center')

    # b: mean of orientations per arm under each gate.
    label(1.5, 42, 'b', fontsize=9, fontweight='semibold', va='center')
    ax = axes(14, 9, 105, 30)
    ax.axhline(0, color='#A7ADB5', lw=.6, ls=(0, (3, 3)), zorder=0)
    for i, arm in enumerate(ARMS):
        for rule in RULES:
            m = np.mean([mean(rule, 'competition', arm, o)['mean'] for o in ORIENT])
            ax.plot([i + offsets[rule] * 2.2], [m], **marker_kwargs(rule, INK, 'D', 3.4))
        b = np.mean([mean('background', 'competition', arm, o)['mean'] for o in ORIENT])
        p = np.mean([mean('primary', 'competition', arm, o)['mean'] for o in ORIENT])
        u = np.mean([mean('uniform', 'competition', arm, o)['mean'] for o in ORIENT])
        ax.text(i, max(b, p, u) + .1, signed(b), ha='center', fontsize=6.2)
    ax.set(xlim=(-.55, 5.55), ylim=(-.3, 1.0), xticks=range(6), xticklabels=[NAMES[a] for a in ARMS],
           ylabel='Mean of orientations')
    tidy(ax)
    ax.tick_params(axis='x', labelsize=6.2)
    ax.text(.99, .97, 'value: adopted gating', transform=ax.transAxes, ha='right', va='top', fontsize=5.8, color=GREY)

    # c: same-genotype controls (GFP balance) under each gate.
    label(124, 42, 'c', fontsize=9, fontweight='semibold', va='center')
    ax = axes(135, 9, 34, 30)
    ax.axhline(0, color='#A7ADB5', lw=.6, ls=(0, (3, 3)), zorder=0)
    for i, (genotype, colour, marker) in enumerate((('CTRL', CTRL_BLUE, 'o'), ('LOY', LOY_RED, 's'))):
        for rule in RULES:
            s = mean(rule, 'same_genotype', genotype, 'GFP fraction')
            x = i + offsets[rule] * 2.2
            ax.plot([x, x], [s['min'], s['max']], color=colour, lw=.6, alpha=.7)
            ax.plot([x], [s['mean']], **marker_kwargs(rule, colour, marker, 3.1))
        ax.text(i, 1.15, '\n'.join(signed(mean(r, 'same_genotype', genotype, 'GFP fraction')['mean']) for r in RULES),
                ha='center', va='top', fontsize=5.6, color=colour, linespacing=1.05)
    ax.text(1.5, 1.15, 'standard\nadopted\nuniform', ha='left', va='top', fontsize=5.6, color=GREY, linespacing=1.05,
            clip_on=False)
    ax.set(xlim=(-.6, 1.45), ylim=(-.3, 1.2), xticks=[0, 1], xticklabels=['CTRL/CTRL', 'LOY/LOY'],
           ylabel='Δ ln(GFP:mCherry)')
    tidy(ax)
    ax.tick_params(axis='x', labelsize=6.2)
    ax.text(.5, 1.03, 'Same-genotype mixtures', transform=ax.transAxes, ha='center', va='bottom', fontsize=6.2, color=GREY)

    stem = OUT / 'ch04_supp_competition_unlabelled_gate_all_arms_draft'
    fig.savefig(stem.with_suffix('.pdf'))
    fig.savefig(stem.with_suffix('.svg'))
    fig.savefig(stem.with_suffix('.png'), dpi=600)
    plt.close(fig)
    return stem


def main():
    df, summary, cuts, hashes, inputs = compute()
    stem = render(df, summary, cuts)
    for rel, h in hashes.items():
        if digest(COMP / rel) != h:
            raise ValueError(f'Raw file changed during rendering: {rel}')
    manifest = dict(generated_utc=datetime.now(timezone.utc).isoformat(), raw_sha256=hashes, input_sha256=inputs,
                    script_sha256=digest(__file__), matplotlib=matplotlib.__version__, figure_mm=[WIDTH, HEIGHT], cuts=cuts,
                    reproduced=['standard-gate per-well delta = D4 wells.csv delta_log_odds_standard_gate', 'background-specific per-well delta = unlabelled_gate_sensitivity.json'],
                    exports={p.name: digest(p) for p in (stem.with_suffix(s) for s in ('.pdf', '.svg', '.png'))})
    (DATA / 'manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    print(summary.round(3).to_string())
    print(stem.with_suffix('.png'))


if __name__ == '__main__':
    main()
