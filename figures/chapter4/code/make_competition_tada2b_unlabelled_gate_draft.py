"""TADA2B D0-to-D4 competition with the day-specific unlabelled-cell gate (supplementary draft).

Reads the hash-verified D0/D4 FCS files through the comp_assay A8 pipeline, reproduces the
saved valley cuts and per-well sensitivity values, writes binned source data, and renders a-e.
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
DATA = OUT / 'source_data/competition_tada2b_unlabelled_gate_draft'
REPORT = COMP / 'analysis/spectral_growth_assay_D4_20261020'
os.environ.setdefault('MPLCONFIGDIR', str(DATA / 'mpl_cache'))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from thesis_style import apply_style
from spectral.a8_pipeline.competition_report import load_cef_gates
from spectral.a8_pipeline.fcs_io import load_fcs
from spectral.a8_pipeline.gating import points_in_polygon
from spectral.a8_pipeline.growth_assay_followup import GATE_NAMES
from spectral.a8_pipeline.growth_assay_unlabelled_gate import VALLEY, find_valley, gated_counts

GREEN, CHERRY, GREY, INK = '#32965B', '#D34F4F', '#777E86', '#222222'
ORIENT = {'LOY_GFP': (GREEN, 'o', 'LOY–GFP'), 'LOY_mCherry': (CHERRY, 's', 'LOY–mCherry')}
WIDTH, HEIGHT = 183, 96
EDGES = np.arange(-6, 12.001, .15)
M_EDGES = np.arange(3, 10.001, .1)


def digest(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def signed(v):
    v = round(float(v), 2)
    return '0.00' if v == 0 else f'{v:+.2f}'.replace('-', '−')


def logit(p):
    return np.log(p / (1 - p))


def extract():
    d4 = json.loads((REPORT / 'results.json').read_text(encoding='utf-8'))
    d0 = json.loads((COMP / d4['baseline_path']).read_text(encoding='utf-8'))
    sens_path = REPORT / 'unlabelled_gate_sensitivity.json'
    sens = json.loads(sens_path.read_text(encoding='utf-8'))
    cuts = sens['rule']['mcherry_below']['TADA2B']
    primary_cut = sens['rule']['primary_mcherry_below']
    weights, intercept = d0['classifier']['weights'], d0['classifier']['intercept']
    gfp_below = d0['negative_rules']['broad']['gfp_below']
    features, hashes = {}, {}
    for day, data in ((0, d0), (4, d4)):
        expected = {m['path']: m['sha256'] for m in data['source_manifest']}
        source = COMP / data['source']
        gates = load_cef_gates(next(source.glob('*.cef')))
        raw = {p.stem: p for p in source.rglob('*.fcs')}
        for w in data['wells']:
            if w['role'] != 'mixture' or w['background'] != 'TADA2B':
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
            features[(day, w['well'], w['orientation'])] = np.arcsinh(events.loc[keep, ['eGFP-A', 'mCherry-A']].to_numpy(float) / 1000)
    hist2d, hist1d, valleys, excluded, wells = [], {}, {}, [], []
    for day in (0, 4):
        pooled = np.vstack([f for (d, _, _), f in features.items() if d == day])
        m = pooled[pooled[:, 0] < VALLEY['gfp_below'], 1]
        valleys[day] = find_valley(m)
        if valleys[day]['valley'] != cuts[f'D{day}']:
            raise ValueError('Valley does not reproduce the saved sensitivity cut')
        hist1d[day] = np.histogram(m, bins=M_EDGES)[0] / len(pooled) * 100
        for orientation in ORIENT:
            f = np.vstack([v for (d, _, o), v in features.items() if d == day and o == orientation])
            hist2d.append(dict(day=day, orientation=orientation, counts=np.histogram2d(f[:, 0], f[:, 1], [EDGES, EDGES])[0]))
            for rule, cut in (('primary', primary_cut), ('unlabelled gate', cuts[f'D{day}'])):
                c = [gated_counts(v, weights, intercept, gfp_below, cut if rule == 'primary' else cuts[f'D{day}'])
                     for (d, _, o), v in features.items() if d == day and o == orientation]
                excluded.append(dict(day=day, orientation=orientation, rule=rule,
                                     excluded_percent=100 * (1 - sum(x['retained'] for x in c) / sum(x['primary_singlets'] for x in c))))
    primary = pd.read_csv(REPORT / 'wells.csv')
    primary = primary[(primary.role == 'mixture') & (primary.group == 'TADA2B')]
    sens_wells = {w['well']: w for w in sens['wells']}
    for r in primary.itertuples():
        s = sens_wells[r.well]
        wells.append(dict(well=r.well, orientation=r.orientation, replicate=r.replicate,
                          d0_loy_primary=r.d0_loy_fraction_standard_gate, d4_loy_primary=r.loy_fraction_standard_gate,
                          delta_primary=r.delta_log_odds_standard_gate,
                          d0_loy_gate=s['d0_loy_fraction'], d4_loy_gate=s['loy_fraction'], delta_gate=s['delta_log_odds']))
    wells = pd.DataFrame(wells)
    wt = next(a for a in d4['dye_agreement'] if a['group'] == 'WT')
    wt_mean = (wt['gfp_delta_log_odds'] + wt['mcherry_delta_log_odds']) / 2
    DATA.mkdir(parents=True, exist_ok=True)
    wells.to_csv(DATA / 'tada2b_wells.csv', index=False)
    pd.DataFrame(excluded).to_csv(DATA / 'excluded_fractions.csv', index=False)
    pd.DataFrame({'mcherry_bin_left': M_EDGES[:-1], 'D0_percent': hist1d[0], 'D4_percent': hist1d[4]}).to_csv(DATA / 'gfp_low_mcherry_histogram.csv', index=False)
    np.savez_compressed(DATA / 'reporter_density_counts.npz', edges=EDGES,
                        **{f"D{h['day']}_{h['orientation']}": h['counts'] for h in hist2d})
    inputs = {p: digest(p) for p in (REPORT / 'results.json', sens_path, REPORT / 'wells.csv', COMP / d4['baseline_path'])}
    return dict(hist2d=hist2d, hist1d=hist1d, valleys=valleys, excluded=pd.DataFrame(excluded), wells=wells,
                cuts=cuts, primary_cut=primary_cut, wt_mean=wt_mean, raw_hashes=hashes,
                input_hashes={str(k): v for k, v in inputs.items()})


def render(d):
    apply_style(7.5)
    plt.rcParams['axes.unicode_minus'] = True
    fig = plt.figure(figsize=(WIDTH / 25.4, HEIGHT / 25.4))
    label = lambda x, y, s, **kw: fig.text(x / WIDTH, y / HEIGHT, s, **kw)
    axes = lambda x, y, w, h: fig.add_axes([x / WIDTH, y / HEIGHT, w / WIDTH, h / HEIGHT])

    def tidy(ax):
        ax.tick_params(labelsize=6.3, length=1.6, pad=1.2, width=.55)
        ax.xaxis.labelpad = 1.5
        ax.yaxis.labelpad = 1.5

    # a: pooled reporter densities with primary and day-specific gates.
    label(1.5, 93.5, 'a', fontsize=9, fontweight='semibold', va='center')
    vmax = max(h['counts'].max() for h in d['hist2d'])
    order = [('LOY_GFP', 0), ('LOY_GFP', 4), ('LOY_mCherry', 0), ('LOY_mCherry', 4)]
    for i, (orientation, day) in enumerate(order):
        ax = axes(12 + 28.5 * i, 53, 26, 34)
        counts = next(h['counts'] for h in d['hist2d'] if h['day'] == day and h['orientation'] == orientation)
        ax.pcolormesh(EDGES, EDGES, np.ma.masked_less(counts.T, 1), norm=LogNorm(1, vmax), cmap='magma_r', rasterized=True)
        cut = d['cuts'][f'D{day}']
        ax.plot([-6, 6, 6], [d['primary_cut']] * 2 + [-6], color=GREY, lw=.7, ls=(0, (3, 2)))
        ax.plot([-6, 6, 6], [cut, cut, -6], color=INK, lw=.8)
        ax.set(xlim=(-6, 12), ylim=(-6, 12), xticks=[-5, 0, 5, 10], yticks=[-5, 0, 5, 10])
        tidy(ax)
        ax.set_xlabel('eGFP', fontsize=6.6)
        if i:
            ax.set_yticklabels([])
        else:
            ax.set_ylabel('mCherry', fontsize=6.6)
        colour = ORIENT[orientation][0]
        ax.text(.03, 1.03, f"{ORIENT[orientation][2]} · D{day}", transform=ax.transAxes, fontsize=6.6,
                fontweight='semibold', color=colour, va='bottom')
    label(12, 93.5, 'arcsinh(unmixed signal / 1000); pooled TADA2B mixture wells', fontsize=6.3, color=GREY, va='center')
    fig.legend(handles=[plt.Line2D([], [], color=GREY, lw=.8, ls=(0, (3, 2)), label=f"Standard gate (mCherry < {d['primary_cut']})"),
                        plt.Line2D([], [], color=INK, lw=.9, label='Day-specific unlabelled-cell gate (adopted)')],
               loc='center right', bbox_to_anchor=(181 / WIDTH, 93.5 / HEIGHT), ncol=2, fontsize=6.3, handlelength=1.8,
               columnspacing=1, frameon=False)

    # b: mCherry of GFP-low events, D0 vs D4, with cuts.
    label(129, 88.5, 'b', fontsize=9, fontweight='semibold', va='center')
    ax = axes(138, 53, 43, 34)
    mid = (M_EDGES[:-1] + M_EDGES[1:]) / 2
    for day, colour in ((0, '#9AA0A8'), (4, INK)):
        v = d['valleys'][day]
        ax.plot(mid, d['hist1d'][day], color=colour, lw=1,
                label=f"D{day}: peak {v['unlabelled_peak']:.1f}, cut {d['cuts'][f'D{day}']:.2f}")
        ax.axvline(d['cuts'][f'D{day}'], color=colour, lw=.8, ymax=.8)
    ax.axvline(d['primary_cut'], color=GREY, lw=.7, ls=(0, (3, 2)), ymax=.8)
    ax.set(xlim=(4, 10), ylim=(0, 3.2), xlabel='mCherry (eGFP < 4)', ylabel='Events (% of singlets)')
    tidy(ax)
    ax.legend(loc='upper left', fontsize=5.8, handlelength=1.2, handletextpad=.4, borderaxespad=.1, labelspacing=.2)

    # c: LOY share D0 -> D4, fixed gate (filled) and primary (faint).
    label(1.5, 44, 'c', fontsize=9, fontweight='semibold', va='center')
    ax = axes(12, 9, 34, 31)
    w = d['wells']
    jitter = np.array([-.08, 0, .08])
    for orientation, (colour, marker, _) in ORIENT.items():
        x = w[w.orientation == orientation].sort_values('replicate')
        for prefix, alpha, mfc, ls in (('primary', .35, 'white', (0, (2, 1.5))), ('gate', 1, None, '-')):
            d0, d4 = 100 * x[f'd0_loy_{prefix}'].iloc[0], 100 * x[f'd4_loy_{prefix}'].to_numpy()
            for xx, yy in zip(1 + jitter, d4):
                ax.plot([0, xx], [d0, yy], color=colour, lw=.5, alpha=.6 * alpha, ls=ls)
            ax.plot([0], [d0], marker, color=colour, ms=3.2, alpha=alpha, mfc=mfc or colour, mew=.6)
            ax.plot(1 + jitter, d4, marker, color=colour, ms=2.8, alpha=alpha, mfc=mfc or 'white', mew=.6)
    ax.set(xlim=(-.35, 1.35), ylim=(20, 70), xticks=[0, 1], xticklabels=['D0', 'D4'], yticks=[20, 40, 60],
           ylabel='LOY of retained cells (%)')
    tidy(ax)

    # d: change in log-odds per well, primary vs unlabelled-cell gate, with WT reference.
    label(55, 44, 'd', fontsize=9, fontweight='semibold', va='center')
    ax = axes(66, 9, 70, 31)
    ax.axhline(0, color='#A7ADB5', lw=.6, ls=(0, (3, 3)), zorder=0)
    ax.axhline(d['wt_mean'], color=GREY, lw=.7, ls=(0, (1, 1.5)), zorder=0)
    ax.text(3.3, d['wt_mean'] + .03, f"WT mean {signed(d['wt_mean'])}", ha='right', va='bottom', fontsize=6, color=GREY)
    means = {}
    for i, (orientation, (colour, marker, name)) in enumerate(ORIENT.items()):
        x = w[w.orientation == orientation]
        for offset, key, mfc, alpha in ((-.17, 'delta_primary', 'white', .45), (.17, 'delta_gate', colour, 1)):
            ax.plot(i + offset + jitter * .5, x[key], marker, color=colour, ms=3, mfc=mfc, mew=.6, alpha=alpha)
            means[(orientation, key)] = x[key].mean()
            ax.plot([i + offset - .1, i + offset + .1], [means[(orientation, key)]] * 2, color=colour, lw=1, alpha=alpha)
        ax.text(i - .17, .82, signed(means[(orientation, 'delta_primary')]), ha='center', fontsize=6.1, color=GREY)
        ax.text(i + .17, .82, signed(means[(orientation, 'delta_gate')]), ha='center', fontsize=6.1, color=colour)
    for offset, key, mfc in ((-.17, 'delta_primary', 'white'), (.17, 'delta_gate', INK)):
        m = (means[('LOY_GFP', key)] + means[('LOY_mCherry', key)]) / 2
        ax.plot([2 + offset], [m], 'D', color=INK, ms=3.4, mfc=mfc, mew=.7)
    m0 = (means[('LOY_GFP', 'delta_primary')] + means[('LOY_mCherry', 'delta_primary')]) / 2
    m1 = (means[('LOY_GFP', 'delta_gate')] + means[('LOY_mCherry', 'delta_gate')]) / 2
    ax.text(2 - .17, .82, signed(m0), ha='center', fontsize=6.1, color=GREY)
    ax.text(2 + .17, .82, signed(m1), ha='center', fontsize=6.1, color=INK)
    ax.set(xlim=(-.55, 3.35), ylim=(-.45, .95), xticks=[0, 1, 2], xticklabels=['LOY–GFP', 'LOY–mCherry', 'Mean of\norientations'],
           ylabel='Δ ln(LOY:CTRL), D4 − D0')
    tidy(ax)
    ax.tick_params(axis='x', labelsize=6.2)
    ax.legend(handles=[plt.Line2D([], [], color=GREY, marker='o', ls='none', ms=3, mfc='white', mew=.6, label='Standard gate'),
                       plt.Line2D([], [], color=GREY, marker='o', ls='none', ms=3, mew=.6, label='Unlabelled-cell gate (adopted)')],
              loc='lower right', fontsize=6, handletextpad=.2, borderaxespad=.2, labelspacing=.25)

    # e: reporter-negative exclusion under each rule.
    label(143, 44, 'e', fontsize=9, fontweight='semibold', va='center')
    ax = axes(153, 9, 28, 31)
    ex = d['excluded']
    xs = {(0, 'LOY_GFP'): 0, (0, 'LOY_mCherry'): 1, (4, 'LOY_GFP'): 2.4, (4, 'LOY_mCherry'): 3.4}
    for (day, orientation), x in xs.items():
        colour = ORIENT[orientation][0]
        for offset, rule, fill in ((-.19, 'primary', 'white'), (.19, 'unlabelled gate', colour)):
            v = ex[(ex.day == day) & (ex.orientation == orientation) & (ex.rule == rule)].excluded_percent.iloc[0]
            ax.bar(x + offset, v, width=.36, color=fill, edgecolor=colour, lw=.6)
    ax.set(xlim=(-.6, 4), ylim=(0, 40), xticks=[.5, 2.9], xticklabels=['D0', 'D4'], ylabel='Excluded (% of singlets)',
           yticks=[0, 10, 20, 30, 40])
    tidy(ax)
    ax.text(.5, 1.03, 'open: standard · filled: adopted', transform=ax.transAxes, fontsize=5.8, ha='center', va='bottom',
            color=GREY)

    stem = OUT / 'ch04_supp_competition_tada2b_unlabelled_gate_draft'
    fig.savefig(stem.with_suffix('.pdf'))
    fig.savefig(stem.with_suffix('.svg'))
    fig.savefig(stem.with_suffix('.png'), dpi=600)
    plt.close(fig)
    return stem


def main():
    d = extract()
    stem = render(d)
    for rel, h in d['raw_hashes'].items():
        if digest(COMP / rel) != h:
            raise ValueError(f'Raw file changed during rendering: {rel}')
    manifest = dict(generated_utc=datetime.now(timezone.utc).isoformat(), raw_sha256=d['raw_hashes'],
                    input_sha256=d['input_hashes'], script_sha256=digest(__file__), matplotlib=matplotlib.__version__,
                    figure_mm=[WIDTH, HEIGHT], cuts=d['cuts'], valleys={f'D{k}': v for k, v in d['valleys'].items()},
                    exports={p.name: digest(p) for p in (stem.with_suffix(s) for s in ('.pdf', '.svg', '.png'))},
                    note='Sensitivity display. Valley cuts reproduce unlabelled_gate_sensitivity.json; per-well values are copied '
                         'from that file and from the standard-gate columns of the D4 report. No inferential test.')
    (DATA / 'manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    print(stem.with_suffix('.png'))


if __name__ == '__main__':
    main()
