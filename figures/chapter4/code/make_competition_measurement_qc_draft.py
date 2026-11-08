"""Supplementary draft: measurement checks behind the competition-assay figure.

a-b: row-C same-genotype identity check (per-cell size/brightness classifier trained on WT
competition cells of known genotype, cross-fitted between rows A and B).
c-d: the TADA2B unlabelled-cell population and the adopted day-specific gate.
e: D0-to-D4 change for every arm under the standard gate, the adopted gating and the
TADA2B cuts applied to every arm.
Reads hash-verified FCS files through the comp_assay A8 pipeline; no inferential test.
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
DATA = OUT / 'source_data/competition_measurement_qc_draft'
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
COLOURS = json.loads((THESIS / 'skills/thesis-figure-style/assets/condition_colours.json').read_text())['palette']
CTRL_BLUE, LOY_RED = COLOURS['Y-retaining'], COLOURS['Y-depleted']
SIZE = ['LightLoss (Violet)-A', 'SSC (Imaging)-A', 'Size (FSC)', 'Size (LightLoss (Imaging))']
ARMS = ('WT', 'WT_C1', 'DDX3X', 'DDX3X_C1', 'SREBF1', 'TADA2B')
NAMES = {'WT': 'WT', 'WT_C1': 'WT C1', 'DDX3X': 'DDX3X', 'DDX3X_C1': 'DDX3X C1', 'SREBF1': 'SREBF1', 'TADA2B': 'TADA2B'}
ORIENT = {'LOY_GFP': (GREEN, 'o', 'LOY–GFP'), 'LOY_mCherry': (CHERRY, 's', 'LOY–mCherry')}
RULES = ('standard', 'adopted', 'uniform')
RULE_LABEL = {'standard': 'Standard gate (mCherry < 6.5)', 'adopted': 'Adopted gating (TADA2B: unlabelled-cell gate)',
              'uniform': 'TADA2B cuts applied to every arm'}
SIZE_EDGES = np.linspace(7.1, 7.9, 41)
REP_EDGES = np.arange(-6, 12.001, .15)
M_EDGES = np.arange(3, 10.001, .1)
WIDTH, HEIGHT = 183, 126


def digest(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def signed(v):
    v = round(float(v), 2)
    return '0.00' if v == 0 else f'{v:+.2f}'.replace('-', '−')


def logit(p):
    return float(np.log(p / (1 - p)))


def lda(a, b):
    ma, mb = a.mean(0), b.mean(0)
    w = np.linalg.solve(np.cov(np.vstack([a - ma, b - mb]).T), ma - mb)
    return w, -w @ (ma + mb) / 2


def load():
    d4 = json.loads((REPORT / 'results.json').read_text(encoding='utf-8'))
    d0 = json.loads((COMP / d4['baseline_path']).read_text(encoding='utf-8'))
    sens = json.loads((REPORT / 'unlabelled_gate_sensitivity.json').read_text(encoding='utf-8'))
    if d4['identity_status'] != 'user_confirmed' or sens['identity_mapping_sha256'] != d4['identity_mapping_sha256']:
        raise ValueError('Expected approved D4 labels and a matching sensitivity file')
    cells, hashes = {}, {}
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
            e = events.loc[keep]
            cells[(day, w['well'])] = dict(meta=w, f=np.arcsinh(e[['eGFP-A', 'mCherry-A']].to_numpy(float) / 1000),
                                           size=e[SIZE].to_numpy(float))
    return d0, d4, sens, cells, hashes


def identity(d0, cells):
    """Row-C identity check: estimated LOY content of each well, as in the approved review."""
    weights, intercept = np.array(d0['classifier']['weights']), d0['classifier']['intercept']
    rule = d0['negative_rules']['broad']
    rows, hist = [], {}
    for day in (0, 4):
        prepared = {}
        for (dd, well), c in cells.items():
            meta = c['meta']
            if dd != day or not (meta['role'] == 'dye_neutral' or meta['group'] in ('WT', 'WT_C1')):
                continue
            f, size = c['f'], c['size']
            ok = (size > 0).all(1)
            neg = (f[:, 0] < rule['gfp_below']) & (f[:, 1] < rule['mcherry_below'])
            gfp = ~neg & (f @ weights + intercept > 0) & ok
            mch = ~neg & ~gfp & ok
            X = np.column_stack([np.log10(np.where(ok[:, None], size, 1)), f])
            prepared[well] = dict(meta=meta, X=X, gfp=gfp, mch=mch)

        def known(row_letter, cls):
            out = {'LOY': [], 'CTRL': []}
            for well, p in prepared.items():
                if well[0] == row_letter and p['meta']['role'] == 'mixture':
                    loy_cls = 'gfp' if p['meta']['orientation'] == 'LOY_GFP' else 'mch'
                    out['LOY' if cls == loy_cls else 'CTRL'].append(p['X'][p[cls]])
            return {k: np.vstack(v) for k, v in out.items()}
        estimates = {}
        for cls in ('gfp', 'mch'):
            cols = [0, 1, 2, 3, 4 if cls == 'gfp' else 5]
            for train, cal in (('A', 'B'), ('B', 'A')):
                tr, ca = known(train, cls), known(cal, cls)
                w, b = lda(tr['LOY'][:, cols], tr['CTRL'][:, cols])
                r_loy, r_ctrl = (ca['LOY'][:, cols] @ w + b > 0).mean(), (ca['CTRL'][:, cols] @ w + b > 0).mean()
                for well, p in prepared.items():
                    if well[0] == 'C':
                        r = (p['X'][p[cls]][:, cols] @ w + b > 0).mean()
                        estimates.setdefault(well, []).append((r - r_ctrl) / (r_loy - r_ctrl))
            if cls == 'gfp':
                both = {g: np.vstack([known('A', cls)[g], known('B', cls)[g]])[:, 0] for g in ('LOY', 'CTRL')}
                both['C01–C03'] = np.concatenate([prepared[x]['X'][prepared[x]['gfp']][:, 0] for x in ('C01', 'C02', 'C03')])
                both['C04–C06'] = np.concatenate([prepared[x]['X'][prepared[x]['gfp']][:, 0] for x in ('C04', 'C05', 'C06')])
                hist[day] = {k: np.histogram(v, SIZE_EDGES)[0] / len(v) * 100 for k, v in both.items()}
        for well, v in sorted(estimates.items()):
            meta = prepared[well]['meta']
            recorded = meta['genotype'] if day == 0 else meta['nominal_genotype']
            rows.append(dict(day=day, well=well, recorded=recorded, analysed=meta['genotype'],
                             loy_content=100 * np.mean(v), low=100 * min(v), high=100 * max(v)))
    return pd.DataFrame(rows), hist


def gates(d0, sens, cells):
    weights, intercept = d0['classifier']['weights'], d0['classifier']['intercept']
    gfp_below = d0['negative_rules']['broad']['gfp_below']
    rule = sens['rule']
    cut = {'standard': lambda bg, day: rule['primary_mcherry_below'],
           'adopted': lambda bg, day: rule['mcherry_below'][bg][f'D{day}'] if bg in rule['adjusted_backgrounds'] else rule['primary_mcherry_below'],
           'uniform': lambda bg, day: rule['mcherry_below']['TADA2B'][f'D{day}']}
    rows = []
    for r in RULES:
        counts = {k: gated_counts(c['f'], weights, intercept, gfp_below, cut[r](c['meta']['background'], k[0]))
                  for k, c in cells.items() if c['meta']['role'] == 'mixture'}
        base = {}
        for b in d0['d0_baselines']:
            wells = b['wells'].split(';')
            kept = sum(counts[(0, x)]['retained'] for x in wells)
            green = sum(counts[(0, x)]['gfp_count'] for x in wells)
            base[(b['group'], b['orientation'])] = (green if b['orientation'] == 'LOY_GFP' else kept - green) / kept
        for (day, well), c in counts.items():
            if day != 4:
                continue
            meta = cells[(day, well)]['meta']
            gfp = c['gfp_count'] / c['retained']
            loy = gfp if meta['orientation'] == 'LOY_GFP' else 1 - gfp
            rows.append(dict(rule=r, group=meta['group'], orientation=meta['orientation'], well=well,
                             delta=logit(loy) - logit(base[(meta['group'], meta['orientation'])])))
    df = pd.DataFrame(rows)
    primary = pd.read_csv(REPORT / 'wells.csv').set_index('well')
    for row in df.itertuples():
        expected = {'standard': primary.loc[row.well, 'delta_log_odds_standard_gate'],
                    'adopted': primary.loc[row.well, 'delta_log_odds']}.get(row.rule)
        if expected is not None and abs(expected - row.delta) > 1e-9:
            raise ValueError(f'{row.rule} result not reproduced for {row.well}')
    return df


def tada2b(sens, cells):
    hist2d, hist1d, valleys = {}, {}, {}
    for day in (0, 4):
        f = np.vstack([c['f'] for (d, _), c in cells.items() if d == day and c['meta']['role'] == 'mixture'
                       and c['meta']['background'] == 'TADA2B'])
        m = f[f[:, 0] < VALLEY['gfp_below'], 1]
        valleys[day] = find_valley(m)
        if valleys[day]['valley'] != sens['rule']['mcherry_below']['TADA2B'][f'D{day}']:
            raise ValueError('TADA2B valley does not reproduce the saved cut')
        hist2d[day] = np.histogram2d(f[:, 0], f[:, 1], [REP_EDGES, REP_EDGES])[0]
        hist1d[day] = np.histogram(m, bins=M_EDGES)[0] / len(f) * 100
    return hist2d, hist1d, valleys


def render(ident, size_hist, gate_df, t2, sens):
    apply_style(7.5)
    plt.rcParams['axes.unicode_minus'] = True
    fig = plt.figure(figsize=(WIDTH / 25.4, HEIGHT / 25.4))
    label = lambda x, y, s, **kw: fig.text(x / WIDTH, y / HEIGHT, s, **kw)
    axes = lambda x, y, w, h: fig.add_axes([x / WIDTH, y / HEIGHT, w / WIDTH, h / HEIGHT])

    def tidy(ax):
        ax.tick_params(labelsize=6.2, length=1.6, pad=1.2, width=.55)
        ax.xaxis.labelpad = 1.5
        ax.yaxis.labelpad = 1.5

    # a: estimated LOY content of the row-C same-genotype wells.
    label(1.5, 122.5, 'a', fontsize=9, fontweight='semibold', va='center')
    ax = axes(14, 94, 74, 24)
    x = 0
    ticks, ticklabels = [], []
    for day in (0, 4):
        for r in ident[ident.day == day].itertuples():
            colour = LOY_RED if r.recorded == 'LOY' else CTRL_BLUE
            ax.bar(x, max(r.loy_content, .8), width=.7, color=colour, edgecolor=colour, lw=.5,
                   hatch='////' if r.recorded != r.analysed else None, alpha=.9)
            ax.plot([x, x], [r.low, r.high], color=INK, lw=.6)
            ticks.append(x)
            ticklabels.append(r.well)
            x += 1
        ax.text(x - 3.5, 118, f'D{day}', ha='center', fontsize=6.6, fontweight='semibold')
        x += .8
    ax.axhline(100, color='#A7ADB5', lw=.5, ls=(0, (3, 3)))
    ax.set(xlim=(-.7, x - .5), ylim=(0, 125), yticks=[0, 50, 100], xticks=ticks, xticklabels=ticklabels,
           ylabel='Estimated LOY content (%)')
    tidy(ax)
    ax.tick_params(axis='x', labelsize=5.8, rotation=0)
    from matplotlib.patches import Patch
    fig.legend(handles=[Patch(color=CTRL_BLUE, label='Recorded CTRL/CTRL'), Patch(color=LOY_RED, label='Recorded LOY/LOY'),
                        Patch(facecolor='white', edgecolor=INK, hatch='////', label='Label corrected')],
               loc='center', bbox_to_anchor=(51 / WIDTH, 84.5 / HEIGHT), ncol=3, fontsize=5.9, handlelength=1.2,
               columnspacing=1, frameon=False)

    # b: GFP+ cell size of row-C wells against known LOY and CTRL cells.
    label(97, 122.5, 'b', fontsize=9, fontweight='semibold', va='center')
    mid = 10 ** ((SIZE_EDGES[:-1] + SIZE_EDGES[1:]) / 2) / 1e6
    for i, day in enumerate((0, 4)):
        ax = axes(107 + 38 * i, 94, 35, 24)
        h = size_hist[day]
        ax.fill_between(mid, h['LOY'], color=LOY_RED, alpha=.25, lw=0, label='Known LOY')
        ax.fill_between(mid, h['CTRL'], color=CTRL_BLUE, alpha=.25, lw=0, label='Known CTRL')
        ax.plot(mid, h['C01–C03'], color=INK, lw=.9, label='C01–C03')
        ax.plot(mid, h['C04–C06'], color=INK, lw=.9, ls=(0, (3, 1.6)), label='C04–C06')
        ax.set(xscale='log', xlim=(mid[0], mid[-1]), xticks=[15, 30, 60], xticklabels=['15', '30', '60'], yticks=[],
               xlabel='GFP⁺ cell size, LightLoss (Violet)-A (×10⁶)' if i == 0 else '')
        if i == 0:
            ax.xaxis.set_label_coords(1.07, -.17)
        ax.xaxis.set_minor_formatter(plt.NullFormatter())
        tidy(ax)
        ax.spines['left'].set_visible(False)
        recorded = 'LOY' if day == 0 else 'CTRL'
        ax.text(.02, 1.03, f'D{day} · C01–C03 recorded {recorded}', transform=ax.transAxes, fontsize=5.9,
                va='bottom', color=GREY)
        if i == 1:
            fig.legend(*ax.get_legend_handles_labels(), loc='center', bbox_to_anchor=(144 / WIDTH, 84.5 / HEIGHT), ncol=4,
                       fontsize=5.9, handlelength=1.4, columnspacing=1, frameon=False)

    # c: TADA2B reporter density with both gates; d: mCherry of GFP-low events.
    label(1.5, 78, 'c', fontsize=9, fontweight='semibold', va='center')
    hist2d, hist1d, valleys = t2
    vmax = max(h.max() for h in hist2d.values())
    cuts = sens['rule']['mcherry_below']['TADA2B']
    std = sens['rule']['primary_mcherry_below']
    for i, day in enumerate((0, 4)):
        ax = axes(14 + 33 * i, 47, 29, 27)
        ax.pcolormesh(REP_EDGES, REP_EDGES, np.ma.masked_less(hist2d[day].T, 1), norm=LogNorm(1, vmax), cmap='magma_r',
                      rasterized=True)
        ax.plot([-6, 6, 6], [std, std, -6], color=GREY, lw=.7, ls=(0, (3, 2)))
        ax.plot([-6, 6, 6], [cuts[f'D{day}']] * 2 + [-6], color=INK, lw=.8)
        ax.set(xlim=(-6, 12), ylim=(-6, 12), xticks=[-5, 0, 5, 10], yticks=[-5, 0, 5, 10], xlabel='eGFP')
        tidy(ax)
        if i:
            ax.set_yticklabels([])
        else:
            ax.set_ylabel('mCherry')
        ax.text(.03, 1.03, f'TADA2B · D{day}', transform=ax.transAxes, fontsize=6.4, fontweight='semibold', va='bottom')
    label(97, 78, 'd', fontsize=9, fontweight='semibold', va='center')
    ax = axes(107, 47, 74, 27)
    m = (M_EDGES[:-1] + M_EDGES[1:]) / 2
    for day, colour in ((0, '#9AA0A8'), (4, INK)):
        v = valleys[day]
        ax.plot(m, hist1d[day], color=colour, lw=1,
                label=f"D{day}: unlabelled peak {v['unlabelled_peak']:.1f}, cut {cuts[f'D{day}']:.2f}")
        ax.axvline(cuts[f'D{day}'], color=colour, lw=.8, ymax=.56)
    ax.axvline(std, color=GREY, lw=.7, ls=(0, (3, 2)), ymax=.56, label=f'Standard cut {std}')
    ax.set(xlim=(4, 10), ylim=(0, 3.3), xlabel='TADA2B mCherry (eGFP < 4)', ylabel='% of singlets')
    tidy(ax)
    ax.legend(loc='upper left', fontsize=5.8, handlelength=1.4, labelspacing=.2, borderaxespad=.1)

    # e: every arm under three gates.
    label(1.5, 37, 'e', fontsize=9, fontweight='semibold', va='center')
    ax = axes(14, 8, 167, 25)
    ax.axhline(0, color='#A7ADB5', lw=.6, ls=(0, (3, 3)), zorder=0)
    offsets = dict(zip(RULES, (-.075, 0, .075)))
    style = {'standard': dict(mfc='white'), 'adopted': dict(), 'uniform': dict(fillstyle='bottom', mfcalt='white')}
    for i, arm in enumerate(ARMS):
        if i:
            ax.axvline(i - .5, color='#E3E5E8', lw=.5, zorder=0)
        for j, (orientation, (colour, marker, _)) in enumerate(ORIENT.items()):
            for r in RULES:
                v = gate_df[(gate_df.rule == r) & (gate_df.group == arm) & (gate_df.orientation == orientation)].delta
                x = i + (-.2 if j == 0 else .2) + offsets[r]
                ax.plot([x, x], [v.min(), v.max()], color=colour, lw=.6, alpha=.7)
                kw = dict(color=colour, ms=3, mew=.6, mfc=colour)
                kw.update(style[r])
                ax.plot([x], [v.mean()], marker, ls='none', **kw)
    ax.set(xlim=(-.55, 5.55), ylim=(-.45, 1.35), xticks=range(6), xticklabels=[NAMES[a] for a in ARMS],
           ylabel='Δ ln(LOY:CTRL), D4 − D0')
    tidy(ax)
    handles = []
    for r in RULES:
        kw = dict(color=GREY, ms=3.2, mew=.6, mfc=GREY)
        kw.update(style[r])
        handles.append(plt.Line2D([], [], marker='o', ls='none', label=RULE_LABEL[r], **kw))
    fig.legend(handles=handles, loc='center right', bbox_to_anchor=(181 / WIDTH, 37 / HEIGHT), ncol=3, fontsize=5.9,
               handletextpad=.3, columnspacing=1.2, frameon=False)

    stem = OUT / 'ch04_supp_competition_measurement_qc_draft'
    fig.savefig(stem.with_suffix('.pdf'))
    fig.savefig(stem.with_suffix('.svg'))
    fig.savefig(stem.with_suffix('.png'), dpi=600)
    plt.close(fig)
    return stem


def main():
    d0, d4, sens, cells, hashes = load()
    ident, size_hist = identity(d0, cells)
    gate_df = gates(d0, sens, cells)
    t2 = tada2b(sens, cells)
    DATA.mkdir(parents=True, exist_ok=True)
    ident.to_csv(DATA / 'row_c_identity.csv', index=False)
    gate_df.to_csv(DATA / 'gate_comparison_wells.csv', index=False)
    pd.DataFrame({'mcherry_bin_left': M_EDGES[:-1], 'D0_percent': t2[1][0], 'D4_percent': t2[1][4]}).to_csv(
        DATA / 'tada2b_gfp_low_mcherry_histogram.csv', index=False)
    np.savez_compressed(DATA / 'binned_distributions.npz', reporter_edges=REP_EDGES, size_edges=SIZE_EDGES,
                        **{f'tada2b_D{k}': v for k, v in t2[0].items()},
                        **{f'size_D{d}_{k}': v for d, h in size_hist.items() for k, v in h.items()})
    stem = render(ident, size_hist, gate_df, t2, sens)
    for rel, h in hashes.items():
        if digest(COMP / rel) != h:
            raise ValueError(f'Raw file changed during rendering: {rel}')
    manifest = dict(generated_utc=datetime.now(timezone.utc).isoformat(), raw_sha256=hashes,
                    input_sha256={n: digest(REPORT / n) for n in ('results.json', 'wells.csv', 'unlabelled_gate_sensitivity.json')},
                    script_sha256=digest(__file__), matplotlib=matplotlib.__version__, figure_mm=[WIDTH, HEIGHT],
                    tada2b_valleys={f'D{k}': v for k, v in t2[2].items()},
                    exports={p.name: digest(p) for p in (stem.with_suffix(s) for s in ('.pdf', '.svg', '.png'))})
    (DATA / 'manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    print(ident.round(1).to_string(index=False))
    print(stem.with_suffix('.png'))


if __name__ == '__main__':
    main()
