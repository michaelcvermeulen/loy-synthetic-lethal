"""Inspect same-genotype controls at D0 and D4 with the approved frozen analysis.

Collect with the flow environment; render with a Matplotlib environment.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
from hashlib import sha256
import json
import math
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
REPORT = REPO / 'analysis/spectral_growth_assay_D4_20261020'
OUTPUT = REPORT / 'same_genotype'
COLORS = ('#6b7280', '#a8329b', '#1f7a4d')  # excluded, mCherry, GFP


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def collect():
    import numpy as np
    from .competition_report import LDAModel, load_cef_gates
    from .fcs_io import load_fcs
    from .gating import points_in_polygon
    from .growth_assay_followup import GATE_NAMES, measure_reporters, read_statistics

    d4path = REPORT / 'results.json'
    d4 = json.loads(d4path.read_text(encoding='utf-8'))
    d0path = REPO / d4['baseline_path']
    d0 = json.loads(d0path.read_text(encoding='utf-8'))
    if digest(d0path) != d4['baseline_sha256'] or d4['identity_status'] != 'user_confirmed':
        raise ValueError('Expected the verified D0 baseline and approved D4 labels')
    model = LDAModel(np.array(d4['classifier']['weights']), d4['classifier']['intercept'])
    source_hashes = {str(p.relative_to(REPO)): digest(p) for p in (d0path, d4path)}
    samples = []
    for day, data in ((0, d0), (4, d4)):
        expected_hashes = {m['path']: m['sha256'] for m in data['source_manifest']}
        source = REPO / data['source']
        cef = next(source.glob('*.cef'))
        stat_file = next(source.glob('*Statistics.csv'))
        paths = [cef, stat_file]
        raw = {p.stem: p for p in source.rglob('*.fcs')}
        wells = sorted((w for w in data['wells'] if w['role'] == 'dye_neutral'),
                       key=lambda w: (w['genotype'], w['replicate']))
        paths.extend(raw[w['sample'] if day == 0 else w['original_sample']] for w in wells)
        for path in paths:
            rel = path.relative_to(REPO).as_posix()
            actual = digest(path)
            if actual != expected_hashes[rel]:
                raise ValueError(f'Input hash differs from approved analysis: {rel}')
            source_hashes[rel] = actual
        gates = load_cef_gates(cef)
        statistics = read_statistics(stat_file)
        for w in wells:
            path = raw[w['sample'] if day == 0 else w['original_sample']]
            events = load_fcs(path).events
            keep = np.ones(len(events), dtype=bool)
            counts = {'All Events': len(events)}
            for name in GATE_NAMES:
                gate = gates[name]
                keep &= points_in_polygon(events[list(gate.parameters)].to_numpy(float), gate.vertices)
                counts[name] = int(keep.sum())
            if counts != statistics[path.stem.split('-', 1)[1]]:
                raise ValueError(f'Instrument gate counts differ: {path.name}')
            features = np.arcsinh(events.loc[keep, ['eGFP-A', 'mCherry-A']].to_numpy(float)/1000)
            measures = measure_reporters(features, d0)
            for key in ('primary_singlets', 'reporter_positive_singlets', 'gfp_count', 'gfp_fraction'):
                if not math.isclose(measures[key], w[key], abs_tol=1e-12):
                    raise ValueError(f'Measurements differ from saved report: {path.name} {key}')
            rule = d0['negative_rules']['broad']
            excluded = (features[:, 0] < rule['gfp_below']) & (features[:, 1] < rule['mcherry_below'])
            classes = np.where(excluded, 0, np.where(model.predict(features), 2, 1))
            # Store class membership before coordinate rounding. Display all gated events.
            scatter = events.loc[keep, list(gates['Scatter'].parameters)].to_numpy(float)/1000
            points = np.column_stack([features, scatter, classes])
            samples.append(dict(day=day, genotype=w['genotype'], replicate=w['replicate'], well=w['well'],
                                sample=w['sample'] if day == 0 else w['analysis_sample'],
                                original_filename=path.name, source_file=path.relative_to(REPO).as_posix(),
                                measures=measures, gate_counts=counts,
                                points=np.round(points, 5).tolist()))
    summaries = []
    for genotype in ('CTRL', 'LOY'):
        initial = [s for s in samples if s['genotype'] == genotype and s['day'] == 0]
        final = [s for s in samples if s['genotype'] == genotype and s['day'] == 4]
        pooled = sum(s['measures']['gfp_count'] for s in initial)/sum(s['measures']['reporter_positive_singlets'] for s in initial)
        values = [s['measures']['gfp_fraction'] for s in final]
        summaries.append(dict(genotype=genotype, d0_gfp_pooled=pooled, d4_gfp_mean=float(np.mean(values)),
                              d4_gfp_sd=float(np.std(values, ddof=1)), change_pp=100*(np.mean(values)-pooled)))
    all_points = np.concatenate([np.array(s['points']) for s in samples])
    reporter_limits = [[math.floor(all_points[:, i].min()), math.ceil(all_points[:, i].max())] for i in (0, 1)]
    scatter_limits = [[min(0, math.floor(all_points[:, i].min()/100)*100), math.ceil(all_points[:, i].max()/100)*100] for i in (2, 3)]
    payload = dict(created_utc=datetime.now(timezone.utc).isoformat(), classifier=d4['classifier'],
                   broad_rule=d0['negative_rules']['broad'], samples=samples, summaries=summaries,
                   reporter_limits=reporter_limits, scatter_limits=scatter_limits,
                   scatter_parameters=list(gates['Scatter'].parameters),
                   point_columns=['asinh_eGFP', 'asinh_mCherry', 'LightLoss_Violet_A_thousands', 'SSC_Imaging_A_thousands', 'class'],
                   class_codes={'0': 'broad_negative_excluded', '1': 'mCherry_like_retained', '2': 'GFP_like_retained'},
                   source_hashes=source_hashes)
    for rel, expected in source_hashes.items():
        if digest(REPO / rel) != expected:
            raise ValueError(f'Source changed during read: {rel}')
    OUTPUT.mkdir(exist_ok=True)
    (OUTPUT / 'data.json').write_text(json.dumps(payload, separators=(',', ':'), allow_nan=False), encoding='utf-8')
    rows = [{**{k: s[k] for k in ('day', 'genotype', 'replicate', 'well', 'sample', 'original_filename')},
             **s['measures']} for s in samples]
    with (OUTPUT / 'wells.csv').open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps(dict(samples=len(samples), singlets=sum(len(s['points']) for s in samples),
                         summaries=summaries, reporter_limits=reporter_limits, scatter_limits=scatter_limits,
                         measurements_reproduced=True, raw_hashes_unchanged=True), indent=2))


def render():
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Rectangle

    path = OUTPUT / 'data.json'
    data = json.loads(path.read_text(encoding='utf-8'))
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'svg.fonttype': 'none', 'pdf.fonttype': 42,
                         'axes.spines.top': False, 'axes.spines.right': False, 'font.size': 9})
    figdir = OUTPUT / 'figures'
    figdir.mkdir(exist_ok=True)
    paths = []
    for mode in ('reporters', 'scatter'):
        limits = data['reporter_limits'] if mode == 'reporters' else data['scatter_limits']
        fig, axes = plt.subplots(4, 3, figsize=(12, 13.8), sharex=True, sharey=True)
        fig.subplots_adjust(left=.075, right=.975, top=.855, bottom=.075, hspace=.39, wspace=.18)
        title = 'GFP–mCherry distributions' if mode == 'reporters' else 'Light loss / side scatter'
        fig.suptitle('Same-genotype mixtures · ' + title, x=.075, y=.98, ha='left', fontsize=20,
                     color='#172d24', fontweight='bold')
        fig.text(.075, .948, 'D0 and D4 · all gated singlets · shared axes · corrected D4 labels (row-C genotype swap applied)', color='#52635b', fontsize=11)
        handles = [Line2D([], [], color=c, marker='o', linestyle='none', markersize=5, label=t)
                   for c, t in zip(COLORS, ('Excluded broad negatives', 'mCherry-like retained', 'GFP-like retained'))]
        if mode == 'reporters':
            handles.append(Line2D([], [], color='#172d24', linestyle='--', label='Frozen D0 classifier'))
        fig.legend(handles=handles, loc='upper left', bbox_to_anchor=(.067, .93), ncol=4, frameon=False)
        for row, (genotype, day) in enumerate((('CTRL', 0), ('CTRL', 4), ('LOY', 0), ('LOY', 4))):
            for col, rep in enumerate((1, 2, 3)):
                s = next(s for s in data['samples'] if (s['genotype'], s['day'], s['replicate']) == (genotype, day, rep))
                ax = axes[row, col]
                points = np.array(s['points'])
                x, y = (0, 1) if mode == 'reporters' else (2, 3)
                ax.set_xlim(*limits[0]); ax.set_ylim(*limits[1])
                if mode == 'reporters':
                    gx, gy = data['broad_rule']['gfp_below'], data['broad_rule']['mcherry_below']
                    ax.add_patch(Rectangle((limits[0][0], limits[1][0]), gx-limits[0][0], gy-limits[1][0],
                                           color='#e2e5e8', alpha=.45, zorder=0))
                    ax.plot([limits[0][0], gx, gx], [gy, gy, limits[1][0]], color='#6b7280', lw=.8)
                for cls in (0, 1, 2):
                    subset = points[points[:, 4] == cls]
                    ax.scatter(subset[:, x], subset[:, y], s=1.8, alpha=.22 if cls else .5,
                               c=COLORS[cls], linewidths=0, rasterized=True)
                if mode == 'reporters':
                    weights, intercept = data['classifier']['weights'], data['classifier']['intercept']
                    ys = np.array(limits[1])
                    ax.plot(-(weights[1]*ys+intercept)/weights[0], ys, color='#172d24', ls='--', lw=1)
                m = s['measures']
                ax.set_title(f'{genotype}/{genotype} · D{day} · R{rep} · {s["well"]}', loc='left',
                             color='#172d24', fontweight='bold', fontsize=10, pad=19)
                ax.text(0, 1.025, f'GFP {100*m["gfp_fraction"]:.1f}%  |  excluded {100*m["excluded_broad_fraction"]:.1f}%  |  n = {m["primary_singlets"]:,}',
                        transform=ax.transAxes, color='#52635b', fontsize=8)
                ax.tick_params(labelbottom=True, labelleft=True, colors='#52635b', labelsize=8)
                if col == 0:
                    ax.set_ylabel('mCherry-A · asinh(x/1000)' if mode == 'reporters' else data['scatter_parameters'][1] + ' (×1,000)')
                if row == 3:
                    ax.set_xlabel('eGFP-A · asinh(x/1000)' if mode == 'reporters' else data['scatter_parameters'][0] + ' (×1,000)')
        fig.text(.075, .038, 'GFP % uses retained singlets; excluded % uses all primary singlets. D0 replicate columns are not paired to individual D4 wells.', fontsize=9, color='#52635b')
        fig.text(.075, .022, 'Detector settings differ between days. Shared axes show the recorded signals; they do not establish stable cross-day fluorescence calibration.', fontsize=8.5, color='#52635b')
        for ext in ('png', 'svg', 'pdf'):
            p = figdir / f'same_genotype_{mode}.{ext}'
            fig.savefig(p, dpi=180, facecolor='white')
            paths.append(p)
        plt.close(fig)
    template = Path(__file__).with_suffix('.html')
    html = template.read_text(encoding='utf-8').replace('__DATA__', json.dumps(data, separators=(',', ':'), allow_nan=False).replace('<', '\\u003c'))
    technical_path = OUTPUT / 'technical_checks.json'
    technical = json.loads(technical_path.read_text(encoding='utf-8')) if technical_path.exists() else None
    html = html.replace('__TECHNICAL__', json.dumps(technical, separators=(',', ':'), allow_nan=False).replace('<', '\\u003c'))
    (OUTPUT / 'report.html').write_text(html, encoding='utf-8')
    manifest = dict(generated_utc=datetime.now(timezone.utc).isoformat(), data_sha256=digest(path),
                    technical_checks_sha256=digest(technical_path) if technical else None,
                    sources_reverified=data['source_hashes'], matplotlib_version=matplotlib.__version__,
                    code={str(p.relative_to(REPO)): digest(p) for p in (Path(__file__), template)},
                    exports={str(p.relative_to(OUTPUT)): digest(p) for p in [*paths, OUTPUT/'report.html', OUTPUT/'wells.csv']})
    (OUTPUT / 'manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    print(str(OUTPUT / 'report.html'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('collect', 'render'))
    {'collect': collect, 'render': render}[parser.parse_args().action]()
