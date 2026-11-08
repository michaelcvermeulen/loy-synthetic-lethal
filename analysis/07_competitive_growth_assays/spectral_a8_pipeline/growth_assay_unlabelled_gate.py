"""Sensitivity analysis: a day-specific reporter-negative gate for unlabelled cells.

Unlabelled TADA2B cells form a GFP-negative cluster whose mCherry-channel position
shifts between D0 and D4, across the primary 6.5 cut. For each knockout background and
day, the mCherry cut is placed at the valley between that cluster and the mCherry-positive
peak of the pooled mixtures. A background keeps the primary cut on both days unless a
separating valley is found on both days, so lines without an unlabelled cluster (and dim
mCherry lines) are not clipped. The GFP cut and the frozen D0 classifier are unchanged.
Primary results are not modified; this writes a separate file that growth_assay_followup
merges as sensitivity columns.

    .\\codex\\flow_env\\python.exe -m spectral.a8_pipeline.growth_assay_unlabelled_gate
"""
from __future__ import annotations

import csv
from datetime import datetime, timezone
from hashlib import sha256
import json
import math
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
REPORT = REPO / 'analysis/spectral_growth_assay_D4_20261020'
OUTPUT = REPORT / 'unlabelled_gate_sensitivity.json'
LABELS = Path(__file__).with_name('growth_assay_D4_labels.json')
VALLEY = dict(gfp_below=4.0, bins=(4.0, 10.0, 0.1), smooth_bins=3, unlabelled_window=(5.0, 7.4), mcherry_window=(7.6, 9.0))


def digest(path: Path) -> str:
    return sha256(Path(path).read_bytes()).hexdigest()


def find_valley(mcherry: np.ndarray, settings: dict = VALLEY) -> dict:
    """Valley between the unlabelled and mCherry-positive peaks of a 1-D mCherry distribution."""
    start, stop, width = settings['bins']
    counts, edges = np.histogram(np.asarray(mcherry, float), bins=np.arange(start, stop + width / 2, width))
    smooth = np.convolve(counts, np.ones(settings['smooth_bins']) / settings['smooth_bins'], 'same')
    mid = (edges[:-1] + edges[1:]) / 2

    def peak(window):
        sel = (mid >= window[0]) & (mid <= window[1])
        return int(np.flatnonzero(sel)[np.argmax(smooth[sel])])

    low, high = peak(settings['unlabelled_window']), peak(settings['mcherry_window'])
    valley = low + int(np.argmin(smooth[low:high + 1]))
    if not low < valley < high or smooth[valley] >= min(smooth[low], smooth[high]):
        raise ValueError('No separating valley between the unlabelled and mCherry peaks')
    return dict(unlabelled_peak=round(float(mid[low]), 4), mcherry_peak=round(float(mid[high]), 4),
                valley=round(float(mid[valley]), 4), events=int(len(mcherry)),
                valley_depth=round(float(smooth[valley] / min(smooth[low], smooth[high])), 4))


def gated_counts(features: np.ndarray, weights, intercept: float, gfp_below: float, mcherry_below: float) -> dict:
    values = np.asarray(features, float)
    excluded = (values[:, 0] < gfp_below) & (values[:, 1] < mcherry_below)
    green = (values @ np.asarray(weights) + intercept > 0) & ~excluded
    retained = int((~excluded).sum())
    return dict(primary_singlets=len(values), retained=retained, gfp_count=int(green.sum()),
                excluded_fraction=float(excluded.mean()) if len(values) else None)


def log_odds(p):
    return math.log(p / (1 - p)) if p is not None and 0 < p < 1 else None


def collect() -> dict:
    from .competition_report import load_cef_gates
    from .fcs_io import load_fcs
    from .gating import points_in_polygon
    from .growth_assay_followup import GATE_NAMES, measure_reporters

    d4path = REPORT / 'results.json'
    d4 = json.loads(d4path.read_text(encoding='utf-8'))
    d0path = REPO / d4['baseline_path']
    d0 = json.loads(d0path.read_text(encoding='utf-8'))
    if digest(d0path) != d4['baseline_sha256'] or d4['identity_status'] != 'user_confirmed':
        raise ValueError('Expected the verified D0 baseline and approved D4 labels')
    if d4['identity_mapping_sha256'] != digest(LABELS):
        raise ValueError('D4 report was not built from the current label file')
    weights, intercept = d0['classifier']['weights'], d0['classifier']['intercept']
    gfp_below = d0['negative_rules']['broad']['gfp_below']
    features, source_hashes = {}, {}
    for day, data in ((0, d0), (4, d4)):
        expected = {m['path']: m['sha256'] for m in data['source_manifest']}
        source = REPO / data['source']
        gates = load_cef_gates(next(source.glob('*.cef')))
        raw = {p.stem: p for p in source.rglob('*.fcs')}
        for w in data['wells']:
            if w['role'] not in ('mixture', 'dye_neutral'):
                continue
            path = raw[w['sample'] if day == 0 else w['original_sample']]
            rel = path.relative_to(REPO).as_posix()
            if digest(path) != expected[rel]:
                raise ValueError(f'Input hash differs from approved analysis: {rel}')
            source_hashes[rel] = expected[rel]
            events = load_fcs(path).events
            keep = np.ones(len(events), dtype=bool)
            for name in GATE_NAMES:
                gate = gates[name]
                keep &= points_in_polygon(events[list(gate.parameters)].to_numpy(float), gate.vertices)
            f = np.arcsinh(events.loc[keep, ['eGFP-A', 'mCherry-A']].to_numpy(float) / 1000)
            primary = measure_reporters(f, d0)
            if (primary['reporter_positive_singlets'], primary['gfp_count']) != (w['reporter_positive_singlets'], w['gfp_count']):
                raise ValueError(f'Primary counts not reproduced: {path.name}')
            features[(day, w['well'])] = f
    # Background- and day-specific cuts from pooled mixtures; primary cut where no valley exists on both days.
    primary_cut = d0['negative_rules']['broad']['mcherry_below']
    background_of = {(day, w['well']): w['background'] for day, data in ((0, d0), (4, d4))
                     for w in data['wells'] if w['role'] in ('mixture', 'dye_neutral')}
    valleys, cut = {}, {}
    for background in sorted(set(background_of.values())):
        found = {}
        for day, data in ((0, d0), (4, d4)):
            wells = [w['well'] for w in data['wells'] if w['role'] == 'mixture' and w['background'] == background]
            pooled = np.vstack([features[(day, well)] for well in wells])
            try:
                found[day] = dict(find_valley(pooled[pooled[:, 0] < VALLEY['gfp_below'], 1]), wells=wells)
            except ValueError:
                found[day] = None
        detected = all(found.values())
        valleys[background] = dict(detected=detected, **{f'D{day}': v for day, v in found.items()})
        for day in (0, 4):
            cut[(day, background)] = found[day]['valley'] if detected else primary_cut

    counts = {key: gated_counts(f, weights, intercept, gfp_below, cut[(key[0], background_of[key])])
              for key, f in features.items()}
    d0_by_well = {w['well']: w for w in d0['wells']}
    baselines = {}
    for b in d0['d0_baselines']:
        wells = b['wells'].split(';')
        retained = sum(counts[(0, x)]['retained'] for x in wells)
        green = sum(counts[(0, x)]['gfp_count'] for x in wells)
        loy = green if b['orientation'] == 'LOY_GFP' else retained - green
        if any(d0_by_well[x]['orientation'] != b['orientation'] for x in wells):
            raise ValueError('D0 baseline wells disagree with their orientation')
        baselines[(b['group'], b['orientation'])] = dict(d0_loy_fraction=loy / retained, retained=retained, loy_count=loy)
    well_rows = []
    for w in d4['wells']:
        c = counts[(4, w['well'])]
        gfp = c['gfp_count'] / c['retained'] if c['retained'] else None
        row = dict(well=w['well'], original_sample=w['original_sample'], role=w['role'], group=w['group'],
                   orientation=w['orientation'], genotype=w['genotype'], replicate=w['replicate'],
                   primary_singlets=c['primary_singlets'], retained=c['retained'], gfp_count=c['gfp_count'],
                   excluded_fraction=c['excluded_fraction'], gfp_fraction=gfp)
        if w['role'] == 'mixture':
            b = baselines[(w['group'], w['orientation'])]
            loy = gfp if w['orientation'] == 'LOY_GFP' else (None if gfp is None else 1 - gfp)
            d0lo = log_odds(b['d0_loy_fraction'])
            row.update(loy_fraction=loy, d0_loy_fraction=b['d0_loy_fraction'],
                       delta_log_odds=None if log_odds(loy) is None or d0lo is None else log_odds(loy) - d0lo)
        well_rows.append(row)
    summaries, agreement = [], []
    groups = list(dict.fromkeys(r['group'] for r in well_rows if r['role'] == 'mixture'))
    for group in groups:
        means = {}
        for orientation in ('LOY_GFP', 'LOY_mCherry'):
            v = [r['delta_log_odds'] for r in well_rows if r['role'] == 'mixture' and r['group'] == group
                 and r['orientation'] == orientation and r['delta_log_odds'] is not None]
            means[orientation] = float(np.mean(v)) if v else None
            summaries.append(dict(group=group, orientation=orientation, n=len(v), delta_log_odds_mean=means[orientation],
                                  delta_log_odds_sd=float(np.std(v, ddof=1)) if len(v) > 1 else None,
                                  d0_loy_fraction=baselines[(group, orientation)]['d0_loy_fraction']))
        g, m = means['LOY_GFP'], means['LOY_mCherry']
        direction = ('unavailable' if None in (g, m) else 'both increase' if g > 0 and m > 0
                     else 'both decrease' if g < 0 and m < 0 else 'opposite or zero directions')
        agreement.append(dict(group=group, gfp_delta_log_odds=g, mcherry_delta_log_odds=m, direction_agreement=direction,
                              orientation_gap_log_odds=None if None in (g, m) else g - m,
                              mean_of_orientations=None if None in (g, m) else (g + m) / 2))
    controls = []
    for genotype in ('CTRL', 'LOY'):
        d0w = [w['well'] for w in d0['wells'] if w['role'] == 'dye_neutral' and w['genotype'] == genotype]
        p0 = sum(counts[(0, x)]['gfp_count'] for x in d0w) / sum(counts[(0, x)]['retained'] for x in d0w)
        d4v = [r['gfp_fraction'] for r in well_rows if r['role'] == 'dye_neutral' and r['genotype'] == genotype]
        controls.append(dict(genotype=genotype, d0_gfp_fraction_pooled=p0, d4_gfp_fraction_mean=float(np.mean(d4v)),
                             mean_delta_log_odds=float(np.mean([log_odds(v) - log_odds(p0) for v in d4v]))))
    for rel, expected in source_hashes.items():
        if digest(REPO / rel) != expected:
            raise ValueError(f'Source changed during analysis: {rel}')
    return dict(schema_version=1, created_utc=datetime.now(timezone.utc).isoformat(),
                purpose='Sensitivity only; primary broad-negative results are unchanged.',
                rule=dict(gfp_below=gfp_below, primary_mcherry_below=primary_cut, valley_settings=VALLEY,
                          mcherry_below={bg: {f'D{day}': cut[(day, bg)] for day in (0, 4)} for bg in valleys},
                          adjusted_backgrounds=[bg for bg, v in valleys.items() if v['detected']],
                          description='Exclude eGFP < 6 and mCherry below the background- and day-specific valley '
                                      'between the unlabelled and mCherry-positive peaks of GFP-low mixture events; '
                                      'backgrounds without a valley on both days keep the primary cut.'),
                valleys=valleys,
                baseline_sha256=d4['baseline_sha256'], identity_mapping_sha256=digest(LABELS),
                source_hashes=source_hashes, script_sha256=digest(Path(__file__)),
                wells=well_rows, orientation_summaries=summaries, dye_agreement=agreement, dye_controls=controls,
                d0_baselines=[dict(group=g, orientation=o, **b) for (g, o), b in baselines.items()])


def main() -> None:
    result = collect()
    OUTPUT.write_text(json.dumps(result, indent=1, allow_nan=False), encoding='utf-8')
    with OUTPUT.with_suffix('.csv').open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(result['dye_agreement'][0]))
        writer.writeheader()
        writer.writerows(result['dye_agreement'])
    print(json.dumps(dict(rule=result['rule']['mcherry_below'], valleys=result['valleys'],
                          agreement=result['dye_agreement'], controls=result['dye_controls']), indent=1))


if __name__ == '__main__':
    main()
