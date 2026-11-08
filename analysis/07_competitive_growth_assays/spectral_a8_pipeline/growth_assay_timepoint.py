r"""Later growth-assay timepoints with technical replicates, measured with the frozen D0 method.

Wells are measured exactly as at D4: the run's own CEF gates, the frozen D0 LDA, the broad
reporter-negative rule, and the day-specific unlabelled-cell gate in the approved backgrounds.
Each biological replicate may be acquired in several adjacent wells. Technical replicates are
averaged within their biological replicate before an arm is summarised, so spread is across
biological replicates. Raw acquisitions are read-only, sources are hashed before and after,
outputs are staged, and existing output directories are refused.

Collect with the flow environment, then render figures with a Matplotlib environment::

    .\codex\flow_env\python.exe -m spectral.a8_pipeline.growth_assay_timepoint --output <new-directory>
    python -m spectral.a8_pipeline.growth_assay_timepoint_plots --report-dir <directory>
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import re
import time

import numpy as np

from .competition_report import LDAModel, fit_lda, load_cef_gates
from .fcs_io import load_fcs
from .gating import points_in_polygon
from .growth_assay_d0_report import BACKGROUNDS
from .growth_assay_followup import GATE_NAMES, analyze, digest, inventory, measure_reporters, write_csv
from .growth_assay_label_review import GROUPS, ORIENTATIONS, direction, stats, strict_baseline
from .growth_assay_unlabelled_gate import VALLEY, find_valley, gated_counts

REPO = Path(__file__).resolve().parents[2]
BASELINE = REPO / 'analysis/spectral_growth_assay_D0_20261015/results.json'
DEFAULT_PREVIOUS = REPO / 'analysis/spectral_growth_assay_D4_20261020'
DEFAULT_INPUT = REPO / '20261029_160937_Export' / 'comp_growth_assay_D14_OCT20 1'
DEFAULT_LABELS = Path(__file__).with_name('growth_assay_D14_labels.json')
DEFAULT_GATING = Path(__file__).with_name('growth_assay_D14_gating.json')
DEFAULT_OUTPUT = REPO / 'analysis/spectral_growth_assay_D14_20261029'
GATES = ('standard', 'unlabelled_cell_gate')
# The C1-labelled arms are second, untreated replicate sets of their background; no inhibitor was
# added (experimenter, 2026-10-29). WT C1 was a separate D0 mix; DDX3X C1 was dispensed from the DDX3X mix.
SECOND_SETS = {'WT_C1': 'WT', 'DDX3X_C1': 'DDX3X'}
# Screening flags only; flagged wells stay in every primary summary. Reporter negatives are flagged when
# they are both common and at least twice their D0 mix, since TADA2B carries unlabelled cells from D0.
FLAGS = dict(excluded_fraction=0.20, excluded_fold_vs_d0=2.0, singlet_yield=0.60, total_events=10000,
             technical_gap_pp=5.0)
SUMMARY_KEYS = ('loy_fraction', 'delta_log_odds', 'delta_log_odds_standard_gate', 'delta_log_odds_unlabelled_gate',
                'delta_log_odds_strict', 'excluded_broad_fraction', 'delta_log_odds_excluded_as_ctrl',
                'delta_log_odds_excluded_as_loy')


def logit(p: float | None) -> float | None:
    return math.log(p / (1 - p)) if p is not None and 0 < p < 1 else None


def change(p: float | None, p0: float | None) -> float | None:
    """Change in log-odds; no pseudocount, so a fraction of 0 or 1 gives no value."""
    a, b = logit(p), logit(p0)
    return None if a is None or b is None else a - b


def mean(values) -> float | None:
    valid = [v for v in values if v is not None]
    return float(np.mean(valid)) if valid else None


def loy_of(gfp: float | None, orientation: str) -> float | None:
    return gfp if gfp is None or orientation == 'LOY_GFP' else 1 - gfp


def well_order(well: str) -> tuple[str, int]:
    return well[0], int(well[1:])


def assign_identities(wells: list[dict], labels: dict) -> None:
    """Recorded labels plus approved overrides; technical replicates follow plate column order."""
    assignments = labels.get('assignments', [])
    overrides = {a['original_sample']: a for a in assignments}
    if len(overrides) != len(assignments):
        raise ValueError('Duplicate identity mapping entry')
    if not set(overrides) <= {w['original_sample'] for w in wells}:
        raise ValueError('Identity mapping contains an unknown recording')
    for w in wells:
        for field in ('role', 'group', 'background', 'condition', 'orientation', 'genotype', 'reporter', 'replicate'):
            w[field] = w['nominal_' + field]
        a = overrides.get(w['original_sample'], {})
        if 'genotype' in a and (w['role'] != 'dye_neutral' or a['genotype'] not in ('CTRL', 'LOY')):
            raise ValueError('Genotype corrections must name CTRL or LOY for same-genotype controls')
        if a.get('replicate', 1) not in (1, 2, 3):
            raise ValueError('Expected replicate 1, 2 or 3')
        w['replicate'] = a.get('replicate', w['replicate'])
        w['genotype'] = a.get('genotype', w['genotype'])
        w['replicate_changed'] = w['replicate'] != w['nominal_replicate']
        w['genotype_changed'] = w['genotype'] != w['nominal_genotype']
        w['identity_status'] = 'corrected' if w['replicate_changed'] or w['genotype_changed'] else 'recorded'
    blocks = defaultdict(list)
    for w in wells:
        blocks[(w['role'], w['group'], w['orientation'], w['genotype'], w['replicate'])].append(w)
    for key, members in blocks.items():
        if len(members) != labels['technical_replicates']:
            raise ValueError(f"Expected {labels['technical_replicates']} technical replicates for {key}; found {len(members)}")
        for technical, w in enumerate(sorted(members, key=lambda w: well_order(w['well'])), 1):
            stem = re.sub(r'_[1-3](?:_\d{3})?$', f"_{w['replicate']}", w['original_sample'])
            if w['genotype_changed']:
                stem = stem.replace(f"-{w['nominal_genotype']}_GFP_CHERRY_", f"-{w['genotype']}_GFP_CHERRY_", 1)
            w['technical_replicate'] = technical
            w['analysis_sample'] = f'{stem}_t{technical}'


def check_design(wells: list[dict], baseline: dict) -> None:
    mixtures = {(w['group'], w['orientation'], w['replicate']) for w in wells if w['role'] == 'mixture'}
    neutral = {(w['genotype'], w['replicate']) for w in wells if w['role'] == 'dye_neutral'}
    if mixtures != {(b['group'], b['orientation'], r) for b in baseline['d0_baselines'] for r in (1, 2, 3)}:
        raise ValueError('Competition design has missing or extra biological replicates')
    if neutral != {(g, r) for g in ('CTRL', 'LOY') for r in (1, 2, 3)}:
        raise ValueError('Same-genotype design has missing or extra biological replicates')


def biological(members: list[dict], keys: tuple[str, ...], readout: str) -> list[dict]:
    """Average technical replicates (equal weight) within each biological replicate."""
    reps = defaultdict(list)
    for w in members:
        reps[w['replicate']].append(w)
    rows = []
    for replicate, ws in sorted(reps.items()):
        ws = sorted(ws, key=lambda w: w['technical_replicate'])
        values = [w[readout] for w in ws if w[readout] is not None]
        rows.append(dict(replicate=replicate, wells=';'.join(w['well'] for w in ws), technical_n=len(ws),
                         technical_gap_pp=100 * (max(values) - min(values)) if len(values) > 1 else None,
                         **{key: mean(w[key] for w in ws) for key in keys}))
    return rows


def largest_gap(rows: list[dict]) -> float | None:
    return max((r['technical_gap_pp'] for r in rows if r['technical_gap_pp'] is not None), default=None)


def valley_cuts(pooled: dict[str, dict[str, np.ndarray]], primary_cut: float) -> tuple[dict, dict]:
    """Per background, the unlabelled/mCherry valley on each day; the primary cut unless found on every day."""
    valleys, cuts = {}, {}
    for background, days in pooled.items():
        found = {}
        for day, features in days.items():
            try:
                found[day] = find_valley(features[features[:, 0] < VALLEY['gfp_below'], 1])
            except ValueError:
                found[day] = None
        detected = all(found.values())
        valleys[background] = dict(detected=detected, **found)
        cuts[background] = {day: found[day]['valley'] if detected else primary_cut for day in found}
    return valleys, cuts


def counts(features: np.ndarray, baseline: dict, mcherry_below: float) -> tuple[int, int, int]:
    c = gated_counts(features, baseline['classifier']['weights'], baseline['classifier']['intercept'],
                     baseline['negative_rules']['broad']['gfp_below'], mcherry_below)
    return c['primary_singlets'], c['retained'], c['gfp_count']


def reference(b: dict, well_counts: dict) -> dict:
    """Pool a D0 mixing batch: LOY among retained cells, and with every excluded event called CTRL or LOY."""
    n, retained, green = (sum(well_counts[x][i] for x in b['wells'].split(';')) for i in range(3))
    loy = green if b['orientation'] == 'LOY_GFP' else retained - green
    return dict(loy_fraction=loy / retained, excluded_fraction=(n - retained) / n, excluded_as_ctrl=loy / n,
                excluded_as_loy=(loy + n - retained) / n)


def gate_events(events, gates) -> np.ndarray:
    keep = np.ones(len(events), dtype=bool)
    for name in GATE_NAMES:
        gate = gates[name]
        keep &= points_in_polygon(events[list(gate.parameters)].to_numpy(float), gate.vertices)
    return keep


def load_d0(baseline: dict) -> tuple[dict, dict]:
    """Hash-verified D0 reporters of mixtures and same-genotype controls, checked against the D0 report."""
    expected = {m['path']: m['sha256'] for m in baseline['source_manifest']}
    source = REPO / baseline['source']
    gates = load_cef_gates(next(source.glob('*.cef')))
    raw = {p.stem: p for p in source.rglob('*.fcs')}
    features, hashes = {}, {}
    for w in baseline['wells']:
        if w['role'] not in ('mixture', 'dye_neutral'):
            continue
        path = raw[w['sample']]
        rel = path.relative_to(REPO).as_posix()
        if digest(path) != expected[rel]:
            raise ValueError(f'D0 input differs from the approved baseline: {rel}')
        events = load_fcs(path).events
        f = np.arcsinh(events.loc[gate_events(events, gates), ['eGFP-A', 'mCherry-A']].to_numpy(float) / 1000)
        m = measure_reporters(f, baseline)
        if (m['primary_singlets'], m['reporter_positive_singlets'], m['gfp_count']) != (
                w['primary_singlets'], w['reporter_positive_singlets'], w['gfp_count']):
            raise ValueError(f'D0 counts not reproduced: {path.name}')
        features[w['well']], hashes[rel] = f, expected[rel]
    return features, hashes


def measure(result: dict, baseline: dict, d0: dict, cells: dict, gating: dict) -> dict:
    """Primary, standard-gate, unlabelled-gate, strict and bounding values for every well."""
    timepoint, wells = result['timepoint'], result['wells']
    standard_cut = baseline['negative_rules']['broad']['mcherry_below']
    d0_background = {w['well']: w['background'] for w in baseline['wells'] if w['well'] in d0}
    d0_mixtures = {w['well'] for w in baseline['wells'] if w['role'] == 'mixture'}
    pooled = {bg: {'D0': np.vstack([d0[x] for x in sorted(d0_mixtures) if d0_background[x] == bg]),
                   timepoint: np.vstack([cells[w['original_sample']]['reporters'] for w in wells
                                         if w['role'] == 'mixture' and w['background'] == bg])}
              for bg in BACKGROUNDS}
    valleys, cuts = valley_cuts(pooled, standard_cut)
    if gating['source'] != result['source']:
        raise ValueError('Gating decision belongs to a different export')
    primary = set(gating['unlabelled_gate_primary_backgrounds'])
    if not primary <= {bg for bg, v in valleys.items() if v['detected']}:
        raise ValueError('The unlabelled-cell gate is primary only where a separating valley is found on both days')

    def gate_of(background):
        return 'unlabelled_cell_gate' if background in primary else 'standard'

    d0_counts = {gate: {well: counts(f, baseline, cuts[d0_background[well]]['D0'] if gate == GATES[1] else standard_cut)
                        for well, f in d0.items()} for gate in GATES}
    references = {}
    for b in baseline['d0_baselines']:
        refs = {gate: reference(b, d0_counts[gate]) for gate in GATES}
        if not math.isclose(refs['standard']['loy_fraction'], b['d0_loy_fraction'], rel_tol=0, abs_tol=1e-12):
            raise ValueError('Standard-gate D0 reference does not reproduce the D0 report')
        references[(b['group'], b['orientation'])] = dict(
            batch=b['mixture_batch'], wells=b['wells'], strict=strict_baseline(baseline, b),
            primary_gate=gate_of(b['group'].split('_')[0]), **refs)
    neutral_d0 = {}
    for genotype in ('CTRL', 'LOY'):
        members = [w['well'] for w in baseline['wells'] if w['role'] == 'dye_neutral' and w['genotype'] == genotype]
        n, retained, green = (sum(d0_counts['standard'][x][i] for x in members) for i in range(3))
        neutral_d0[genotype] = dict(wells=';'.join(members), gfp_fraction=green / retained,
                                    excluded_fraction=(n - retained) / n)
    for w in wells:
        features = cells[w['original_sample']]['reporters']
        if counts(features, baseline, standard_cut) != (w['primary_singlets'], w['reporter_positive_singlets'], w['gfp_count']):
            raise ValueError(f"Standard-gate counts not reproduced: {w['well']}")
        n, retained, green = counts(features, baseline, cuts[w['background']][timepoint])
        gfp_unlabelled = green / retained if retained else None
        w['primary_gate'] = gate_of(w['background'])
        w['excluded_fraction_unlabelled_gate'] = (n - retained) / n if n else None
        w['gfp_fraction_unlabelled_gate'] = gfp_unlabelled
        for key in ('excluded_broad_fraction', 'gfp_fraction'):
            w[key + '_standard_gate'] = w[key]
        if w['primary_gate'] == 'unlabelled_cell_gate':
            w.update(reporter_positive_singlets=retained, excluded_broad_count=n - retained,
                     excluded_broad_fraction=(n - retained) / n if n else None, gfp_count=green,
                     mcherry_count=retained - green, gfp_fraction=gfp_unlabelled,
                     mcherry_fraction=None if gfp_unlabelled is None else 1 - gfp_unlabelled)
        if w['role'] == 'dye_neutral':
            ref = neutral_d0[w['genotype']]
            w.update(d0_gfp_fraction=ref['gfp_fraction'], d0_wells=ref['wells'], d0_excluded_fraction=ref['excluded_fraction'],
                     gfp_delta_log_odds=change(w['gfp_fraction'], ref['gfp_fraction']))
            continue
        ref, o = references[(w['group'], w['orientation'])], w['orientation']
        w.update(d0_batch=ref['batch'], d0_wells=ref['wells'], d0_excluded_fraction=ref[w['primary_gate']]['excluded_fraction'])
        for suffix, gfp, gate in (('_standard_gate', w['gfp_fraction_standard_gate'], 'standard'),
                                  ('_unlabelled_gate', gfp_unlabelled, 'unlabelled_cell_gate')):
            w['loy_fraction' + suffix] = loy_of(gfp, o)
            w['d0_loy_fraction' + suffix] = ref[gate]['loy_fraction']
            w['delta_log_odds' + suffix] = change(loy_of(gfp, o), ref[gate]['loy_fraction'])
        chosen = ref[w['primary_gate']]
        w['loy_fraction'], w['d0_loy_fraction'] = loy_of(w['gfp_fraction'], o), chosen['loy_fraction']
        w['delta_log_odds'] = change(w['loy_fraction'], w['d0_loy_fraction'])
        w['loy_count'] = w['gfp_count'] if o == 'LOY_GFP' else w['reporter_positive_singlets'] - w['gfp_count']
        n, excluded = w['primary_singlets'], w['primary_singlets'] - w['reporter_positive_singlets']
        w['loy_fraction_excluded_as_ctrl'] = w['loy_count'] / n if n else None
        w['loy_fraction_excluded_as_loy'] = (w['loy_count'] + excluded) / n if n else None
        w['delta_log_odds_excluded_as_ctrl'] = change(w['loy_fraction_excluded_as_ctrl'], chosen['excluded_as_ctrl'])
        w['delta_log_odds_excluded_as_loy'] = change(w['loy_fraction_excluded_as_loy'], chosen['excluded_as_loy'])
        w['loy_fraction_strict'] = loy_of(w['gfp_fraction_strict'], o)
        w['d0_loy_fraction_strict'] = ref['strict']
        w['delta_log_odds_strict'] = change(w['loy_fraction_strict'], ref['strict'])
    return dict(valleys=valleys, cuts=cuts, primary_backgrounds=sorted(primary),
                references=[dict(group=g, orientation=o, **{k: v for k, v in r.items()}) for (g, o), r in references.items()],
                same_genotype_references=neutral_d0)


def summarise(wells: list[dict], previous: dict) -> dict:
    """Arm summaries across biological replicates, dye agreement, C1 contrasts and same-genotype controls."""
    before = {(s['group'], s['orientation']): s for s in previous['orientation_summaries']}
    summaries, replicates = [], []
    for group in GROUPS:
        for orientation in ORIENTATIONS:
            members = [w for w in wells if w['role'] == 'mixture' and w['group'] == group and w['orientation'] == orientation]
            bio = biological(members, SUMMARY_KEYS, 'loy_fraction')
            replicates += [dict(group=group, orientation=orientation, **b) for b in bio]
            s = dict(group=group, orientation=orientation, primary_gate=members[0]['primary_gate'],
                     d0_loy_fraction=members[0]['d0_loy_fraction'], d0_batch=members[0]['d0_batch'],
                     d0_wells=members[0]['d0_wells'],
                     wells_by_replicate=' | '.join(f"{b['replicate']}: {b['wells']}" for b in bio))
            for key in SUMMARY_KEYS:
                s.update({f'{key}_{k}': v for k, v in stats([b[key] for b in bio]).items()})
            s['technical_gap_pp_max'] = largest_gap(bio)
            p = before[(group, orientation)]
            s.update(previous_loy_fraction_mean=p['loy_fraction_mean'], previous_delta_log_odds_mean=p['delta_log_odds_mean'],
                     previous_delta_log_odds_standard_gate_mean=p['delta_log_odds_standard_gate_mean'],
                     change_since_previous_log_odds=s['delta_log_odds_mean'] - p['delta_log_odds_mean'])
            summaries.append(s)
    by = {(s['group'], s['orientation']): s for s in summaries}
    agreement = []
    for group in GROUPS:
        g, m = by[(group, 'LOY_GFP')], by[(group, 'LOY_mCherry')]
        values = [g['delta_log_odds_mean'], m['delta_log_odds_mean']]
        standard = [g['delta_log_odds_standard_gate_mean'], m['delta_log_odds_standard_gate_mean']]
        earlier = [g['previous_delta_log_odds_mean'], m['previous_delta_log_odds_mean']]
        agreement.append(dict(group=group, primary_gate=g['primary_gate'], gfp_delta_log_odds=values[0],
                              mcherry_delta_log_odds=values[1], direction_agreement=direction(values),
                              orientation_gap_log_odds=values[0] - values[1], mean_of_orientations=mean(values),
                              gfp_delta_log_odds_standard_gate=standard[0], mcherry_delta_log_odds_standard_gate=standard[1],
                              standard_gate_direction_agreement=direction(standard),
                              mean_of_orientations_standard_gate=mean(standard),
                              previous_mean_of_orientations=mean(earlier),
                              change_since_previous_mean_of_orientations=mean(values) - mean(earlier)))
    comparisons = []
    for second, first in SECOND_SETS.items():
        for orientation in (*ORIENTATIONS, 'mean_of_orientations'):
            def value(group, key):
                if orientation == 'mean_of_orientations':
                    return mean(by[(group, o)][key] for o in ORIENTATIONS)
                return by[(group, orientation)][key]
            comparisons.append(dict(
                background=first, first_set=first, second_set=second, orientation=orientation,
                shared_d0_mix=by[(first, 'LOY_GFP')]['d0_batch'] == by[(second, 'LOY_GFP')]['d0_batch'],
                first_delta_log_odds=value(first, 'delta_log_odds_mean'),
                second_delta_log_odds=value(second, 'delta_log_odds_mean'),
                difference_since_d0=value(second, 'delta_log_odds_mean') - value(first, 'delta_log_odds_mean'),
                first_change_since_previous=value(first, 'change_since_previous_log_odds'),
                second_change_since_previous=value(second, 'change_since_previous_log_odds'),
                difference_since_previous=value(second, 'change_since_previous_log_odds')
                - value(first, 'change_since_previous_log_odds'),
                difference_since_d0_standard_gate=value(second, 'delta_log_odds_standard_gate_mean')
                - value(first, 'delta_log_odds_standard_gate_mean')))
    # Each background over all of its replicate sets: biological replicates pooled, and the mean of set means.
    backgrounds = []
    for background in BACKGROUNDS:
        sets = [g for g in GROUPS if g.split('_')[0] == background]
        for orientation in ORIENTATIONS:
            bio = [r for r in replicates if r['group'] in sets and r['orientation'] == orientation]
            row = dict(background=background, orientation=orientation, sets=';'.join(sets),
                       d0_mixes=';'.join(dict.fromkeys(by[(g, orientation)]['d0_batch'] for g in sets)),
                       set_means='; '.join(f"{g}: {by[(g, orientation)]['delta_log_odds_mean']:+.3f}" for g in sets))
            row.update({f'delta_log_odds_{k}': v for k, v in stats([r['delta_log_odds'] for r in bio]).items()})
            row.update(mean_of_set_means=mean(by[(g, orientation)]['delta_log_odds_mean'] for g in sets),
                       previous_mean_of_set_means=mean(by[(g, orientation)]['previous_delta_log_odds_mean'] for g in sets))
            backgrounds.append(row)
    earlier = {c['genotype']: c for c in previous['dye_controls']}
    tp = previous['timepoint'].lower()
    controls = []
    for genotype in ('CTRL', 'LOY'):
        members = [w for w in wells if w['role'] == 'dye_neutral' and w['genotype'] == genotype]
        bio = biological(members, ('gfp_fraction', 'gfp_delta_log_odds'), 'gfp_fraction')
        replicates += [dict(group='same_genotype', orientation=genotype, **b) for b in bio]
        p = earlier[genotype]
        prior_gfp = p.get('gfp_fraction_mean', p.get(f'{tp}_gfp_fraction_mean'))
        s = dict(genotype=genotype, d0_gfp_fraction_pooled=members[0]['d0_gfp_fraction'], d0_wells=members[0]['d0_wells'],
                 wells_by_replicate=' | '.join(f"{b['replicate']}: {b['wells']}" for b in bio),
                 previous_gfp_fraction_mean=prior_gfp,
                 previous_delta_log_odds=p.get('gfp_delta_log_odds_mean', p.get('mean_delta_log_odds')))
        for key in ('gfp_fraction', 'gfp_delta_log_odds'):
            s.update({f'{key}_{k}': v for k, v in stats([b[key] for b in bio]).items()})
        s['technical_gap_pp_max'] = largest_gap(bio)
        s['change_since_previous_log_odds'] = s['gfp_delta_log_odds_mean'] - s['previous_delta_log_odds']
        controls.append(s)
    return dict(orientation_summaries=summaries, biological_replicates=replicates, dye_agreement=agreement,
                replicate_set_comparisons=comparisons, background_summaries=backgrounds, dye_controls=controls)


def flag_wells(wells: list[dict], replicates: list[dict]) -> None:
    gaps = {well: r['technical_gap_pp'] for r in replicates for well in r['wells'].split(';')}
    for w in wells:
        flags, excluded, before = [], w['excluded_broad_fraction'], w['d0_excluded_fraction']
        if excluded is not None and excluded > FLAGS['excluded_fraction'] and excluded > FLAGS['excluded_fold_vs_d0'] * before:
            flags.append(f'reporter-negative {100 * excluded:.0f}% (D0 mix {100 * before:.0f}%)')
        if w['singlet_yield'] is not None and w['singlet_yield'] < FLAGS['singlet_yield']:
            flags.append(f"singlet yield < {FLAGS['singlet_yield']:.0%}")
        if w['total_events'] < FLAGS['total_events']:
            flags.append(f"fewer than {FLAGS['total_events']:,} events")
        if (gaps.get(w['well']) or 0) > FLAGS['technical_gap_pp']:
            flags.append(f"technical replicates differ by > {FLAGS['technical_gap_pp']:g} pp")
        w['qc_flags'] = '; '.join(flags)


def identity_checks(wells: list[dict], cells: dict, baseline: dict) -> dict:
    """WT-background fingerprints from per-cell LightLoss (Violet)-A, SSC (Imaging)-A and reporter brightness."""
    model = LDAModel(np.array(baseline['classifier']['weights']), baseline['classifier']['intercept'])
    rule = baseline['negative_rules']['broad']

    def classes(w):  # standard rule: 0 excluded, 1 mCherry-like, 2 GFP-like
        f = cells[w['original_sample']]['reporters']
        excluded = (f[:, 0] < rule['gfp_below']) & (f[:, 1] < rule['mcherry_below'])
        return np.where(excluded, 0, np.where(model.predict(f) == 1, 2, 1))

    def features(w, cls):  # log scatter and own-reporter brightness of one reporter class
        c, f = classes(w), cells[w['original_sample']]
        own = f['reporters'][:, 0] if cls == 2 else f['reporters'][:, 1]
        return np.column_stack([np.log10(np.clip(f['scatter'][c == cls], 1, None)), own[c == cls]])

    size = []
    for w in wells:
        if w['background'] != 'WT':
            continue
        c, light = classes(w), cells[w['original_sample']]['scatter'][:, 0]
        g, m = float(np.median(light[c == 2])), float(np.median(light[c == 1]))
        row = dict(well=w['well'], analysis_sample=w['analysis_sample'], role=w['role'], group=w['group'],
                   orientation=w['orientation'], genotype=w['genotype'], gfp_class_lightloss_median=g,
                   mcherry_class_lightloss_median=m, gfp_to_mcherry_ratio=g / m)
        if w['role'] == 'mixture':
            row.update(expected='GFP class larger' if w['orientation'] == 'LOY_GFP' else 'mCherry class larger',
                       consistent=bool((g > m) == (w['orientation'] == 'LOY_GFP')))
        size.append(row)
    arms = {group: [w for w in wells if w['role'] == 'mixture' and w['group'] == group] for group in ('WT', 'WT_C1')}
    neutral = sorted((w for w in wells if w['role'] == 'dye_neutral'), key=lambda w: well_order(w['well']))
    estimates, calibration = defaultdict(dict), []
    for train, test in (('WT', 'WT_C1'), ('WT_C1', 'WT')):
        for cls, name in ((2, 'gfp'), (1, 'mcherry')):
            data = {}
            for arm in (train, test):
                X = [features(w, cls) for w in arms[arm]]
                loy = [cls == (2 if w['orientation'] == 'LOY_GFP' else 1) for w in arms[arm]]
                data[arm] = np.vstack(X), np.concatenate([np.full(len(x), int(l)) for x, l in zip(X, loy)])
            fitted = fit_lda(*data[train])
            called = fitted.predict(data[test][0]).astype(bool)
            tpr, fpr = float(called[data[test][1] == 1].mean()), float(called[data[test][1] == 0].mean())
            calibration.append(dict(train=train, test=test, reporter_class=name, held_out_tpr=tpr, held_out_fpr=fpr))
            for w in neutral:
                share = float(fitted.predict(features(w, cls)).mean())
                estimates[w['well']][f'loy_content_{train}_to_{test}_{name}'] = (share - fpr) / (tpr - fpr)
    genotype = []
    for w in neutral:
        e = estimates[w['well']]
        low, high = min(e.values()), max(e.values())
        call = 'LOY' if low > .5 else 'CTRL' if high < .5 else 'unclear'
        genotype.append(dict(well=w['well'], analysis_sample=w['analysis_sample'], recorded_genotype=w['nominal_genotype'],
                             analysis_genotype=w['genotype'], loy_content_minimum=low, loy_content_maximum=high,
                             call=call, consistent=call == w['genotype'], **e))
    return dict(size_fingerprint=size, same_genotype=genotype, calibration=calibration,
                description='Size: median LightLoss (Violet)-A of GFP-like and mCherry-like retained cells; LOY cells '
                            'are larger in WT backgrounds. Same-genotype: per reporter class, an equal-prior LDA on '
                            'log10 LightLoss (Violet)-A, log10 SSC (Imaging)-A and own-reporter brightness, trained on '
                            'the WT or WT C1 competition cells of known genotype on this plate and calibrated on the '
                            'other arm; LOY content = (share called LOY - FPR) / (TPR - FPR). The imaging size '
                            'parameters used at D4 changed scale at D14 and are not used.')


def gate_changes(current: dict, previous: dict) -> dict:
    """Each polygon relative to the previous export: unchanged, translated (with its shift) or reshaped."""
    out = {}
    for name, g in current.items():
        p = previous[name]
        if g['parameters'] != p['parameters'] or len(g['vertices']) != len(p['vertices']):
            out[name] = dict(status='reshaped')
            continue
        dx = np.array([a['x'] - b['x'] for a, b in zip(g['vertices'], p['vertices'])])
        dy = np.array([a['y'] - b['y'] for a, b in zip(g['vertices'], p['vertices'])])
        # Vertices are stored at float32 precision (spacing 8 near 1e8), so allow a few ulps of jitter.
        scale = max(abs(v[k]) for v in g['vertices'] + p['vertices'] for k in ('x', 'y'))
        if not dx.any() and not dy.any():
            out[name] = dict(status='unchanged')
        elif max(np.ptp(dx), np.ptp(dy)) <= 1 + 1e-6 * scale:
            out[name] = dict(status='translated', parameters=g['parameters'], shift=[float(dx.mean()), float(dy.mean())])
        else:
            out[name] = dict(status='reshaped')
    return out


def reference_gate_check(wells: list[dict], baseline: dict, cef: Path) -> dict:
    """Re-gate every recording with an earlier export's polygons; standard reporter rule throughout."""
    gates = load_cef_gates(cef)
    rows = []
    for w in wells:
        events = load_fcs(REPO / w['source_file']).events
        f = np.arcsinh(events.loc[gate_events(events, gates), ['eGFP-A', 'mCherry-A']].to_numpy(float) / 1000)
        m = measure_reporters(f, baseline)
        rows.append(dict(well=w['well'], primary_singlets=m['primary_singlets'], retained=m['reporter_positive_singlets'],
                         gfp_count=m['gfp_count'], gfp_fraction=m['gfp_fraction'],
                         gfp_fraction_difference=None if m['gfp_fraction'] is None else
                         m['gfp_fraction'] - w['gfp_fraction_standard_gate']))
    by_well = {w['well']: w for w in wells}
    return dict(cef=cef.relative_to(REPO).as_posix(), wells=rows,
                singlet_count_differences=sum(r['primary_singlets'] != by_well[r['well']]['primary_singlets'] for r in rows),
                max_abs_gfp_fraction_difference=max(abs(r['gfp_fraction_difference'] or 0) for r in rows))


def pct(v: float | None, digits: int = 1) -> str:
    return 'NA' if v is None else f'{100 * v:.{digits}f}'


def num(v: float | None, digits: int = 3) -> str:
    return 'NA' if v is None else f'{v:+.{digits}f}'


def plain(v: float | None, digits: int = 1) -> str:
    return 'NA' if v is None else f'{v:.{digits}f}'


def mean_sd(s: dict, key: str, scale: float = 1, digits: int = 3, signed: bool = True) -> str:
    m, sd = s[f'{key}_mean'], s[f'{key}_sd']
    if m is None:
        return 'NA'
    text = f'{scale * m:+.{digits}f}' if signed else f'{scale * m:.{digits}f}'
    return text + ('' if sd is None else f' ± {scale * sd:.{digits}f}')


def report_markdown(result: dict) -> list[str]:
    qc, tp, prev = result['qc'], result['timepoint'], result['previous_timepoint']
    meta = result['acquisition_metadata']
    gates = result['gate_changes']
    moved = {k: v for k, v in gates.items() if v['status'] != 'unchanged'}
    check = result.get('reference_gate_check')
    ids = result['identity_checks']
    size_mix = [r for r in ids['size_fingerprint'] if r['role'] == 'mixture']
    rule = result['unlabelled_gate']['rule']
    primary = result['unlabelled_gate']['primary_backgrounds']
    lines = [
        f"# {tp} competition assay — {'approved labels' if result['identity_status'] == 'user_confirmed' else 'recorded labels'}, "
        f'compared with D0 and {prev}', '',
        f"Source: `{result['source']}`; acquired {', '.join(meta[f'{tp.lower()}_fcs_dates'])}. "
        f"The D0 plate was acquired {meta['d0_fcs_date']}, {result['elapsed_days_from_d0']} days earlier; "
        f"the timepoint is labelled {tp} ({result['identity_mapping']['day_label']}).",
        f"Loaded {qc['recordings']} recordings / {qc['total_events']:,} events; {qc['primary_singlets']:,} primary singlets; "
        f"{qc['retained_events']:,} retained under the standard rule.", '',
        f"**Identity status: {result['identity_status']}.** {result['identity_mapping']['identity_basis']}", '',
        '## Method and checks', '',
        '- Run CEF Scatter → SSC Singlets → FSC Singlets polygons, verified against each recording ZIP; every Statistics '
        f"CSV gate count reproduced ({qc['gate_counts_matched']}/{qc['recordings']}).",
        '- Frozen D0 LDA (no retraining); arcsinh(value / 1000); broad negatives excluded (eGFP < 6 and mCherry < 6.5); '
        'GFP/mCherry fractions use retained singlets.',
        f"- Same unmixing payload and optical configuration as D0; {result['settings_summary']}",
        *([f"- Gate change since {prev}: " + '; '.join(
            f"{k} {v['status']}" + (f" by ({v['shift'][0]:,.0f}, {v['shift'][1]:,.0f}) on {' / '.join(v['parameters'])}"
                                    if v['status'] == 'translated' else '') for k, v in moved.items())
           + (f". Re-gating every recording with the {prev} polygons changes {check['singlet_count_differences']} singlet "
              f"counts; the largest GFP-fraction difference is {100 * check['max_abs_gfp_fraction_difference']:.3f} pp."
              if check else '.')] if moved else [f'- All gate polygons are identical to {prev}.']),
        f"- Technical replicates: {result['identity_mapping']['technical_replicates']} wells per biological replicate, "
        'numbered by plate column. They are averaged with equal weight within each biological replicate; arm means '
        'and SDs are across the three biological replicates.',
        '- All export file sizes and SHA-256 hashes were recorded and rechecked after analysis, as were the D0 sources.', '',
        '## Identity checks', '',
        f"- WT size fingerprint: {sum(r['consistent'] for r in size_mix)}/{len(size_mix)} WT and WT C1 competition wells "
        'have the recorded LOY reporter on the larger cells (GFP:mCherry LightLoss ratio '
        f"{min(r['gfp_to_mcherry_ratio'] for r in size_mix if r['orientation'] == 'LOY_GFP'):.2f}–"
        f"{max(r['gfp_to_mcherry_ratio'] for r in size_mix if r['orientation'] == 'LOY_GFP'):.2f} in LOY-GFP wells, "
        f"{min(r['gfp_to_mcherry_ratio'] for r in size_mix if r['orientation'] == 'LOY_mCherry'):.2f}–"
        f"{max(r['gfp_to_mcherry_ratio'] for r in size_mix if r['orientation'] == 'LOY_mCherry'):.2f} in LOY-mCherry wells).",
        *[f"- Same-genotype {g}/{g} wells ({', '.join(r['well'] for r in ids['same_genotype'] if r['analysis_genotype'] == g)}): "
          f"estimated LOY content {100 * min(r['loy_content_minimum'] for r in ids['same_genotype'] if r['analysis_genotype'] == g):.0f}"
          f"–{100 * max(r['loy_content_maximum'] for r in ids['same_genotype'] if r['analysis_genotype'] == g):.0f}%; "
          f"{sum(r['consistent'] for r in ids['same_genotype'] if r['analysis_genotype'] == g)}/"
          f"{sum(r['analysis_genotype'] == g for r in ids['same_genotype'])} called {g}." for g in ('LOY', 'CTRL')],
        '- Knockout backgrounds have no reliable LOY/CTRL size difference, so rows with knockout arms are not '
        'independently verified; their recorded labels are used.',
        f"- Method: {ids['description']}", '',
        '## Gating', '',
        '| Background | D0 valley | ' + f'{tp} valley | Valley on both days | Primary gate | D0 cut | {tp} cut |',
        '|---|---:|---:|---|---|---:|---:|',
        *[f"| {bg} | {v['D0']['valley'] if v['D0'] else '—'} | {v[tp]['valley'] if v[tp] else '—'} | "
          f"{'yes' if v['detected'] else 'no'} | {'unlabelled-cell' if bg in primary else 'standard'} | "
          f"{rule['mcherry_below'][bg]['D0'] if bg in primary else rule['primary_mcherry_below']} | "
          f"{rule['mcherry_below'][bg][tp] if bg in primary else rule['primary_mcherry_below']} |"
          for bg, v in result['unlabelled_gate']['valleys'].items()], '',
        f"Primary gating: {result['gating']['approval']} Standard-gate values are kept for every arm as a sensitivity.", '',
        '## Competition arms', '',
        f'Change in log-odds = ln(LOY/CTRL) at the timepoint minus ln(LOY/CTRL) in the pooled D0 mix. {tp} values are '
        'mean ± SD across biological replicates.', '',
        f'| Arm | LOY reporter | Gate | D0 LOY % | {prev} LOY % | {tp} LOY % | Δ log-odds D0→{prev} | '
        f'Δ log-odds D0→{tp} | Change {prev}→{tp} | Excluded % | Max technical gap (pp) |',
        '|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|',
        *[f"| {s['group']} | {'GFP' if s['orientation'] == 'LOY_GFP' else 'mCherry'} | "
          f"{'unlabelled' if s['primary_gate'] != 'standard' else 'standard'} | {pct(s['d0_loy_fraction'])} | "
          f"{pct(s['previous_loy_fraction_mean'])} | {mean_sd(s, 'loy_fraction', 100, 1, False)} | "
          f"{num(s['previous_delta_log_odds_mean'])} | {mean_sd(s, 'delta_log_odds')} | {num(s['change_since_previous_log_odds'])} | "
          f"{pct(s['excluded_broad_fraction_mean'])} | {plain(s['technical_gap_pp_max'])} |" for s in result['orientation_summaries']], '',
        '| Arm | Δ LOY-GFP | Δ LOY-mCherry | Direction | Orientation gap | Mean of orientations | '
        f'{prev} mean of orientations | Standard gate: LOY-GFP / LOY-mCherry |',
        '|---|---:|---:|---|---:|---:|---:|---:|',
        *[f"| {a['group']} | {num(a['gfp_delta_log_odds'])} | {num(a['mcherry_delta_log_odds'])} | {a['direction_agreement']} | "
          f"{num(a['orientation_gap_log_odds'])} | {num(a['mean_of_orientations'])} | {num(a['previous_mean_of_orientations'])} | "
          f"{num(a['gfp_delta_log_odds_standard_gate'])} / {num(a['mcherry_delta_log_odds_standard_gate'])} |"
          for a in result['dye_agreement']], '',
        '## Replicate sets', '',
        'The C1-labelled arms received no drug; they are second, untreated replicate sets of their background. WT C1 '
        'was a separate D0 mix, so WT and WT C1 are independent runs; DDX3X C1 was dispensed from the DDX3X mix. '
        'Differences are second set minus first set, in log-odds.', '',
        f'| Background | Orientation | Shared D0 mix | First set | Second set | Difference since D0 | '
        f'Difference since {prev} |',
        '|---|---|---|---:|---:|---:|---:|',
        *[f"| {c['background']} | {c['orientation'].replace('_', ' ')} | {'yes' if c['shared_d0_mix'] else 'no'} | "
          f"{num(c['first_delta_log_odds'])} | {num(c['second_delta_log_odds'])} | {num(c['difference_since_d0'])} | "
          f"{num(c['difference_since_previous'])} |" for c in result['replicate_set_comparisons']], '',
        'Each background over all of its replicate sets. The pooled SD is across biological replicates of every set; '
        'the mean of set means weights each set equally.', '',
        f'| Background | LOY reporter | D0 mixes | Set means | Pooled mean ± SD (n) | Mean of set means | '
        f'{prev} mean of set means |',
        '|---|---|---|---|---:|---:|---:|',
        *[f"| {b['background']} | {'GFP' if b['orientation'] == 'LOY_GFP' else 'mCherry'} | {b['d0_mixes']} | "
          f"{b['set_means']} | {mean_sd(b, 'delta_log_odds')} ({b['delta_log_odds_n']}) | {num(b['mean_of_set_means'])} | "
          f"{num(b['previous_mean_of_set_means'])} |" for b in result['background_summaries']], '',
        '## Same-genotype controls (WT background)', '',
        f'| Mixture | D0 GFP % | {prev} GFP % | {tp} GFP % | Δ GFP log-odds D0→{tp} | Change {prev}→{tp} |',
        '|---|---:|---:|---:|---:|---:|',
        *[f"| {c['genotype']}/{c['genotype']} | {pct(c['d0_gfp_fraction_pooled'])} | {pct(c['previous_gfp_fraction_mean'])} | "
          f"{mean_sd(c, 'gfp_fraction', 100, 1, False)} | {mean_sd(c, 'gfp_delta_log_odds')} | {num(c['change_since_previous_log_odds'])} |"
          for c in result['dye_controls']], '',
        'A GFP shift in LOY/LOY or CTRL/CTRL measures reporter-subline growth within one genotype. Such shifts move the '
        'two competition orientations in opposite directions and cancel in the mean of orientations if they are additive.', '',
        '## Reporter-negative bounds', '',
        'Arms whose mean excluded fraction exceeds 20% of singlets. Bounds assign every excluded event to CTRL or to LOY '
        'on both days (same gate).', '',
        '| Arm | LOY reporter | Excluded % | Δ log-odds (primary) | All excluded = CTRL | All excluded = LOY |',
        '|---|---|---:|---:|---:|---:|',
        *[f"| {s['group']} | {'GFP' if s['orientation'] == 'LOY_GFP' else 'mCherry'} | {pct(s['excluded_broad_fraction_mean'])} | "
          f"{num(s['delta_log_odds_mean'])} | {num(s['delta_log_odds_excluded_as_ctrl_mean'])} | {num(s['delta_log_odds_excluded_as_loy_mean'])} |"
          for s in result['orientation_summaries'] if (s['excluded_broad_fraction_mean'] or 0) > FLAGS['excluded_fraction']], '',
        '## Quality flags', '',
        f"Screening flags only; every flagged well stays in the primary summaries. Rules: reporter-negative fraction "
        f"above {FLAGS['excluded_fraction']:.0%} of singlets and more than {FLAGS['excluded_fold_vs_d0']:g}× its D0 mix "
        f"(same gate); singlet yield below {FLAGS['singlet_yield']:.0%}; fewer than {FLAGS['total_events']:,} events; "
        f"technical replicates more than {FLAGS['technical_gap_pp']:g} pp apart.", '',
        *[f"- {w['well']} ({w['analysis_sample']}): {w['qc_flags']}" for w in result['wells'] if w['qc_flags']], '',
        '## Interpretation limits', '',
        *['- ' + t for t in result['interpretation']], '',
        '## Reproduce', '',
        '```powershell',
        '.\\codex\\flow_env\\python.exe -m spectral.a8_pipeline.growth_assay_timepoint --output <new-directory>',
        'python -m spectral.a8_pipeline.growth_assay_timepoint_plots --report-dir <new-directory>',
        '```', '',
        '| Well | Analysis identity | Recorded filename | Tech. rep | LOY or GFP % | Δ log-odds | Excluded % | Flags |',
        '|---|---|---|---:|---:|---:|---:|---|',
    ]
    for w in result['wells']:
        mixture = w['role'] == 'mixture'
        value = w['loy_fraction'] if mixture else w['gfp_fraction']
        delta = w['delta_log_odds'] if mixture else w['gfp_delta_log_odds']
        lines.append(f"| {w['well']} | {w['analysis_sample']} | {w['original_filename']} | {w['technical_replicate']} | "
                     f"{'LOY ' if mixture else 'GFP '}{pct(value, 2)} | {num(delta)} | {pct(w['excluded_broad_fraction'], 2)} | "
                     f"{w['qc_flags'] or ''} |")
    return lines


def write_outputs(result: dict, output: Path, same_genotype_rows: list[dict]) -> None:
    output = output.resolve()
    staging = output.with_name(output.name + '_staging')
    if output.exists() or staging.exists():
        raise FileExistsError(f'Choose a new output directory: {output}')
    (staging / 'same_genotype').mkdir(parents=True)
    (staging / 'results.json').write_text(json.dumps(result, ensure_ascii=False, allow_nan=False), encoding='utf-8')
    write_csv(staging / 'wells.csv', result['wells'])
    write_csv(staging / 'biological_replicates.csv', result['biological_replicates'])
    write_csv(staging / 'orientation_summaries.csv', result['orientation_summaries'])
    write_csv(staging / 'dye_agreement.csv', result['dye_agreement'])
    write_csv(staging / 'replicate_set_comparisons.csv', result['replicate_set_comparisons'])
    write_csv(staging / 'background_summaries.csv', result['background_summaries'])
    write_csv(staging / 'dye_controls.csv', result['dye_controls'])
    write_csv(staging / 'source_manifest.csv', result['source_manifest'])
    write_csv(staging / 'identity_review.csv', [
        {k: w.get(k) for k in ('source_file', 'original_filename', 'plate', 'well', 'nominal_role', 'nominal_group',
                               'nominal_orientation', 'nominal_genotype', 'nominal_replicate', 'analysis_sample',
                               'replicate', 'technical_replicate', 'genotype', 'identity_status')} for w in result['wells']])
    write_csv(staging / 'identity_size_fingerprint.csv', result['identity_checks']['size_fingerprint'])
    write_csv(staging / 'identity_same_genotype.csv', result['identity_checks']['same_genotype'])
    # Same layout as the D4 same_genotype/wells.csv: D0 rows and this timepoint's rows.
    write_csv(staging / 'same_genotype' / 'wells.csv', same_genotype_rows)
    (staging / 'settings_comparison.json').write_text(json.dumps(result['settings_comparison'], indent=2), encoding='utf-8')
    (staging / 'report.md').write_text('\n'.join(report_markdown(result)) + '\n', encoding='utf-8')
    here = Path(__file__)
    provenance = {str(p.relative_to(REPO)): digest(p) for p in
                  (here, here.with_name('growth_assay_followup.py'), here.with_name('growth_assay_d0_report.py'),
                   here.with_name('growth_assay_label_review.py'), here.with_name('growth_assay_unlabelled_gate.py'),
                   here.with_name('competition_report.py'), here.with_name('gating.py'), here.with_name('fcs_io.py'),
                   here.with_name('longitudinal_report.py'))}
    (staging / 'code_manifest.json').write_text(json.dumps(provenance, indent=2), encoding='utf-8')
    # Dropbox or antivirus can hold a new directory on Windows; never replace an existing report.
    for attempt in range(4):
        if output.exists():
            raise FileExistsError(output)
        try:
            staging.rename(output)
            return
        except PermissionError:
            time.sleep(0.5 * (attempt + 1))
    # The directory itself stayed locked: move its entries into a new directory instead.
    output.mkdir()
    for entry in staging.iterdir():
        entry.rename(output / entry.name)
    for attempt in range(6):
        try:
            staging.rmdir()
            return
        except PermissionError:
            time.sleep(1 + attempt)
    print(f'Outputs are complete; remove the empty folder {staging} before the next run')


def collect(args: argparse.Namespace) -> dict:
    if args.output.exists() or args.output.with_name(args.output.name + '_staging').exists():
        raise FileExistsError('Output already exists; use --output with a new directory')
    labels = json.loads(args.labels.read_text(encoding='utf-8'))
    gating = json.loads(args.gating.read_text(encoding='utf-8'))
    baseline = json.loads(BASELINE.read_text(encoding='utf-8'))
    previous_path = args.previous / 'results.json'
    previous = json.loads(previous_path.read_text(encoding='utf-8'))
    if previous['baseline_sha256'] != digest(BASELINE) or previous['identity_status'] != 'user_confirmed':
        raise ValueError('The previous timepoint must use the same D0 baseline and approved labels')
    previous_primary = set(previous.get('gating', {}).get('unlabelled_gate_primary_backgrounds', []))
    if previous_primary != set(gating['unlabelled_gate_primary_backgrounds']):
        raise ValueError('Primary gating differs from the previous timepoint; comparisons would mix gates')
    cells = {}
    result = analyze(args.input, BASELINE, labels['timepoint'], cells, labels['technical_replicates'])
    if result['source'] != labels['source']:
        raise ValueError('Identity mapping belongs to a different export')
    tp = result['timepoint']
    assign_identities(result['wells'], labels)
    check_design(result['wells'], baseline)
    result['wells'].sort(key=lambda w: well_order(w['well']))
    d0, d0_hashes = load_d0(baseline)
    gate = measure(result, baseline, d0, cells, gating)
    for bg in gate['primary_backgrounds']:
        recorded = previous['unlabelled_gate']['rule']['mcherry_below'][bg]['D0']
        if recorded != gate['cuts'][bg]['D0']:
            raise ValueError(f'D0 valley for {bg} differs from the previous timepoint ({recorded})')
    result.update(summarise(result['wells'], previous))
    flag_wells(result['wells'], result['biological_replicates'])
    result['identity_checks'] = identity_checks(result['wells'], cells, baseline)
    previous_gates = previous['gate_definitions']
    result['gate_changes'] = gate_changes(result['gate_definitions'], previous_gates)
    if any(v['status'] != 'unchanged' for v in result['gate_changes'].values()):
        result['reference_gate_check'] = reference_gate_check(
            result['wells'], baseline, next((REPO / previous['source']).glob('*.cef')))
    d0_date = datetime.strptime(result['acquisition_metadata']['d0_fcs_date'], '%d-%b-%Y')
    dates = {datetime.strptime(d, '%d-%b-%Y') for d in result['acquisition_metadata'][f'{tp.lower()}_fcs_dates']}
    if len(dates) != 1:
        raise ValueError('Expected a single acquisition date')
    changed = {c['parameter'].rsplit('-', 1)[0] for s in result['settings_comparison'] for c in s['differences']}
    lightloss = next((c for s in result['settings_comparison'] for c in s['differences']
                      if c['parameter'] == 'LightLoss (Violet)-A'), None)
    approved = bool(labels.get('approval'))
    result.update(
        schema_version=3, title=f'{tp} competition assay — technical replicates, D0 and {previous["timepoint"]} comparison',
        identity_status='user_confirmed' if approved else 'recorded_labels_checked',
        identity_mapping=labels, identity_mapping_path=args.labels.resolve().relative_to(REPO).as_posix(),
        identity_mapping_sha256=digest(args.labels), elapsed_days_from_d0=(dates.pop() - d0_date).days,
        previous_timepoint=previous['timepoint'], previous_report=previous_path.relative_to(REPO).as_posix(),
        previous_sha256=digest(previous_path), d0_source_hashes=d0_hashes,
        settings_summary=(f'{len(changed)} detector settings differ from D0 in the metadata'
                          + (f" (e.g. LightLoss (Violet) {lightloss['d0_value']} → {lightloss[f'{tp.lower()}_value']})"
                             if lightloss else '')
                          + '; see settings_comparison.json.'),
        unlabelled_gate=dict(rule=dict(gfp_below=baseline['negative_rules']['broad']['gfp_below'],
                                       primary_mcherry_below=baseline['negative_rules']['broad']['mcherry_below'],
                                       valley_settings=VALLEY, mcherry_below=gate['cuts'],
                                       adjusted_backgrounds=[bg for bg, v in gate['valleys'].items() if v['detected']]),
                             valleys=gate['valleys'], cuts=gate['cuts'], primary_backgrounds=gate['primary_backgrounds']),
        gating=dict(path=args.gating.resolve().relative_to(REPO).as_posix(), sha256=digest(args.gating),
                    approval=gating['approval'], unlabelled_gate_primary_backgrounds=gate['primary_backgrounds']),
        d0_references=gate['references'], same_genotype_references=gate['same_genotype_references'],
        flag_thresholds=FLAGS,
        interpretation=[
            'Labels are the recorded acquisition names unless the identity mapping lists an approved correction. '
            'Fluorescence does not identify replicate numbers; replicate k is linked to the earlier replicate k by number.',
            'Change in log-odds is measured against the pooled D0 mix of each mixing batch. Biological replicates are '
            'post-split cultures of one D0 mix, so their spread does not include D0 mixing error.',
            'Technical replicates measure the same culture; they are averaged before any spread is computed.',
            'Direction agreement and replicate-set comparisons are descriptive. The C1-labelled arms received no drug '
            '(experimenter, 2026-10-29); they are second replicate sets, and only WT C1 is an independent D0 mix.',
            'Same-genotype controls (WT background) track reporter-subline growth; their shifts are not subtracted.',
            'No single-colour or unstained controls were acquired at this timepoint; reporter negativity does not '
            'establish viability. Detector settings differ from D0 in the metadata.'])
    result['qc'].update(identity_checks_consistent=all(r['consistent'] for r in result['identity_checks']['same_genotype'])
                        and all(r['consistent'] for r in result['identity_checks']['size_fingerprint'] if r['role'] == 'mixture'),
                        technical_replicates=labels['technical_replicates'],
                        flagged_wells=[w['well'] for w in result['wells'] if w['qc_flags']])
    rows = [dict(day=0, genotype=w['genotype'], replicate=w['replicate'], technical_replicate=1, well=w['well'],
                 sample=w['sample'], original_filename=f"{w['sample']}.fcs",
                 **measure_reporters(d0[w['well']], baseline))
            for w in baseline['wells'] if w['role'] == 'dye_neutral']
    rows += [dict(day=int(tp[1:]), genotype=w['genotype'], replicate=w['replicate'],
                  technical_replicate=w['technical_replicate'], well=w['well'], sample=w['analysis_sample'],
                  original_filename=w['original_filename'], **measure_reporters(cells[w['original_sample']]['reporters'], baseline))
             for w in result['wells'] if w['role'] == 'dye_neutral']
    # Every source, including D0, must be unchanged before anything is published.
    source = REPO / result['source']
    if inventory([p for p in source.rglob('*') if p.is_file()], REPO) != result['source_manifest']:
        raise ValueError('An input changed during processing; results were not published')
    if any(digest(REPO / rel) != h for rel, h in d0_hashes.items()) or digest(previous_path) != result['previous_sha256']:
        raise ValueError('A D0 source or the previous report changed during processing')
    write_outputs(result, args.output, rows)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--input', type=Path, default=DEFAULT_INPUT)
    parser.add_argument('--labels', type=Path, default=DEFAULT_LABELS)
    parser.add_argument('--gating', type=Path, default=DEFAULT_GATING)
    parser.add_argument('--previous', type=Path, default=DEFAULT_PREVIOUS, help='Approved report of the previous timepoint')
    parser.add_argument('--output', type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    result = collect(args)
    print(json.dumps(result['qc'], indent=2))
    print(args.output.resolve() / 'report.md')


if __name__ == '__main__':
    main()
