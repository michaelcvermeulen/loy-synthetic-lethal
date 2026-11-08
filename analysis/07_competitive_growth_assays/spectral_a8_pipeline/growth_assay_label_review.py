"""Apply the experimenter's D4 identity record and assess descriptive consistency.

Replicate numbers cannot be inferred from fluorescence. Keep raw nominal fields
and every measurement intact, and compare each orientation with its shared D0 mix.
Same-genotype controls may also carry an approved genotype correction.
"""
from __future__ import annotations

from collections import Counter
from copy import deepcopy
from datetime import datetime, timezone
import math
import re

import numpy as np

GROUPS = ('WT', 'WT_C1', 'DDX3X', 'DDX3X_C1', 'SREBF1', 'TADA2B')
ORIENTATIONS = ('LOY_GFP', 'LOY_mCherry')


def log_odds_change(p: float | None, p0: float | None) -> float | None:
    # No pseudocount or clipped infinite effect is silently introduced.
    if p is None or p0 is None or not 0 < p < 1 or not 0 < p0 < 1:
        return None
    return math.log(p / (1 - p)) - math.log(p0 / (1 - p0))


def stats(values: list[float | None]) -> dict:
    valid = [x for x in values if x is not None]
    return dict(n=len(valid), mean=float(np.mean(valid)) if valid else None,
                sd=float(np.std(valid, ddof=1)) if len(valid) > 1 else None,
                minimum=min(valid) if valid else None, maximum=max(valid) if valid else None)


def strict_baseline(baseline: dict, b: dict) -> float:
    wells = [w for w in baseline['wells'] if w['well'] in b['wells'].split(';')]
    if len(wells) != len(b['wells'].split(';')):
        raise ValueError('Missing D0 source well')
    denominators = [round(w['primary_singlets'] * (1 - w['negative_strict_fraction'])) for w in wells]
    return sum(n * w['loy_fraction_excluding_strict_negatives'] for n, w in zip(denominators, wells)) / sum(denominators)


def direction(values: list[float | None]) -> str:
    if any(v is None for v in values):
        return 'unavailable'
    if all(v > 0 for v in values):
        return 'both increase'
    if all(v < 0 for v in values):
        return 'both decrease'
    return 'opposite or zero directions'


def summarise_arms(rows: list[dict]) -> tuple[list[dict], list[dict]]:
    """Orientation summaries and dye agreement from the current per-well fields."""
    summaries, agreement = [], []
    for group in GROUPS:
        arms = []
        for orientation in ORIENTATIONS:
            ws = sorted((w for w in rows if w['role'] == 'mixture' and w['group'] == group
                         and w['orientation'] == orientation), key=lambda w: w['replicate'])
            s = dict(group=group, orientation=orientation,
                     wells_by_replicate=';'.join(w['well'] for w in ws),
                     d0_loy_fraction=ws[0]['d0_loy_fraction'],
                     d0_loy_fraction_strict=ws[0]['d0_loy_fraction_strict'])
            for key in ('loy_fraction', 'excluded_broad_fraction', 'delta_log_odds', 'delta_log_odds_strict'):
                s.update({key + '_' + k: v for k, v in stats([w[key] for w in ws]).items()})
            s['replicate_range_pp'] = 100 * (s['loy_fraction_maximum'] - s['loy_fraction_minimum']) if s['loy_fraction_n'] else None
            summaries.append(s)
            arms.append(s)
        deltas = [a['delta_log_odds_mean'] for a in arms]
        dstrict = [a['delta_log_odds_strict_mean'] for a in arms]
        agreement.append(dict(group=group, gfp_delta_log_odds=deltas[0], mcherry_delta_log_odds=deltas[1],
                              direction_agreement=direction(deltas), strict_direction_agreement=direction(dstrict),
                              orientation_gap_log_odds=abs(deltas[0] - deltas[1]) if None not in deltas else None))
    return summaries, agreement


def apply_labels_and_compare(imported: dict, baseline: dict, labels: dict) -> dict:
    if imported['source'] != labels['source']:
        raise ValueError('Identity mapping belongs to a different export')
    if imported['classifier'] != {k: baseline['classifier'][k] for k in ('weights', 'intercept')}:
        raise ValueError('Cached reporter classifier differs from D0')
    if imported['negative_rules'] != baseline['negative_rules'] or imported['primary_negative_rule'] != 'broad':
        raise ValueError('Cached negative rules differ from approved D0 rules')
    result = deepcopy(imported)
    rows = result['wells']
    assignments = {a['original_sample']: a['replicate'] for a in labels['assignments']}
    genotypes = {a['original_sample']: a['genotype'] for a in labels['assignments'] if 'genotype' in a}
    if len(assignments) != len(labels['assignments']):
        raise ValueError('Duplicate identity mapping entry')
    if not set(assignments) <= {w['original_sample'] for w in rows}:
        raise ValueError('Identity mapping contains an unknown recording')
    if any(type(r) is not int or r not in (1, 2, 3) for r in assignments.values()):
        raise ValueError('Expected replicate 1, 2 or 3')
    # Genotype corrections apply only to same-genotype controls; a competition well's genotype is its orientation.
    neutral_samples = {w['original_sample'] for w in rows if w['nominal_role'] == 'dye_neutral'}
    if not set(genotypes) <= neutral_samples or not set(genotypes.values()) <= {'CTRL', 'LOY'}:
        raise ValueError('Genotype corrections must name CTRL or LOY for same-genotype controls')
    baselines = {(b['group'], b['orientation']): b for b in baseline['d0_baselines']}
    strict = {key: strict_baseline(baseline, b) for key, b in baselines.items()}
    for w in rows:
        for field in ('role', 'group', 'background', 'condition', 'orientation', 'genotype', 'reporter', 'replicate'):
            w[field] = w['nominal_' + field]
        w['replicate'] = assignments.get(w['original_sample'], w['nominal_replicate'])
        w['replicate_changed'] = w['replicate'] != w['nominal_replicate']
        w['genotype'] = genotypes.get(w['original_sample'], w['nominal_genotype'])
        w['genotype_changed'] = w['genotype'] != w['nominal_genotype']
        w['identity_status'] = 'corrected' if w['replicate_changed'] or w['genotype_changed'] else 'user_confirmed'
        w['analysis_sample'] = re.sub(r'_[1-3](?:_\d{3})?$', '_' + str(w['replicate']), w['original_sample'])
        if w['genotype_changed']:
            w['analysis_sample'] = w['analysis_sample'].replace(
                f"-{w['nominal_genotype']}_GFP_CHERRY_", f"-{w['genotype']}_GFP_CHERRY_", 1)
            w['correction_note'] = (f"Genotype corrected from {w['nominal_genotype']} to {w['genotype']} (approved swap of "
                                    'the row-C same-genotype blocks); replicate follows the physical well.')
        else:
            w['correction_note'] = ('Replicate reassigned from experimental record; group identity unchanged.'
                                    if w['replicate_changed'] else 'Recorded identity accepted by experimenter.')
        if w['role'] == 'mixture':
            b = baselines[(w['group'], w['orientation'])]
            w['d0_loy_fraction'] = b['d0_loy_fraction']
            w['d0_batch'] = b['mixture_batch']
            w['d0_wells'] = b['wells']
            w['loy_fraction'] = w['gfp_fraction'] if w['orientation'] == 'LOY_GFP' else w['mcherry_fraction']
            p = w['gfp_fraction_strict']
            w['loy_fraction_strict'] = p if p is None or w['orientation'] == 'LOY_GFP' else 1 - p
            w['d0_loy_fraction_strict'] = strict[(w['group'], w['orientation'])]
            w['delta_log_odds'] = log_odds_change(w['loy_fraction'], w['d0_loy_fraction'])
            w['delta_log_odds_strict'] = log_odds_change(w['loy_fraction_strict'], w['d0_loy_fraction_strict'])
    design = Counter((w['group'], w['orientation'], w['replicate']) for w in rows if w['role'] == 'mixture')
    expected = {(g, o, r) for g, o in baselines for r in (1, 2, 3)}
    neutral = Counter((w['genotype'], w['replicate']) for w in rows if w['role'] == 'dye_neutral')
    if set(design) != expected or set(design.values()) != {1}:
        raise ValueError('Corrected competition design has missing or duplicate replicates')
    if set(neutral) != {(g, r) for g in ('CTRL', 'LOY') for r in (1, 2, 3)} or set(neutral.values()) != {1}:
        raise ValueError('Corrected dye-control design has missing or duplicate replicates')

    summaries, agreement = summarise_arms(rows)
    controls = []
    for genotype in ('CTRL', 'LOY'):
        d0 = [w for w in baseline['wells'] if w['role'] == 'dye_neutral' and w['genotype'] == genotype]
        d4 = [w for w in rows if w['role'] == 'dye_neutral' and w['genotype'] == genotype]
        p0 = sum(w['gfp_count'] for w in d0) / sum(w['reporter_positive_singlets'] for w in d0)
        s = dict(genotype=genotype, d0_gfp_fraction_pooled=p0, d0_wells=';'.join(w['well'] for w in d0))
        s.update({'d4_gfp_fraction_' + k: v for k, v in stats([w['gfp_fraction'] for w in d4]).items()})
        s['mean_delta_log_odds'] = stats([log_odds_change(w['gfp_fraction'], p0) for w in d4])['mean']
        controls.append(s)
    result.update(schema_version=2, title='D4 competition assay — corrected labels and D0 comparison',
                  reviewed_utc=datetime.now(timezone.utc).isoformat(), identity_status='user_confirmed',
                  identity_mapping=deepcopy(labels), orientation_summaries=summaries,
                  dye_agreement=agreement, dye_controls=controls)
    genotype_changes = sum(w['genotype_changed'] for w in rows)
    result['qc'].update(corrected_replicate_count=sum(w['replicate_changed'] for w in rows),
                        corrected_genotype_count=genotype_changes,
                        accepted_identities=len(rows) - genotype_changes, corrected_design_complete=True)
    result['interpretation'] = [
        'Replicate numbering follows the experimenter\'s record. Fluorescence does not independently identify replicate 1, 2 or 3.',
        'Replicate corrections reorder wells within groups, so group means, replicate spread and dye-orientation comparisons are unchanged.',
        *(['Same-genotype genotype corrections follow the approved identity check: per-cell size and reporter-brightness '
           'profiles of the relabelled wells match known LOY or CTRL cells on the same plate. Competition labels are unchanged.']
          if genotype_changes else []),
        'D0 references pool retained counts within each original mixing batch and dye orientation. D4 summaries weight each replicate well equally.',
        'Change in log-odds = ln(LOY/CTRL) at D4 minus ln(LOY/CTRL) at D0. Positive values mean increased LOY representation among retained cells.',
        'Direction agreement is descriptive, not a statistical test or proof of identity. D0 splits share a mixture; their uncertainty is not independent for each D4 well.',
        'Same-genotype controls track reporter balance in WT backgrounds only. Their shifts are not automatically subtracted from knockout groups.',
        'No D4 single-colour or unstained controls were acquired. Detector settings differ from D0, and reporter negativity does not establish viability.',
        'C1 is an arm label. Treatment timing and dose are not supplied, so comparisons do not establish an inhibitor effect.'
    ]
    return result
