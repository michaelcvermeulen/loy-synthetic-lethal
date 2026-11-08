"""Load a growth-assay follow-up with the frozen D0 reporter model.

Sample identities are assigned from an explicit experimenter-approved mapping.
Raw acquisitions are never renamed or rewritten. Outputs are staged, sources are
hashed before and after analysis, and existing output directories are refused.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
import time
import zipfile

import numpy as np

from .competition_report import LDAModel, load_cef_gates
from .fcs_io import load_fcs
from .gating import points_in_polygon
from .growth_assay_d0_report import parse_sample
from .longitudinal_report import verify_unmixing_compatibility
from .growth_assay_label_review import apply_labels_and_compare, direction, summarise_arms

REPO = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = REPO / '20261020_173018_Export' / 'comp_growth_assay_D4_OCT20'
DEFAULT_BASELINE = REPO / 'analysis/spectral_growth_assay_D0_20261015/results.json'
DEFAULT_OUTPUT = REPO / 'analysis/spectral_growth_assay_D4_20261020'
GATE_NAMES = ('Scatter', 'SSC Singlets', 'FSC Singlets')
# Scatter parameters kept per cell for identity checks; their scale is stable between days.
CELL_SCATTER = ('LightLoss (Violet)-A', 'SSC (Imaging)-A')
TEMPLATE = Path(__file__).with_name('growth_assay_followup.html')
DEFAULT_LABELS = Path(__file__).with_name('growth_assay_D4_labels.json')
DEFAULT_GATING = Path(__file__).with_name('growth_assay_D4_gating.json')


def digest(path: Path) -> str:
    h = sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def inventory(paths: list[Path], root: Path) -> list[dict]:
    return [dict(path=p.relative_to(root).as_posix(), size_bytes=p.stat().st_size,
                 sha256=digest(p)) for p in sorted(paths)]


def read_statistics(path: Path) -> dict:
    result, current = {}, None
    with path.open(encoding='utf-8-sig', newline='') as handle:
        for row in csv.reader(handle):
            if row and row[0] == 'Recording Name':
                current = row[1]
                if current in result:
                    raise ValueError(f'Duplicate statistics recording: {current}')
                result[current] = {}
            elif current and row and row[0] in ('All Events', *GATE_NAMES):
                result[current][row[0]] = int(row[1])
    return result


def measure_reporters(features: np.ndarray, baseline: dict) -> dict:
    """Fractions have an explicit retained-event denominator, including at zero."""
    values = np.asarray(features, dtype=float)
    if values.ndim != 2 or values.shape[1] != 2 or not np.isfinite(values).all():
        raise ValueError('Reporter features must be finite N x 2 values')
    if baseline['primary_negative_rule'] != 'broad':
        raise ValueError('Expected the approved broad-negative D0 analysis')
    rule = baseline['negative_rules']['broad']
    excluded = (values[:, 0] < rule['gfp_below']) & (values[:, 1] < rule['mcherry_below'])
    model = LDAModel(np.array(baseline['classifier']['weights']), baseline['classifier']['intercept'])
    classes = model.predict(values).astype(bool)
    kept = ~excluded
    total, retained = len(values), int(kept.sum())
    green = int((classes & kept).sum())
    strict = baseline['negative_rules']['strict']
    keep_strict = ~((values[:, 0] < strict['gfp_below']) & (values[:, 1] < strict['mcherry_below']))
    return dict(primary_singlets=total, reporter_positive_singlets=retained,
                excluded_broad_count=int(excluded.sum()),
                excluded_broad_fraction=float(excluded.mean()) if total else None,
                gfp_count=green, mcherry_count=retained - green,
                gfp_fraction=green / retained if retained else None,
                mcherry_fraction=(retained - green) / retained if retained else None,
                gfp_fraction_strict=float(classes[keep_strict].mean()) if keep_strict.any() else None,
                gfp_fraction_all_singlets=float(classes.mean()) if total else None)


def provisional_identity(stem: str, baseline: dict) -> dict:
    """Match by the recorded biological label; physical D0/D4 wells differ."""
    sample = parse_sample(stem)
    result = dict(identity_status='pending_correction', nominal_role=sample.role,
                  nominal_group=sample.group, nominal_background=sample.background,
                  nominal_condition=sample.condition, nominal_orientation=sample.orientation,
                  nominal_genotype=sample.genotype, nominal_reporter=sample.dye,
                  nominal_replicate=sample.replicate)
    key = (sample.group, sample.orientation)
    lookup = {(b['group'], b['orientation']): b for b in baseline['d0_baselines']}
    if sample.role == 'mixture':
        b = lookup[key]
        result.update(d0_loy_fraction_provisional=b['d0_loy_fraction'],
                      d0_batch_provisional=b['mixture_batch'], d0_wells_provisional=b['wells'])
    return result


def gate_signature(gates: list[dict]) -> dict:
    return {g['name']: dict(parameters=[p['name'] for p in g['parameters']],
                           vertices=g['vertices'], parent=g['parentPopulationId'],
                           gate_id=g['gateId'], children=g['children'])
            for g in gates if g['name'] in GATE_NAMES and g['gateKind'] == 'Polygon'}


def voltages(metadata: dict) -> dict:
    return {metadata[f'p{i}n']: metadata.get(f'p{i}v')
            for i in range(1, int(metadata['par']) + 1)}


def analyze(export: Path, baseline_path: Path, timepoint: str = 'D4', cells: dict | None = None,
            technical_replicates: int = 1) -> dict:
    """Import one follow-up export. ``cells``, when given, receives each recording's gated per-cell
    reporter and scatter values keyed by FCS stem; ``technical_replicates`` is the expected number of
    wells per mixture replicate."""
    export, baseline_path = export.resolve(), baseline_path.resolve()
    day = timepoint.lower()
    files = sorted(export.rglob('*.fcs'))
    if not files:
        raise ValueError(f'No FCS recordings in {export}')
    cef_paths, stat_paths = list(export.glob('*.cef')), list(export.glob('*Statistics.csv'))
    if len(cef_paths) != 1 or len(stat_paths) != 1:
        raise ValueError('Expected exactly one CEF and Statistics CSV')
    cef, stats_path = cef_paths[0], stat_paths[0]
    sources_before = inventory([p for p in export.rglob('*') if p.is_file()], REPO)
    baseline_hash = digest(baseline_path)
    baseline = json.loads(baseline_path.read_text(encoding='utf-8'))
    reference_file = next(REPO.joinpath(baseline['source']).rglob('*.fcs'))
    reference_hash = digest(reference_file)
    expected_hashes = {m['path']: m['sha256'] for m in baseline['source_manifest']}
    if expected_hashes.get(reference_file.relative_to(REPO).as_posix()) != reference_hash:
        raise ValueError('D0 reference acquisition does not match its recorded hash')
    reference = load_fcs(reference_file).metadata
    gates = load_cef_gates(cef)
    with zipfile.ZipFile(cef) as archive:
        panel = json.loads(archive.read('experiment.json'))['experiment']['panels'][0]
    raw_gates = panel['analysis']['gates']
    signature = gate_signature(raw_gates)
    if set(signature) != set(GATE_NAMES):
        raise ValueError('Missing expected scatter/singlet polygon')
    parent = '0-1'
    for name in GATE_NAMES:
        if signature[name]['parent'] != parent or len(signature[name]['children']) != 1:
            raise ValueError(f'Unexpected gate hierarchy at {name}')
        parent = signature[name]['children'][0]['populationId']
    stats = read_statistics(stats_path)
    if set(stats) != {p.stem.split('-', 1)[1] for p in files}:
        raise ValueError('Statistics recordings do not match the FCS inventory')
    wells, settings, points = [], [], {}
    d0_voltages = voltages(reference)
    for path in files:
        print(f'Loading {path.name}', flush=True)
        sample = parse_sample(path.stem)
        with zipfile.ZipFile(path.with_suffix('.zip')) as archive:
            recorded_gates = json.loads(archive.read('recording-data/analysis-model.json'))['analysisModel']['gates']
        if gate_signature(recorded_gates) != signature:
            raise ValueError(f'Recording and CEF gates differ: {path.name}')
        data = load_fcs(path)
        metadata, events = data.metadata, data.events
        for key in ('bdspectral matrix name', 'spill', 'cytometer configuration name'):
            if not reference.get(key) or not metadata.get(key):
                raise ValueError(f'Missing compatibility metadata: {key}')
        verify_unmixing_compatibility(reference, metadata)
        if metadata['cytometer configuration name'] != reference['cytometer configuration name']:
            raise ValueError(f'Optical configuration differs: {path.name}')
        keep = np.ones(len(events), dtype=bool)
        counts = {'All Events': len(events)}
        for name in GATE_NAMES:
            g = gates[name]
            keep &= points_in_polygon(events[list(g.parameters)].to_numpy(float), g.vertices)
            counts[name] = int(keep.sum())
        if counts != stats[path.stem.split('-', 1)[1]]:
            raise ValueError(f'Instrument gate-count mismatch: {path.name}')
        features = np.arcsinh(events.loc[keep, ['eGFP-A', 'mCherry-A']].to_numpy(float) / 1000.0)
        if cells is not None:
            cells[path.stem] = dict(reporters=features, scatter=events.loc[keep, list(CELL_SCATTER)].to_numpy(float))
        row = dict(source_file=path.relative_to(REPO).as_posix(), original_filename=path.name,
                   original_sample=path.stem, plate=path.parent.name, well=sample.well,
                   timepoint=timepoint, acquired_date=metadata.get('date'),
                   acquisition_start_utc=metadata.get('begindatetime'),
                   total_events=len(events), scatter_count=counts['Scatter'],
                   ssc_singlet_count=counts['SSC Singlets'], instrument_singlets=counts['FSC Singlets'],
                   gate_counts_match=True, recording_gates_match=True, unmixing_matches_d0=True)
        row.update(measure_reporters(features, baseline))
        row.update(provisional_identity(path.stem, baseline))
        if row['nominal_role'] == 'mixture':
            key = 'gfp_fraction' if row['nominal_orientation'] == 'LOY_GFP' else 'mcherry_fraction'
            row[f'{day}_loy_fraction_provisional'] = row[key]
        row['saturated_fraction'] = float((events.loc[keep, 'Saturated'] > 0).mean()) if keep.any() else None
        row['singlet_yield'] = counts['FSC Singlets'] / len(events) if len(events) else None
        v = voltages(metadata)
        changes = [{'parameter': n, 'd0_value': d0_voltages.get(n), f'{day}_value': value}
                   for n, value in v.items() if value != d0_voltages.get(n)]
        row['voltage_metadata_differences'] = len(changes)
        settings.append(dict(original_filename=path.name, differences=changes))
        # Values are un-clipped; subsampling is solely for interactive display.
        indexes = np.linspace(0, len(features) - 1, min(1500, len(features)), dtype=int)
        points[path.stem] = np.round(features[indexes], 4).tolist()
        wells.append(row)
    if len({(w['plate'], w['well']) for w in wells}) != len(wells):
        raise ValueError('Duplicate physical well identifiers')
    # Verify run completeness from label roles, without asserting their identities are correct.
    groups = Counter((w['nominal_group'], w['nominal_orientation'], w['nominal_replicate'])
                     for w in wells if w['nominal_role'] == 'mixture')
    expected = {(b['group'], b['orientation'], rep) for b in baseline['d0_baselines'] for rep in (1, 2, 3)}
    complete = set(groups) == expected and all(n == technical_replicates for n in groups.values())
    sources_after = inventory([p for p in export.rglob('*') if p.is_file()], REPO)
    if sources_before != sources_after or digest(baseline_path) != baseline_hash or digest(reference_file) != reference_hash:
        raise ValueError('An input changed during processing; results were not published')
    return dict(schema_version=1, title=f'{timepoint} competition assay — identity review',
                created_utc=datetime.now(timezone.utc).isoformat(),
                source=export.relative_to(REPO).as_posix(), identity_status='pending_correction',
                timepoint=timepoint, source_manifest=sources_before, source_hashes_verified=True,
                baseline_path=baseline_path.relative_to(REPO).as_posix(), baseline_sha256=baseline_hash,
                reference_file=reference_file.relative_to(REPO).as_posix(), reference_sha256=reference_hash,
                classifier={k: baseline['classifier'][k] for k in ('weights', 'intercept')},
                classifier_training=baseline['classifier']['training'],
                transform='arcsinh(value / 1000)', negative_rules=baseline['negative_rules'],
                primary_negative_rule=baseline['primary_negative_rule'],
                unmixing_matrix=reference['bdspectral matrix name'],
                spill_sha256=sha256(reference['spill'].encode()).hexdigest(),
                optical_configuration=reference['cytometer configuration name'],
                gate_definitions=signature, wells=wells, plot_points=points,
                d0_baselines=baseline['d0_baselines'],
                acquisition_metadata={'d0_fcs_date': reference.get('date'),
                                      f'{day}_fcs_dates': sorted({w['acquired_date'] for w in wells})},
                settings_comparison=settings,
                qc=dict(recordings=len(wells), total_events=sum(w['total_events'] for w in wells),
                        primary_singlets=sum(w['primary_singlets'] for w in wells),
                        retained_events=sum(w['reporter_positive_singlets'] for w in wells),
                        role_counts=dict(Counter(w['nominal_role'] for w in wells)),
                        gate_counts_matched=len(wells), recorded_gates_matched=len(wells),
                        unmixing_matched=len(wells), nominal_mixture_design_complete=complete,
                        empty_retained_wells=[w['well'] for w in wells if not w['reporter_positive_singlets']]))


def attach_unlabelled_gate(result: dict, labels_path: Path, path: Path, gating_path: Path | None = None) -> None:
    """Merge the day-specific unlabelled-cell gate.

    Backgrounds named in the approved gating file take the unlabelled-cell gate as their primary
    rule; standard-gate values are kept for every well as a sensitivity. Without a gating file the
    gate is merged as sensitivity fields only.
    """
    if not path.exists():
        return
    sensitivity = json.loads(path.read_text(encoding='utf-8'))
    gating = json.loads(gating_path.read_text(encoding='utf-8')) if gating_path and gating_path.exists() else None
    manifest = {m['path']: m['sha256'] for m in result['source_manifest']}
    d4_sources = {k: v for k, v in sensitivity['source_hashes'].items() if k.startswith(result['source'] + '/')}
    if (sensitivity['identity_mapping_sha256'] != digest(labels_path)
            or sensitivity['baseline_sha256'] != result['baseline_sha256']
            or not d4_sources or any(manifest.get(k) != v for k, v in d4_sources.items())):
        raise ValueError('Unlabelled-gate sensitivity was built from different labels, baseline or sources; rerun it')
    if any('primary_gate' in w for w in result['wells']):
        raise ValueError('Gating already applied; start from the raw measurement cache')
    rule = sensitivity['rule']
    primary_backgrounds = set(gating['unlabelled_gate_primary_backgrounds']) if gating else set()
    if gating and gating['source'] != result['source']:
        raise ValueError('Gating decision belongs to a different export')
    if not primary_backgrounds <= set(rule['adjusted_backgrounds']):
        raise ValueError('Unlabelled-cell gate is primary only where a separating valley was detected on both days')
    by_well = {w['well']: w for w in sensitivity['wells']}
    for w in result['wells']:
        s = by_well[w['well']]
        if s['original_sample'] != w['original_sample'] or s['primary_singlets'] != w['primary_singlets']:
            raise ValueError('Unlabelled-gate sensitivity wells do not match this export')
        w['excluded_fraction_unlabelled_gate'] = s['excluded_fraction']
        w['delta_log_odds_unlabelled_gate'] = s.get('delta_log_odds')
        for key in ('excluded_broad_fraction', 'gfp_fraction', 'loy_fraction', 'd0_loy_fraction', 'delta_log_odds'):
            if key in w:
                w[key + '_standard_gate'] = w[key]
        w['primary_gate'] = 'unlabelled_cell_gate' if w['background'] in primary_backgrounds else 'standard'
        if w['primary_gate'] == 'unlabelled_cell_gate':
            w.update(reporter_positive_singlets=s['retained'], excluded_broad_count=s['primary_singlets'] - s['retained'],
                     excluded_broad_fraction=s['excluded_fraction'], gfp_count=s['gfp_count'],
                     mcherry_count=s['retained'] - s['gfp_count'], gfp_fraction=s['gfp_fraction'],
                     mcherry_fraction=None if s['gfp_fraction'] is None else 1 - s['gfp_fraction'])
            if w['role'] == 'mixture':
                w.update(loy_fraction=s['loy_fraction'], d0_loy_fraction=s['d0_loy_fraction'],
                         delta_log_odds=s['delta_log_odds'])
    gated_baselines = {(b['group'], b['orientation']): b for b in sensitivity['d0_baselines']}
    for b in result['d0_baselines']:
        b['d0_loy_fraction_standard_gate'] = b['d0_loy_fraction']
        b['primary_gate'] = 'unlabelled_cell_gate' if b['group'].split('_')[0] in primary_backgrounds else 'standard'
        if b['primary_gate'] == 'unlabelled_cell_gate':
            g = gated_baselines[(b['group'], b['orientation'])]
            b.update(d0_loy_fraction=g['d0_loy_fraction'], loy_count=g['loy_count'], reporter_positive_singlets=g['retained'])
    summaries, agreement = summarise_arms(result['wells'])
    standard = {}
    for s in summaries:
        values = [w['delta_log_odds_standard_gate'] for w in result['wells'] if w['role'] == 'mixture'
                  and w['group'] == s['group'] and w['orientation'] == s['orientation']]
        s['delta_log_odds_standard_gate_mean'] = float(np.mean(values))
        standard[(s['group'], s['orientation'])] = s['delta_log_odds_standard_gate_mean']
    gate = {a['group']: a for a in sensitivity['dye_agreement']}
    for a in agreement:
        g, m = standard[(a['group'], 'LOY_GFP')], standard[(a['group'], 'LOY_mCherry')]
        a.update(gfp_delta_log_odds_standard_gate=g, mcherry_delta_log_odds_standard_gate=m,
                 standard_gate_direction_agreement=direction([g, m]),
                 gfp_delta_log_odds_unlabelled_gate=gate[a['group']]['gfp_delta_log_odds'],
                 mcherry_delta_log_odds_unlabelled_gate=gate[a['group']]['mcherry_delta_log_odds'],
                 primary_gate='unlabelled_cell_gate' if a['group'].split('_')[0] in primary_backgrounds else 'standard')
    result.update(orientation_summaries=summaries, dye_agreement=agreement)
    result['unlabelled_gate'] = dict(path=path.resolve().relative_to(REPO).as_posix(), sha256=digest(path),
                                     rule=rule, valleys=sensitivity['valleys'],
                                     primary_backgrounds=sorted(primary_backgrounds))
    if gating:
        result['gating'] = dict(path=gating_path.resolve().relative_to(REPO).as_posix(), sha256=digest(gating_path),
                                approval=gating['approval'], approved_date=gating['approved_date'],
                                unlabelled_gate_primary_backgrounds=sorted(primary_backgrounds))
    cuts = '; '.join(f"{bg} mCherry < {c['D0']} at D0 and < {c['D4']} at D4" for bg, c in rule['mcherry_below'].items()
                     if bg in rule['adjusted_backgrounds'])
    if primary_backgrounds:
        result['interpretation'].append(
            f"Primary gating is background-specific: {', '.join(sorted(primary_backgrounds))} uses the day-specific "
            f'unlabelled-cell gate ({cuts}), approved {gating["approved_date"]}, because an unlabelled population crosses '
            'the standard mCherry cut between D0 and D4. All other arms use the standard gate (eGFP < 6 and mCherry < 6.5). '
            'Standard-gate values for every arm are retained as a sensitivity.')
    else:
        result['interpretation'].append(
            'Sensitivity: a day-specific unlabelled-cell gate moves the mCherry cut to the valley between unlabelled and '
            f"mCherry-positive cells where such a cluster is detected on both days ({cuts or 'none detected'}); other "
            'backgrounds keep the primary cut. It is not the primary analysis.')


def write_csv(path: Path, rows: list[dict]) -> None:
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def unlabelled_gate_table(result: dict) -> list[str]:
    rule = result['unlabelled_gate']['rule']
    primary = result['unlabelled_gate']['primary_backgrounds']
    cuts = '; '.join(f"{bg}: D0 {c['D0']}, D4 {c['D4']}" for bg, c in rule['mcherry_below'].items()
                     if bg in rule['adjusted_backgrounds'])
    note = (f"Primary columns use the unlabelled-cell gate for {', '.join(primary)} ({cuts}) and the standard gate "
            'for every other arm. Sensitivity columns apply the standard gate (eGFP < 6 and mCherry < 6.5) to every arm.'
            if primary else
            f'Primary columns use the standard gate. The day-specific unlabelled-cell gate ({cuts}) is in '
            '`unlabelled_gate_sensitivity.json`.')
    return [
        '| Arm | Primary gate | Change in log-odds, LOY-GFP | Change in log-odds, LOY-mCherry | Direction agreement | '
        'Standard gate (sensitivity), LOY-GFP | Standard gate (sensitivity), LOY-mCherry |',
        '|---|---|---:|---:|---|---:|---:|',
        *[f"| {a['group']} | {'unlabelled-cell' if a['primary_gate'] == 'unlabelled_cell_gate' else 'standard'} | "
          f"{a['gfp_delta_log_odds']:.3f} | {a['mcherry_delta_log_odds']:.3f} | {a['direction_agreement']} | "
          f"{a['gfp_delta_log_odds_standard_gate']:.3f} | {a['mcherry_delta_log_odds_standard_gate']:.3f} |"
          for a in result['dye_agreement']], '',
        note,
    ]


def write_outputs(result: dict, output: Path) -> None:
    output = output.resolve()
    staging = output.with_name(output.name + '_staging')
    if output.exists() or staging.exists():
        raise FileExistsError(f'Choose a new output directory: {output}')
    staging.mkdir(parents=True)
    wells, qc = result['wells'], result['qc']
    payload = json.dumps(result, ensure_ascii=False, allow_nan=False)
    (staging / 'results.json').write_text(payload, encoding='utf-8')
    write_csv(staging / 'wells.csv', wells)
    write_csv(staging / 'source_manifest.csv', result['source_manifest'])
    write_csv(staging / 'd0_baselines.csv', result['d0_baselines'])
    write_csv(staging / 'identity_review.csv', [
        {**{k: w.get(k) for k in ('source_file', 'original_filename', 'plate', 'well',
                                'nominal_role', 'nominal_group', 'nominal_orientation',
                                'nominal_genotype', 'nominal_replicate', 'identity_status')},
         **{k: w[k] for k in ('analysis_sample', 'replicate', 'replicate_changed', 'genotype', 'genotype_changed',
                              'correction_note')}} for w in wells])
    write_csv(staging / 'orientation_summaries.csv', result['orientation_summaries'])
    write_csv(staging / 'dye_agreement.csv', result['dye_agreement'])
    write_csv(staging / 'dye_controls.csv', result['dye_controls'])
    (staging / 'identity_mapping.json').write_text(json.dumps(result['identity_mapping'], indent=2), encoding='utf-8')
    (staging / 'settings_comparison.json').write_text(
        json.dumps(result['settings_comparison'], indent=2), encoding='utf-8')
    html = TEMPLATE.read_text(encoding='utf-8').replace('__DATA__', payload.replace('<', '\\u003c'))
    (staging / 'report.html').write_text(html, encoding='utf-8')
    markdown = [
        '# D4 competition assay — corrected labels and D0 comparison', '',
        (f"**Experimenter-approved labels applied: {qc['corrected_replicate_count']} replicate numbers corrected; "
         f"{qc['corrected_genotype_count']} same-genotype control genotypes corrected (row C swap); "
         f"the other {qc['accepted_identities']} group identities accepted as recorded.**"
         if qc['corrected_genotype_count'] else
         f"**Experimenter-approved labels applied: {qc['corrected_replicate_count']} replicate numbers corrected; all {qc['accepted_identities']} group identities accepted.**"), '',
        f"Source: `{result['source']}`; acquired {', '.join(result['acquisition_metadata']['d4_fcs_dates'])}.",
        f"Loaded {qc['recordings']} recordings / {qc['total_events']:,} events; "
        f"{qc['primary_singlets']:,} primary singlets; {qc['retained_events']:,} retained after broad exclusion.", '',
        '- D4 CEF Scatter → SSC Singlets → FSC Singlets, verified against each recording ZIP.',
        f"- All {qc['gate_counts_matched']} recordings reproduce every gate count in the Statistics CSV.",
        '- Saved D0 LDA is frozen; no D4 retraining. Transform: arcsinh(value / 1000).',
        '- Broad negatives excluded: eGFP < 6 AND mCherry < 6.5 in transformed units.',
        '- GFP/mCherry fractions use retained singlets as denominator; exclusions use all primary singlets.',
        '- Same unmixing payload and optical configuration as D0; voltage metadata differ (see settings_comparison.json).',
        '- No D4 single-colour or unstained controls are present; reporter calibration is not independently validated.',
        '- No viability measurement is inferred from reporter negativity or raw detector intensity.',
        '- All export file sizes and SHA-256 hashes were recorded and rechecked after analysis.', '',
        '## Identity review', '',
        'Use report.html to compare D0 and D4, browse physical wells and inspect reporter plots. '
        'identity_mapping.json records the approved mapping; identity_review.csv preserves recorded and effective labels. '
        'nominal_* and *_provisional fields preserve the original import; effective labels and comparisons use the new unsuffixed fields.', '',
        'D0 matches use condition and dye orientation, not physical well position. DDX3X arms share their D0 mix; '
        'WT arms retain separate D0 mixtures. C1 remains the planned inhibitor arm; treatment timing/dose are '
        'not encoded by the acquisition label.', '',
        *(['Genotype evidence for the row-C swap: ' + result['identity_mapping']['genotype_evidence'], '']
          if 'genotype_evidence' in result['identity_mapping'] else []),
        '## Interpretation', '', *['- ' + t for t in result['interpretation']], '',
        *(unlabelled_gate_table(result) if 'unlabelled_gate' in result else [
            '| Arm | Change in log-odds, LOY-GFP | Change in log-odds, LOY-mCherry | Direction agreement |',
            '|---|---:|---:|---|',
            *[f"| {a['group']} | {a['gfp_delta_log_odds']:.3f} | {a['mcherry_delta_log_odds']:.3f} | {a['direction_agreement']} |"
              for a in result['dye_agreement']]]), '',
        'WT and WT C1 increase in both orientations; DDX3X arms show small decreases in both. '
        'SREBF1 and TADA2B have opposing orientation trends. Replicate relabeling does not resolve those disagreements. '
        + ('TADA2B uses the day-specific unlabelled-cell gate as its primary rule because an unlabelled TADA2B '
           'population crosses the standard mCherry cut between D0 and D4; the standard-gate columns keep the earlier '
           'TADA2B values. SREBF1 remains opposed under every gate. See `srebf1_tada2b_diagnosis/README.md`. '
           if result.get('unlabelled_gate', {}).get('primary_backgrounds') else
           'See `srebf1_tada2b_diagnosis/README.md` for the gate diagnosis. ' if 'unlabelled_gate' in result else '')
        + 'Inspect the strict-exclusion sensitivity and same-genotype controls in the HTML page.', '',
        '## Reproduce', '',
        '```powershell',
        '.\\codex\\flow_env\\python.exe -m spectral.a8_pipeline.growth_assay_followup --output <new-directory>',
        '```', '',
        '| Well | Analysis identity | Recorded filename | GFP % retained | Excluded % singlets |',
        '|---|---|---|---:|---:|',
    ]
    for w in wells:
        g = 'NA' if w['gfp_fraction'] is None else f"{100*w['gfp_fraction']:.2f}"
        e = 'NA' if w['excluded_broad_fraction'] is None else f"{100*w['excluded_broad_fraction']:.2f}"
        markdown.append(f"| {w['well']} | {w['analysis_sample']} | {w['original_filename']} | {g} | {e} |")
    (staging / 'report.md').write_text('\n'.join(markdown) + '\n', encoding='utf-8')
    provenance = {str(p.relative_to(REPO)): digest(p) for p in
                  (Path(__file__), TEMPLATE, Path(__file__).with_name('growth_assay_d0_report.py'),
                   Path(__file__).with_name('competition_report.py'), Path(__file__).with_name('gating.py'),
                   Path(__file__).with_name('fcs_io.py'), Path(__file__).with_name('longitudinal_report.py'),
                   Path(__file__).with_name('growth_assay_label_review.py'),
                   Path(__file__).with_name('growth_assay_unlabelled_gate.py'))}
    (staging / 'code_manifest.json').write_text(json.dumps(provenance, indent=2), encoding='utf-8')
    # Dropbox or antivirus can briefly hold a newly written directory on Windows.
    # Both paths are siblings; never replace an existing report when retrying.
    for attempt in range(4):
        if output.exists():
            raise FileExistsError(output)
        try:
            staging.rename(output)
            break
        except PermissionError:
            if attempt == 3:
                raise
            time.sleep(0.5 * (attempt + 1))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=DEFAULT_INPUT)
    parser.add_argument('--baseline', type=Path, default=DEFAULT_BASELINE)
    parser.add_argument('--output', type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument('--labels', type=Path, default=DEFAULT_LABELS)
    parser.add_argument('--unlabelled-gate', type=Path, default=DEFAULT_OUTPUT / 'unlabelled_gate_sensitivity.json',
                        help='Sensitivity file from growth_assay_unlabelled_gate; merged when present')
    parser.add_argument('--gating', type=Path, default=DEFAULT_GATING,
                        help='Approved backgrounds that use the unlabelled-cell gate as their primary rule')
    parser.add_argument('--cached', type=Path, help='Reuse measured counts and points from this results.json; recheck all source hashes')
    args = parser.parse_args()
    if args.output.exists() or args.output.with_name(args.output.name + '_staging').exists():
        raise FileExistsError('Output already exists; use --output with a new directory')
    baseline = json.loads(args.baseline.read_text(encoding='utf-8'))
    labels = json.loads(args.labels.read_text(encoding='utf-8'))
    if args.cached:
        result = json.loads(args.cached.read_text(encoding='utf-8'))
        if result['baseline_sha256'] != digest(args.baseline):
            raise ValueError('Cached D0 baseline hash changed')
        source = REPO / result['source']
        if result['source_manifest'] != inventory([p for p in source.rglob('*') if p.is_file()], REPO):
            raise ValueError('Cached source manifest does not match current raw files')
        if result['reference_sha256'] != digest(REPO / result['reference_file']):
            raise ValueError('D0 reference acquisition hash changed')
        result['measurement_cache'] = dict(path=args.cached.resolve().relative_to(REPO).as_posix(),
                                         sha256=digest(args.cached))
        result['cached_sources_reverified_utc'] = datetime.now(timezone.utc).isoformat()
    else:
        result = analyze(args.input, args.baseline)
    result = apply_labels_and_compare(result, baseline, labels)
    result['identity_mapping_path'] = args.labels.resolve().relative_to(REPO).as_posix()
    result['identity_mapping_sha256'] = digest(args.labels)
    attach_unlabelled_gate(result, args.labels, args.unlabelled_gate, args.gating)
    write_outputs(result, args.output)
    print(json.dumps(result['qc'], indent=2))
    print(args.output.resolve() / 'report.html')


if __name__ == '__main__':
    main()
