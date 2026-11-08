"""Exploratory same-genotype sensitivity checks; do not alter primary results."""
from __future__ import annotations

from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
import zipfile

import numpy as np

from .competition_report import load_cef_gates
from .fcs_io import load_fcs
from .gating import points_in_polygon
from .saturation import DETECTOR_NAMES, MASK_COLUMNS, decode_mask_matrix, decode_saturation

ROOT = Path(__file__).resolve().parents[2]
DIRECTORY = ROOT / 'analysis/spectral_growth_assay_D4_20261020/same_genotype'


def main():
    path = DIRECTORY / 'data.json'
    data = json.loads(path.read_text(encoding='utf-8'))
    approved_d4 = json.loads((DIRECTORY.parent / 'results.json').read_text(encoding='utf-8'))
    approved_d0 = json.loads((ROOT / approved_d4['baseline_path']).read_text(encoding='utf-8'))
    approved_sources = {s['path']: s['sha256'] for d in (approved_d0, approved_d4) for s in d['source_manifest']}
    hashes = dict(data['source_hashes'])
    hashes[path.relative_to(ROOT).as_posix()] = sha256(path.read_bytes()).hexdigest()
    for rel, expected in hashes.items():
        if sha256((ROOT / rel).read_bytes()).hexdigest() != expected:
            raise ValueError(f'Input changed: {rel}')
    weights = np.array(data['classifier']['weights'])
    intercept = data['classifier']['intercept']
    gates = {day: load_cef_gates(next((ROOT / source).rglob('*.cef'))) for day, source in
             [(0, '20261016_171604_Export'), (4, '20261020_173018_Export')]}
    same_gates = {name: gates[0][name].parameters == gates[4][name].parameters and
                 np.array_equal(gates[0][name].vertices, gates[4][name].vertices) for name in gates[0]}
    wells = []
    for sample in data['samples']:
        recording = (ROOT / sample['source_file']).with_suffix('.zip')
        recording_rel = recording.relative_to(ROOT).as_posix()
        recording_hash = sha256(recording.read_bytes()).hexdigest()
        if approved_sources[recording_rel] != recording_hash:
            raise ValueError(f'Recording configuration hash changed: {recording_rel}')
        hashes[recording_rel] = recording_hash
        with zipfile.ZipFile(recording) as archive:
            configuration = json.loads(archive.read('recording-data/configuration.json'))
        detector_names = tuple(s['detectorName'] for s in configuration['detectorSettings'])
        if detector_names != DETECTOR_NAMES:
            raise ValueError('Recording detector order differs from the validated saturation decoder')
        fcs = load_fcs(ROOT / sample['source_file'])
        events = fcs.events
        flags = decode_saturation(events)
        if not np.array_equal(flags['any_saturated'].to_numpy(), events['Saturated'].to_numpy() != 0):
            raise ValueError('Decoded saturation flags disagree with the instrument flag')
        detector_flags = decode_mask_matrix(events[list(MASK_COLUMNS)].to_numpy())
        if detector_flags[:, len(detector_names):].any():
            raise ValueError('Unmapped saturation bits are active')
        features = np.arcsinh(events[['eGFP-A', 'mCherry-A']].to_numpy(float)/1000)
        distance = (features @ weights + intercept)/weights[0]
        positive = ~((features[:, 0] < 6) & (features[:, 1] < 6.5))
        def count(mask, included=positive, boundary_shift=0):
            selected = mask & included
            n = int(selected.sum())
            green = int((selected & (distance > boundary_shift)).sum())
            return dict(gfp=green, n=n, fraction=green/n if n else None)
        keep = np.ones(len(events), bool)
        variants = {}
        for name in ('Scatter', 'SSC Singlets', 'FSC Singlets'):
            gate = gates[sample['day']][name]
            keep &= points_in_polygon(events[list(gate.parameters)].to_numpy(float), gate.vertices)
            variants[name] = count(keep)
        variants['exclude_saturated'] = count(keep & (events['Saturated'].to_numpy() == 0))
        variants['exclude_spectral_saturation'] = count(keep & ~flags['spectral_saturated'].to_numpy())
        variants['strict_negatives'] = count(keep, ~((features[:, 0] < 6) & (features[:, 1] < 2)))
        variants['no_negative_exclusion'] = count(keep, np.ones(len(events), bool))
        for shift in (-.5, .5):
            variants[f'boundary_shift_{shift:+.1f}'] = count(keep, boundary_shift=shift)
        for gfp in (5.5, 6, 6.5):
            for mch in (6, 6.5, 7):
                variants[f'negative_gate_{gfp}_{mch}'] = count(keep, ~((features[:, 0] < gfp) & (features[:, 1] < mch)))
        quarters = []
        for indices in np.array_split(np.flatnonzero(keep), 4):
            mask = np.zeros(len(events), bool)
            mask[indices] = True
            quarters.append(count(mask))
        metrics = dict(near_boundary_retained_fraction=float((np.abs(distance[keep & positive]) < .5).mean()),
                       saturated_fraction=float((events.loc[keep, 'Saturated'] > 0).mean()),
                       spectral_saturated_fraction=float(flags.loc[keep, 'spectral_saturated'].mean()),
                       imgb3_saturated_fraction=float(flags.loc[keep, 'imgb3_saturated'].mean()),
                       imgb3_only_saturated_fraction=float(flags.loc[keep, 'imgb3_only_saturated'].mean()),
                       scatter_lightloss_saturated_fraction=float(flags.loc[keep, 'scatter_lightloss_saturated'].mean()),
                       gfp_median_asinh=float(np.median(features[keep & positive & (distance > 0), 0])),
                       mcherry_median_asinh=float(np.median(features[keep & positive & (distance <= 0), 1])))
        wells.append(dict(day=sample['day'], genotype=sample['genotype'], replicate=sample['replicate'],
                          well=sample['well'], source_file=sample['source_file'], variants=variants,
                          acquisition_order_quarters=quarters, metrics=metrics,
                          saturation_decoder_verified=True,
                          saturation_by_detector={name: int(detector_flags[keep, i].sum())
                                                  for i, name in enumerate(detector_names) if detector_flags[keep, i].any()}))
    summaries = []
    for genotype in ('CTRL', 'LOY'):
        for variant in wells[0]['variants']:
            initial = [s['variants'][variant] for s in wells if s['day'] == 0 and s['genotype'] == genotype]
            final = [s['variants'][variant] for s in wells if s['day'] == 4 and s['genotype'] == genotype]
            p0 = sum(s['gfp'] for s in initial)/sum(s['n'] for s in initial)
            p4 = float(np.mean([s['fraction'] for s in final]))
            summaries.append(dict(genotype=genotype, variant=variant, d0_pooled_gfp=p0, d4_mean_gfp=p4,
                                  change_pp=100*(p4-p0)))
    for rel, expected in hashes.items():
        if sha256((ROOT / rel).read_bytes()).hexdigest() != expected:
            raise ValueError(f'Input changed during checks: {rel}')
    output = dict(generated_utc=datetime.now(timezone.utc).isoformat(), input_hashes=hashes,
                  script_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
                  purpose='Exploratory sensitivities only; approved primary results unchanged.',
                  identical_gate_parameters_and_vertices=same_gates, summaries=summaries, wells=wells,
                  notes=['Acquisition quarters use event order, not equal elapsed-time bins.',
                         'Boundary shifts move the classifier horizontally in transformed GFP units at both days.',
                         'Earlier gate stages can contain doublets; they are sensitivity checks, not alternative validated analyses.',
                         'Saturated is an any-detector flag. SaturatedChannels masks identify the detectors.',
                         'Detector order was verified against every recording configuration; decoded masks match the Saturated flag event-for-event.'])
    (DIRECTORY / 'technical_checks.json').write_text(json.dumps(output, indent=2), encoding='utf-8')
    print('Identical polygons:', same_gates)
    for genotype in ('CTRL', 'LOY'):
        print(genotype)
        for s in summaries:
            if s['genotype'] == genotype and not s['variant'].startswith('negative_gate'):
                print(s['variant'], 'D0', round(100*s['d0_pooled_gfp'], 2), 'D4', round(100*s['d4_mean_gfp'], 2), 'change pp', round(s['change_pp'], 2))
        grid = [s['change_pp'] for s in summaries if s['genotype'] == genotype and s['variant'].startswith('negative_gate')]
        print('Negative threshold grid change range:', min(grid), max(grid))
        for s in wells:
            if s['genotype'] == genotype:
                print('D'+str(s['day']), s['well'], 'GFP% across acquisition quarters', [round(100*q['fraction'], 1) for q in s['acquisition_order_quarters']],
                      'saturated %', round(100*s['metrics']['saturated_fraction'], 1))
        for day in (0, 4):
            selected = [s for s in wells if s['genotype'] == genotype and s['day'] == day]
            print('D'+str(day), 'mean detector saturation %',
                  {key: round(100*float(np.mean([s['metrics'][key] for s in selected])), 4)
                   for key in ('saturated_fraction', 'spectral_saturated_fraction', 'imgb3_saturated_fraction',
                               'imgb3_only_saturated_fraction', 'scatter_lightloss_saturated_fraction')})


if __name__ == '__main__':
    main()
