from copy import deepcopy
import json
import math
import unittest

import numpy as np

from spectral.a8_pipeline.growth_assay_d0_report import parse_sample
from spectral.a8_pipeline.growth_assay_timepoint import (
    DEFAULT_INPUT, DEFAULT_LABELS, GROUPS, ORIENTATIONS, assign_identities, biological, change, check_design,
    flag_wells, gate_changes, loy_of, reference, summarise, valley_cuts,
)
from spectral.a8_pipeline.growth_assay_followup import DEFAULT_BASELINE


def nominal(stem):
    s = parse_sample(stem)
    return dict(original_sample=stem, well=s.well, nominal_role=s.role, nominal_group=s.group,
                nominal_background=s.background, nominal_condition=s.condition, nominal_orientation=s.orientation,
                nominal_genotype=s.genotype, nominal_reporter=s.dye, nominal_replicate=s.replicate)


class IdentityTests(unittest.TestCase):
    def setUp(self):
        self.labels = json.loads(DEFAULT_LABELS.read_text(encoding='utf-8'))
        self.baseline = json.loads(DEFAULT_BASELINE.read_text(encoding='utf-8'))
        self.wells = [nominal(p.stem) for p in sorted(DEFAULT_INPUT.rglob('*.fcs'))]

    def test_d14_export_has_two_technical_replicates_per_biological_replicate(self):
        assign_identities(self.wells, self.labels)
        check_design(self.wells, self.baseline)
        by_well = {w['well']: w for w in self.wells}
        # Instrument counters (_002, _001) do not decide the technical replicate; column order does.
        self.assertEqual((by_well['B01']['technical_replicate'], by_well['B02']['technical_replicate']), (1, 2))
        self.assertEqual(by_well['B01']['analysis_sample'], 'B01-LOY_GFP_WT_C1_1_t1')
        self.assertEqual((by_well['C07']['genotype'], by_well['C07']['replicate']), ('CTRL', 1))
        self.assertTrue(all(w['identity_status'] == 'recorded' for w in self.wells))

    def test_missing_technical_replicate_is_rejected(self):
        with self.assertRaises(ValueError):
            assign_identities([w for w in self.wells if w['well'] != 'A02'], self.labels)

    def test_genotype_override_only_for_same_genotype_controls(self):
        labels = deepcopy(self.labels)
        labels['assignments'] = [dict(original_sample='A01-LOY_GFP_WT_1_001', replicate=1, genotype='CTRL')]
        with self.assertRaises(ValueError):
            assign_identities(self.wells, labels)

    def test_genotype_override_renames_both_technical_wells(self):
        labels = deepcopy(self.labels)
        labels['assignments'] = [dict(original_sample=s, genotype='CTRL') for s in
                                 ('C01-LOY_GFP_CHERRY_1_001', 'C02-LOY_GFP_CHERRY_1_003')]
        labels['assignments'] += [dict(original_sample=s, genotype='LOY') for s in
                                  ('C07-CTRL_GFP_CHERRY_1_003', 'C08-CTRL_GFP_CHERRY_1_004')]
        assign_identities(self.wells, labels)
        by_well = {w['well']: w for w in self.wells}
        self.assertEqual(by_well['C02']['analysis_sample'], 'C02-CTRL_GFP_CHERRY_1_t2')
        self.assertEqual(by_well['C07']['identity_status'], 'corrected')


class ArithmeticTests(unittest.TestCase):
    def test_change_has_no_pseudocount(self):
        self.assertIsNone(change(1.0, 0.5))
        self.assertIsNone(change(None, 0.5))
        self.assertAlmostEqual(change(0.75, 0.5), math.log(3))

    def test_loy_fraction_follows_orientation(self):
        self.assertEqual(loy_of(0.3, 'LOY_GFP'), 0.3)
        self.assertAlmostEqual(loy_of(0.3, 'LOY_mCherry'), 0.7)
        self.assertIsNone(loy_of(None, 'LOY_mCherry'))

    def test_technical_replicates_are_averaged_before_spread(self):
        rows = [dict(replicate=r, technical_replicate=t, well=f'A{2 * r + t - 2:02d}', loy_fraction=v, delta_log_odds=d)
                for r, t, v, d in ((1, 1, .40, .1), (1, 2, .46, .3), (2, 1, .50, .5), (2, 2, .50, .5))]
        bio = biological(rows, ('loy_fraction', 'delta_log_odds'), 'loy_fraction')
        self.assertEqual([b['wells'] for b in bio], ['A01;A02', 'A03;A04'])
        self.assertAlmostEqual(bio[0]['delta_log_odds'], .2)
        self.assertAlmostEqual(bio[0]['technical_gap_pp'], 6.0)
        self.assertEqual(bio[1]['technical_gap_pp'], 0)

    def test_reference_pools_batch_and_bounds_excluded_events(self):
        # Well counts are (primary singlets, retained, GFP among retained).
        counts = {'A01': (100, 80, 20), 'A02': (100, 90, 30)}
        gfp = reference(dict(wells='A01;A02', orientation='LOY_GFP'), counts)
        self.assertAlmostEqual(gfp['loy_fraction'], 50 / 170)
        self.assertAlmostEqual(gfp['excluded_as_ctrl'], 50 / 200)
        self.assertAlmostEqual(gfp['excluded_as_loy'], 80 / 200)
        mcherry = reference(dict(wells='A01;A02', orientation='LOY_mCherry'), counts)
        self.assertAlmostEqual(mcherry['loy_fraction'], 120 / 170)

    def test_gate_translation_is_detected(self):
        poly = lambda dx, dy: dict(parameters=['x', 'y'], vertices=[dict(x=a + dx, y=b + dy) for a, b in ((0, 0), (5, 0), (0, 5))])
        result = gate_changes(dict(A=poly(-4.25, 2), B=poly(0, 0)), dict(A=poly(0, 0), B=poly(0, 0)))
        self.assertEqual(result['A']['status'], 'translated')
        self.assertEqual(result['A']['shift'], [-4.25, 2.0])
        self.assertEqual(result['B']['status'], 'unchanged')

    def test_float32_jitter_is_still_a_translation_but_a_stretch_is_not(self):
        base = [(16356742.0, 25652368.0), (5719142.5, 7831994.5), (124298824.0, 2668817.0)]
        gate = lambda vs: dict(parameters=['x', 'y'], vertices=[dict(x=x, y=y) for x, y in vs])
        # The D14 Scatter polygon: every vertex moved by about (-4416136, +1938708), with float32 rounding.
        moved = [(11940606.0, 27591076.0), (1303007.25, 9770702.0), (119882696.0, 4607523.0)]
        self.assertEqual(gate_changes(dict(S=gate(moved)), dict(S=gate(base)))['S']['status'], 'translated')
        stretched = [(x * 1.01, y) for x, y in base]
        self.assertEqual(gate_changes(dict(S=gate(stretched)), dict(S=gate(base)))['S']['status'], 'reshaped')

    def test_reporter_negative_flag_is_relative_to_the_d0_mix(self):
        wells = [dict(well=w, excluded_broad_fraction=e, d0_excluded_fraction=d, singlet_yield=.8, total_events=10000)
                 for w, e, d in (('G01', .35, .25), ('F01', .65, .15), ('A01', .05, .01))]
        flag_wells(wells, [dict(wells='G01;F01;A01', technical_gap_pp=1.0)])
        self.assertEqual([w['qc_flags'] for w in wells], ['', 'reporter-negative 65% (D0 mix 15%)', ''])

    def test_valley_needs_both_days(self):
        rng = np.random.default_rng(0)
        two_peaks = np.column_stack([np.zeros(4000), np.r_[rng.normal(6.2, .2, 2000), rng.normal(8.3, .2, 2000)]])
        one_peak = np.column_stack([np.zeros(2000), rng.normal(8.3, .2, 2000)])
        valleys, cuts = valley_cuts({'A': {'D0': two_peaks, 'D14': two_peaks}, 'B': {'D0': one_peak, 'D14': two_peaks}}, 6.5)
        self.assertTrue(valleys['A']['detected'])
        self.assertTrue(6.5 < cuts['A']['D14'] < 8)
        self.assertFalse(valleys['B']['detected'])
        self.assertEqual(cuts['B'], {'D0': 6.5, 'D14': 6.5})


class SummaryTests(unittest.TestCase):
    def test_spread_is_across_biological_replicates_and_contrasts_use_arm_means(self):
        wells, before = [], []
        for g, group in enumerate(GROUPS):
            for o, orientation in enumerate(ORIENTATIONS):
                before.append(dict(group=group, orientation=orientation, loy_fraction_mean=.5, delta_log_odds_mean=.1 * g,
                                   delta_log_odds_standard_gate_mean=.1 * g))
                for r in (1, 2, 3):
                    for t in (1, 2):
                        d = .1 * g + .01 * r + (.02 if t == 2 else 0)
                        wells.append(dict(role='mixture', group=group, orientation=orientation, replicate=r,
                                          technical_replicate=t, well=f'{"ABCDEFG"[g]}{6 * o + 2 * r + t - 2:02d}',
                                          primary_gate='standard', d0_loy_fraction=.5, d0_wells='X',
                                          d0_batch='DDX3X' if group.startswith('DDX3X') else group, loy_fraction=.5 + d / 10,
                                          delta_log_odds=d, delta_log_odds_standard_gate=d, delta_log_odds_unlabelled_gate=d,
                                          delta_log_odds_strict=d, excluded_broad_fraction=.05,
                                          delta_log_odds_excluded_as_ctrl=d, delta_log_odds_excluded_as_loy=d))
        for genotype in ('CTRL', 'LOY'):
            for r in (1, 2, 3):
                for t in (1, 2):
                    wells.append(dict(role='dye_neutral', genotype=genotype, replicate=r, technical_replicate=t,
                                      well=f'C{r}{t}', d0_gfp_fraction=.5, d0_wells='C', gfp_fraction=.55,
                                      gfp_delta_log_odds=change(.55, .5)))
        previous = dict(timepoint='D4', orientation_summaries=before,
                        dye_controls=[dict(genotype=g, d4_gfp_fraction_mean=.5, mean_delta_log_odds=0) for g in ('CTRL', 'LOY')])
        out = summarise(wells, previous)
        wt = next(s for s in out['orientation_summaries'] if (s['group'], s['orientation']) == ('WT', 'LOY_GFP'))
        self.assertEqual(wt['delta_log_odds_n'], 3)
        self.assertAlmostEqual(wt['delta_log_odds_mean'], .03)
        self.assertAlmostEqual(wt['delta_log_odds_sd'], .01)
        self.assertAlmostEqual(wt['technical_gap_pp_max'], .2)
        wt = next(c for c in out['replicate_set_comparisons'] if c['second_set'] == 'WT_C1' and c['orientation'] == 'LOY_GFP')
        self.assertAlmostEqual(wt['difference_since_d0'], .1)
        self.assertAlmostEqual(wt['difference_since_previous'], 0)
        self.assertFalse(wt['shared_d0_mix'])
        ddx = next(c for c in out['replicate_set_comparisons'] if c['second_set'] == 'DDX3X_C1')
        self.assertTrue(ddx['shared_d0_mix'])
        # WT pooled over both sets: six biological replicates (.02-.04 and .12-.14); mean of set means .08.
        pooled = next(b for b in out['background_summaries'] if (b['background'], b['orientation']) == ('WT', 'LOY_GFP'))
        self.assertEqual((pooled['sets'], pooled['delta_log_odds_n']), ('WT;WT_C1', 6))
        self.assertAlmostEqual(pooled['mean_of_set_means'], .08)
        self.assertAlmostEqual(pooled['delta_log_odds_mean'], .08)
        single = next(b for b in out['background_summaries'] if b['background'] == 'SREBF1')
        self.assertEqual(single['delta_log_odds_n'], 3)
        ctrl = next(c for c in out['dye_controls'] if c['genotype'] == 'CTRL')
        self.assertAlmostEqual(ctrl['previous_gfp_fraction_mean'], .5)
        self.assertAlmostEqual(ctrl['change_since_previous_log_odds'], change(.55, .5))


if __name__ == '__main__':
    unittest.main()
