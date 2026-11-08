from pathlib import Path
from copy import deepcopy
import json
import unittest

import numpy as np

from spectral.a8_pipeline.growth_assay_followup import (
    DEFAULT_BASELINE, DEFAULT_LABELS, DEFAULT_OUTPUT, measure_reporters, provisional_identity,
)
from spectral.a8_pipeline.growth_assay_label_review import apply_labels_and_compare, log_odds_change


class FollowupTests(unittest.TestCase):
    def setUp(self):
        self.baseline = json.loads(DEFAULT_BASELINE.read_text(encoding='utf-8'))

    def test_d0_matching_uses_identity_not_physical_well(self):
        # A01 was LOY-mCherry at D0, but is named LOY-GFP at D4.
        row = provisional_identity('A01-LOY_GFP_WT_1_001', self.baseline)
        self.assertEqual(row['d0_wells_provisional'], 'A04;A05;A06')
        self.assertEqual(row['identity_status'], 'pending_correction')

    def test_preserves_out_of_order_replicates_in_recorded_names(self):
        row = provisional_identity('B04-LOY_CHERRY_WT_C1_2_002', self.baseline)
        self.assertEqual(row['nominal_replicate'], 2)
        self.assertEqual(row['nominal_group'], 'WT_C1')

    def test_dye_neutral_does_not_get_a_loy_fraction_or_baseline(self):
        row = provisional_identity('C05-LOY_GFP_CHERRY_3_002', self.baseline)
        self.assertEqual(row['nominal_genotype'], 'LOY')
        self.assertEqual(row['nominal_replicate'], 3)
        self.assertNotIn('d0_loy_fraction_provisional', row)

    def test_shared_ddx3x_baseline_and_separate_wt_baselines(self):
        def baseline(stem):
            return provisional_identity(stem, self.baseline)['d0_loy_fraction_provisional']
        self.assertEqual(baseline('D01-LOY_GFP_DDX3X_1_002'),
                         baseline('E01-LOY_GFP_DDX3X_C1_1_001'))
        self.assertNotEqual(baseline('A01-LOY_GFP_WT_1_001'),
                            baseline('B01-LOY_GFP_WT_C1_1_002'))

    def test_frozen_model_exclusion_counts_and_denominators(self):
        features = np.array([[10., -4.], [10., -4.], [2., 10.], [0., -4.], [0., 5.]])
        row = measure_reporters(features, self.baseline)
        self.assertEqual((row['primary_singlets'], row['reporter_positive_singlets']), (5, 3))
        self.assertEqual((row['gfp_count'], row['mcherry_count']), (2, 1))
        self.assertAlmostEqual(row['gfp_fraction'], 2/3)
        self.assertAlmostEqual(row['excluded_broad_fraction'], 2/5)

    def test_exact_cutoff_is_retained(self):
        row = measure_reporters(np.array([[6., 0.], [0., 6.5], [5.9, 6.4]]), self.baseline)
        self.assertEqual(row['reporter_positive_singlets'], 2)

    def test_empty_retained_population_is_missing_not_zero_percent(self):
        row = measure_reporters(np.array([[0., 0.], [0., 5.]]), self.baseline)
        self.assertIsNone(row['gfp_fraction'])
        self.assertIsNone(row['mcherry_fraction'])
        self.assertEqual(row['reporter_positive_singlets'], 0)

    def test_nonfinite_reporters_fail_instead_of_silent_classification(self):
        with self.assertRaises(ValueError):
            measure_reporters(np.array([[np.nan, 10.]]), self.baseline)


class ApprovedLabelTests(unittest.TestCase):
    def setUp(self):
        self.baseline = json.loads(DEFAULT_BASELINE.read_text(encoding='utf-8'))
        # The report's results.json carries applied gating; label review starts from the raw measurement cache.
        report = json.loads((DEFAULT_OUTPUT / 'results.json').read_text(encoding='utf-8'))
        cache = DEFAULT_OUTPUT.parents[1] / report['measurement_cache']['path']
        self.imported = json.loads(cache.read_text(encoding='utf-8'))
        self.labels = json.loads(DEFAULT_LABELS.read_text(encoding='utf-8'))

    def test_approved_triplets_and_original_measurements_are_preserved(self):
        before = deepcopy(self.imported)
        reviewed = apply_labels_and_compare(self.imported, self.baseline, self.labels)
        self.assertEqual(self.imported, before)
        by_well = {w['well']: w for w in reviewed['wells']}
        for prefix, expected in [('B0', [2, 1, 3, 1, 2, 3]), ('C0', [1, 3, 2, 1, 2, 3])]:
            self.assertEqual([by_well[prefix + str(i)]['replicate'] for i in range(1, 7)], expected)
        self.assertEqual(reviewed['qc']['corrected_replicate_count'], 8)
        self.assertEqual([by_well['C0' + str(i)]['genotype'] for i in range(1, 7)], ['LOY'] * 3 + ['CTRL'] * 3)
        self.assertEqual([by_well['C0' + str(i)]['nominal_genotype'] for i in range(1, 7)], ['CTRL'] * 3 + ['LOY'] * 3)
        self.assertEqual(by_well['C02']['analysis_sample'], 'C02-LOY_GFP_CHERRY_3')
        self.assertEqual(reviewed['qc']['corrected_genotype_count'], 6)
        self.assertFalse(any(w['genotype_changed'] for w in reviewed['wells'] if w['role'] == 'mixture'))
        for original, row in zip(before['wells'], reviewed['wells']):
            for key in ('source_file', 'original_filename', 'original_sample', 'nominal_replicate',
                        'gfp_count', 'mcherry_count', 'gfp_fraction', 'primary_singlets', 'reporter_positive_singlets'):
                self.assertEqual(row[key], original[key])
            self.assertEqual(row['orientation'], original['nominal_orientation'])
        self.assertEqual(reviewed['plot_points'], before['plot_points'])
        self.assertEqual(reviewed['source_manifest'], before['source_manifest'])

    def test_replicate_collision_is_rejected(self):
        self.labels['assignments'][0]['replicate'] = 1
        with self.assertRaisesRegex(ValueError, 'duplicate replicates'):
            apply_labels_and_compare(self.imported, self.baseline, self.labels)

    def test_wrong_export_or_unknown_recording_is_rejected(self):
        for bad in ('source', 'recording'):
            labels = deepcopy(self.labels)
            if bad == 'source':
                labels['source'] = 'some_other_export'
            else:
                labels['assignments'][0]['original_sample'] = 'unknown'
            with self.assertRaises(ValueError):
                apply_labels_and_compare(self.imported, self.baseline, labels)

    def test_genotype_correction_is_limited_to_same_genotype_controls(self):
        for sample, genotype in (('B01-LOY_GFP_WT_C1_1_002', 'CTRL'), ('C01-CTRL_GFP_CHERRY_1_001', 'WT')):
            labels = deepcopy(self.labels)
            next(a for a in labels['assignments'] if a['original_sample'] == sample)['genotype'] = genotype
            with self.assertRaisesRegex(ValueError, 'Genotype corrections'):
                apply_labels_and_compare(self.imported, self.baseline, labels)

    def test_one_sided_genotype_swap_is_rejected(self):
        labels = deepcopy(self.labels)
        for a in labels['assignments']:
            if a['original_sample'].startswith(('C04', 'C05', 'C06')):
                del a['genotype']
        with self.assertRaisesRegex(ValueError, 'duplicate replicates'):
            apply_labels_and_compare(self.imported, self.baseline, labels)

    def test_shared_baseline_and_equal_weighted_well_effects(self):
        reviewed = apply_labels_and_compare(self.imported, self.baseline, self.labels)
        ws = [w for w in reviewed['wells'] if w['role'] == 'mixture' and w['group'] == 'WT_C1' and w['orientation'] == 'LOY_GFP']
        b = next(b for b in self.baseline['d0_baselines'] if b['group'] == 'WT_C1' and b['orientation'] == 'LOY_GFP')
        self.assertTrue(all(w['d0_loy_fraction'] == b['d0_loy_fraction'] for w in ws))
        summary = next(s for s in reviewed['orientation_summaries'] if s['group'] == 'WT_C1' and s['orientation'] == 'LOY_GFP')
        expected = np.mean([np.log(w['gfp_count'] / w['mcherry_count']) - np.log(b['d0_loy_fraction'] / (1-b['d0_loy_fraction'])) for w in ws])
        self.assertAlmostEqual(summary['delta_log_odds_mean'], expected)
        self.assertEqual(summary['wells_by_replicate'], 'B02;B01;B03')

    def test_no_change_doubling_and_boundary_log_odds(self):
        self.assertEqual(log_odds_change(0.25, 0.25), 0)
        self.assertAlmostEqual(log_odds_change(0.4, 0.25), np.log(2))
        self.assertIsNone(log_odds_change(0, 0.25))
        self.assertIsNone(log_odds_change(None, 0.25))


if __name__ == '__main__':
    unittest.main()
