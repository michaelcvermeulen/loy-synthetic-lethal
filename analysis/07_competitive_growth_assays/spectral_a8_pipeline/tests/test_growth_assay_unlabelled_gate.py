import json
import tempfile
import unittest
from copy import deepcopy
from pathlib import Path

import numpy as np

from spectral.a8_pipeline.growth_assay_followup import (DEFAULT_BASELINE, DEFAULT_GATING, DEFAULT_LABELS, DEFAULT_OUTPUT,
                                                        REPO, attach_unlabelled_gate)
from spectral.a8_pipeline.growth_assay_label_review import apply_labels_and_compare
from spectral.a8_pipeline.growth_assay_unlabelled_gate import find_valley, gated_counts

SENSITIVITY = DEFAULT_OUTPUT / 'unlabelled_gate_sensitivity.json'


def reviewed_from_cache():
    """Label-reviewed result rebuilt from the raw measurement cache, before any gating is applied."""
    report = json.loads((DEFAULT_OUTPUT / 'results.json').read_text(encoding='utf-8'))
    imported = json.loads((REPO / report['measurement_cache']['path']).read_text(encoding='utf-8'))
    baseline = json.loads(DEFAULT_BASELINE.read_text(encoding='utf-8'))
    labels = json.loads(DEFAULT_LABELS.read_text(encoding='utf-8'))
    return apply_labels_and_compare(imported, baseline, labels)


class ValleyTests(unittest.TestCase):
    def test_bimodal_valley_lies_between_peaks(self):
        rng = np.random.default_rng(1)
        m = np.concatenate([rng.normal(6.1, .35, 8000), rng.normal(8.2, .3, 10000)])
        v = find_valley(m)
        self.assertAlmostEqual(v['unlabelled_peak'], 6.1, delta=.15)
        self.assertAlmostEqual(v['mcherry_peak'], 8.2, delta=.15)
        self.assertTrue(6.8 < v['valley'] < 7.6)
        self.assertLess(v['valley_depth'], 1)

    def test_unimodal_distribution_has_no_valley(self):
        rng = np.random.default_rng(2)
        with self.assertRaisesRegex(ValueError, 'No separating valley'):
            find_valley(rng.normal(8.2, .5, 20000))

    def test_gate_counts_follow_cut_and_classifier(self):
        f = np.array([[9., 1.], [1., 9.], [1., 7.], [1., 6.]])
        weights, intercept = [1., 0.], -5.
        self.assertEqual(gated_counts(f, weights, intercept, 6, 6.5)['retained'], 3)
        c = gated_counts(f, weights, intercept, 6, 7.5)
        self.assertEqual((c['retained'], c['gfp_count']), (2, 1))
        self.assertAlmostEqual(c['excluded_fraction'], .5)


@unittest.skipUnless(SENSITIVITY.exists(), 'sensitivity output not generated')
class MergeTests(unittest.TestCase):
    def setUp(self):
        self.result = reviewed_from_cache()
        self.sensitivity = json.loads(SENSITIVITY.read_text(encoding='utf-8'))

    def test_sensitivity_only_merge_keeps_primary_fields(self):
        before = deepcopy(self.result)
        attach_unlabelled_gate(self.result, DEFAULT_LABELS, SENSITIVITY)
        for old, new in zip(before['dye_agreement'], self.result['dye_agreement']):
            self.assertAlmostEqual(old['gfp_delta_log_odds'], new['gfp_delta_log_odds'])
            self.assertAlmostEqual(old['gfp_delta_log_odds'], new['gfp_delta_log_odds_standard_gate'])
        for old, new in zip(before['wells'], self.result['wells']):
            self.assertEqual(old.get('delta_log_odds'), new.get('delta_log_odds'))
            self.assertEqual(new['primary_gate'], 'standard')

    def test_approved_gating_makes_tada2b_primary_and_keeps_standard_values(self):
        before = {w['well']: w for w in deepcopy(self.result)['wells']}
        attach_unlabelled_gate(self.result, DEFAULT_LABELS, SENSITIVITY, DEFAULT_GATING)
        gate = {w['well']: w for w in self.sensitivity['wells']}
        for w in self.result['wells']:
            old = before[w['well']]
            if w['role'] != 'mixture':
                continue
            self.assertEqual(w['delta_log_odds_standard_gate'], old['delta_log_odds'])
            if w['background'] == 'TADA2B':
                self.assertEqual(w['primary_gate'], 'unlabelled_cell_gate')
                self.assertEqual(w['delta_log_odds'], gate[w['well']]['delta_log_odds'])
                self.assertEqual(w['gfp_count'], gate[w['well']]['gfp_count'])
            else:
                self.assertEqual(w['primary_gate'], 'standard')
                self.assertEqual(w['delta_log_odds'], old['delta_log_odds'])
        tada2b = next(a for a in self.result['dye_agreement'] if a['group'] == 'TADA2B')
        self.assertAlmostEqual(tada2b['gfp_delta_log_odds'], tada2b['gfp_delta_log_odds_unlabelled_gate'])
        self.assertLess(tada2b['gfp_delta_log_odds_standard_gate'], 0)
        baseline = next(b for b in self.result['d0_baselines'] if b['group'] == 'TADA2B' and b['orientation'] == 'LOY_GFP')
        self.assertNotEqual(baseline['d0_loy_fraction'], baseline['d0_loy_fraction_standard_gate'])

    def test_gating_requires_a_detected_valley(self):
        gating = json.loads(DEFAULT_GATING.read_text(encoding='utf-8'))
        gating['unlabelled_gate_primary_backgrounds'] = ['SREBF1']
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / 'gating.json'
            path.write_text(json.dumps(gating), encoding='utf-8')
            with self.assertRaisesRegex(ValueError, 'separating valley'):
                attach_unlabelled_gate(self.result, DEFAULT_LABELS, SENSITIVITY, path)

    def test_gating_cannot_be_applied_twice(self):
        attach_unlabelled_gate(self.result, DEFAULT_LABELS, SENSITIVITY, DEFAULT_GATING)
        with self.assertRaisesRegex(ValueError, 'already applied'):
            attach_unlabelled_gate(self.result, DEFAULT_LABELS, SENSITIVITY, DEFAULT_GATING)

    def test_stale_sensitivity_file_is_rejected(self):
        data = deepcopy(self.sensitivity)
        data['identity_mapping_sha256'] = '0' * 64
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / 'stale.json'
            path.write_text(json.dumps(data), encoding='utf-8')
            with self.assertRaisesRegex(ValueError, 'rerun it'):
                attach_unlabelled_gate(self.result, DEFAULT_LABELS, path)

    def test_missing_file_leaves_result_unchanged(self):
        before = deepcopy(self.result)
        attach_unlabelled_gate(self.result, DEFAULT_LABELS, DEFAULT_OUTPUT / 'does_not_exist.json')
        self.assertEqual(self.result, before)


if __name__ == '__main__':
    unittest.main()
