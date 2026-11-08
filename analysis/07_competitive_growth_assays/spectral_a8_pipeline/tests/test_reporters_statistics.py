from __future__ import annotations

import unittest

import numpy as np

from spectral.a8_pipeline.reporters import (
    classify_reporters,
    equal_posterior_threshold,
    fit_two_gaussian,
)
from spectral.a8_pipeline.statistics import competition_metrics, wilson_interval


class ReporterAndStatisticsTests(unittest.TestCase):
    def test_threshold_lies_between_separated_modes(self) -> None:
        rng = np.random.default_rng(7)
        values = np.r_[rng.normal(-2, 0.3, 2000), rng.normal(3, 0.4, 2000)]
        model = fit_two_gaussian(values, seed=20260828, n_init=8)
        threshold = equal_posterior_threshold(model)
        self.assertLess(model.means[0], threshold)
        self.assertLess(threshold, model.means[1])

    def test_invalid_overlapping_model_is_rejected(self) -> None:
        rng = np.random.default_rng(9)
        values = rng.normal(0, 1, 4000)
        model = fit_two_gaussian(values, seed=20260828, n_init=8)
        with self.assertRaisesRegex(ValueError, "insufficient component separation"):
            equal_posterior_threshold(model, minimum_separation=2.0)

    def test_quadrants_and_dye_orientation(self) -> None:
        classes = classify_reporters(
            np.array([10, 10, 0, 0], dtype=float),
            np.array([0, 10, 10, 0], dtype=float),
            gfp_threshold=5,
            mcherry_threshold=5,
        )
        self.assertEqual(
            classes.tolist(),
            ["GFP-only", "Double-positive", "mCherry-only", "Double-negative"],
        )
        f11 = competition_metrics(classes, loy_reporter="GFP")
        g11 = competition_metrics(classes, loy_reporter="mCherry")
        self.assertEqual((f11.loy, f11.ctrl), (1, 1))
        self.assertEqual((g11.loy, g11.ctrl), (1, 1))

    def test_wilson_interval_has_expected_counting_bounds(self) -> None:
        low, high = wilson_interval(50, 100)
        self.assertAlmostEqual(low, 0.4038, places=4)
        self.assertAlmostEqual(high, 0.5962, places=4)


if __name__ == "__main__":
    unittest.main()
