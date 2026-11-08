from __future__ import annotations

import math
import unittest
from pathlib import Path

import numpy as np

from spectral.a8_pipeline.competition_report import LDAModel
from spectral.a8_pipeline.growth_assay_d0_report import (
    LOY_KEYS,
    build_baselines,
    build_design,
    measure_well,
    negative_mask,
    parse_sample,
    summarise_conditions,
    unstained_check,
)

# GFP class whenever eGFP > 6 (weights on eGFP only).
GFP_AXIS = LDAModel(np.array([1.0, 0.0]), -6.0)

EXPORT = (
    Path(__file__).resolve().parents[3]
    / "20261016_171604_Export"
    / "comp_growth_assay_D0_OCT_15"
    / "Plate_004"
)


class ParseSampleTests(unittest.TestCase):
    def test_mixture_dye_names_the_loy_reporter(self) -> None:
        """Swapping the orientation would invert every LOY fraction."""
        sample = parse_sample("D01-LOY_CHERRY_DDX3X_1_003")
        self.assertEqual(
            (sample.role, sample.background, sample.orientation, sample.replicate),
            ("mixture", "DDX3X", "LOY_mCherry", 1),
        )

    def test_c1_is_a_separate_condition_not_a_replicate(self) -> None:
        """C1 wells must never be pooled with the standard WT wells."""
        sample = parse_sample("B04-LOY_GFP_WT_C1_1_002")
        self.assertEqual((sample.group, sample.replicate), ("WT_C1", 1))

    def test_rerecord_counter_is_not_read_as_replicate(self) -> None:
        """A01 has no counter; A02 has _001 — both are replicates 1 and 2."""
        self.assertEqual(parse_sample("A01-LOY_CHERRY_WT_1").replicate, 1)
        self.assertEqual(parse_sample("A02-LOY_CHERRY_WT_2_001").replicate, 2)

    def test_solo_and_dye_neutral_wells(self) -> None:
        solo = parse_sample("H09-CTRL_CHERRY_SOLO_TADA2B_004")
        neutral = parse_sample("C04-CTRL_GFP_CHERRY_1_001")
        self.assertEqual((solo.role, solo.genotype, solo.dye), ("solo", "CTRL", "mCherry"))
        self.assertEqual((neutral.role, neutral.genotype), ("dye_neutral", "CTRL"))

    def test_unknown_names_are_rejected(self) -> None:
        with self.assertRaises(ValueError):
            parse_sample("A01-LOY_CHERRY_WT_4")

    @unittest.skipUnless(EXPORT.exists(), "export not present")
    def test_real_plate_is_complete(self) -> None:
        samples = build_design(path.stem for path in EXPORT.glob("*.fcs"))
        roles = [sample.role for sample in samples]
        self.assertEqual(
            (len(samples), roles.count("mixture"), roles.count("solo"), roles.count("dye_neutral")),
            (58, 36, 16, 6),
        )


class SummaryTests(unittest.TestCase):
    def test_balanced_estimate_is_logit_mean_of_orientations(self) -> None:
        wells = []
        for orientation, values in (("LOY_GFP", (0.30, 0.30)), ("LOY_mCherry", (0.50, 0.50))):
            for index, value in enumerate(values):
                wells.append(
                    {
                        "role": "mixture", "group": "WT", "background": "WT", "condition": "standard",
                        "orientation": orientation, "well": f"{orientation}{index}", "primary_singlets": 100,
                        "reporter_positive_singlets": 100, "negative_broad_fraction": 0.0,
                        **{key: value for key in LOY_KEYS},
                    }
                )
        _orientation, balanced = summarise_conditions(wells)
        expected = 1 / (1 + math.exp(-0.5 * math.log(0.3 / 0.7)))
        self.assertAlmostEqual(balanced[0]["loy_fraction_balanced"], expected)
        self.assertAlmostEqual(balanced[0]["loy_fraction_dye_swap_pp"], 20.0)

    def test_ddx3x_arms_share_a_baseline_but_wt_arms_do_not(self) -> None:
        """DDX3X/DDX3X C1 were mixed once; WT/WT C1 were two separate batches."""
        wells = []
        for group, loy in (("WT", 30), ("WT_C1", 40), ("DDX3X", 50), ("DDX3X_C1", 52)):
            for orientation in ("LOY_GFP", "LOY_mCherry"):
                wells.append(
                    {"role": "mixture", "group": group, "orientation": orientation, "well": group,
                     "loy_count": loy, "reporter_positive_singlets": 100, "loy_fraction": loy / 100}
                )
        baselines = {
            (row["group"], row["orientation"]): row["d0_loy_fraction"] for row in build_baselines(wells)
        }
        self.assertAlmostEqual(baselines[("DDX3X", "LOY_GFP")], 0.51)
        self.assertAlmostEqual(baselines[("DDX3X_C1", "LOY_GFP")], 0.51)
        self.assertAlmostEqual(baselines[("WT", "LOY_GFP")], 0.30)
        self.assertAlmostEqual(baselines[("WT_C1", "LOY_GFP")], 0.40)

    def test_negative_rules_nest(self) -> None:
        """Every strict negative must also be a broad negative."""
        points = np.array([[0.0, -5.0], [0.0, 5.0], [9.0, -5.0], [3.0, 9.0]])
        strict, broad = negative_mask(points, "strict"), negative_mask(points, "broad")
        self.assertEqual(strict.tolist(), [True, False, False, False])
        self.assertEqual(broad.tolist(), [True, True, False, False])
        self.assertTrue(np.all(broad[strict]))

    def test_primary_fraction_excludes_broad_negatives(self) -> None:
        """Unlabelled LOY cells must not be counted as mCherry-class CTRL cells."""
        sample = parse_sample("A04-LOY_GFP_WT_1_001")
        points = np.array([[9.0, -5.0]] * 3 + [[-2.0, 10.0]] * 3 + [[0.0, -5.0], [0.0, 5.0]])
        row = measure_well(points, GFP_AXIS, GFP_AXIS, sample)
        self.assertEqual((row["primary_singlets"], row["reporter_positive_singlets"], row["loy_count"]), (8, 6, 3))
        self.assertAlmostEqual(row["loy_fraction"], 0.5)
        self.assertAlmostEqual(row["loy_fraction_all_events"], 3 / 8)
        self.assertAlmostEqual(row["negative_broad_fraction"], 0.25)

    def test_unstained_check_separates_ambiguous_from_kept_mcherry(self) -> None:
        points = np.array([[0.0, 5.0], [0.0, 5.0], [0.0, 9.0], [0.0, 9.0], [0.0, -5.0]])
        raw = np.array([[0, 50_000], [0, 300_000], [0, 600_000], [0, 100_000], [0, 1_000]], dtype=float)
        check = unstained_check(points, raw, GFP_AXIS, {"mcherry_threshold": 186_360.0})
        self.assertAlmostEqual(check["ambiguous_fraction"], 0.4)
        self.assertAlmostEqual(check["ambiguous_below_unstained_fraction"], 0.5)
        self.assertAlmostEqual(check["kept_mcherry_above_unstained_fraction"], 0.5)
        self.assertEqual(check["raw_mcherry_median_strict_negatives"], 1_000.0)


if __name__ == "__main__":
    unittest.main()
