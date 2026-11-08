from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from spectral.a8_pipeline.report import write_report
from spectral.a8_pipeline.sensitivity import (
    Scenario,
    build_sensitivity_scenarios,
    default_sensitivity_config,
    summarise_families,
)
from spectral.a8_pipeline.svg import (
    write_reporter_density_svg,
    write_saturation_composition_svg,
    write_sensitivity_delta_svg,
    write_time_qc_svg,
)

FAMILIES = {"time", "segment", "cell", "singlet", "saturation", "threshold", "cofactor"}


class SensitivityScenarioTests(unittest.TestCase):
    def test_sensitivity_contains_predefined_families(self) -> None:
        scenarios = build_sensitivity_scenarios(default_sensitivity_config())
        self.assertEqual({row.family for row in scenarios}, FAMILIES)

    def test_exactly_one_baseline_scenario(self) -> None:
        scenarios = build_sensitivity_scenarios(default_sensitivity_config())
        baselines = [row for row in scenarios if row.is_baseline]
        self.assertEqual(len(baselines), 1)
        base = baselines[0]
        self.assertTrue(base.time_qc)
        self.assertEqual(base.saturation_rule, "spectral")
        self.assertEqual(base.polygon_scale, 1.0)
        self.assertEqual(base.singlet_mad, 3.5)
        self.assertEqual(base.cofactor, 10000.0)
        self.assertEqual((base.gfp_shift, base.mcherry_shift), (0.0, 0.0))

    def test_threshold_family_covers_all_nine_shift_combinations(self) -> None:
        scenarios = build_sensitivity_scenarios(default_sensitivity_config())
        shifts = {
            (row.gfp_shift, row.mcherry_shift)
            for row in scenarios
            if row.family == "threshold"
        }
        self.assertEqual(len(shifts), 9)

    def test_summarise_reports_maximum_absolute_change_per_family(self) -> None:
        frame = pd.DataFrame(
            {
                "family": ["time", "time", "singlet", "singlet"],
                "sample": ["F11"] * 4,
                "is_baseline": [True, False, False, False],
                "loy_fraction": [0.50, 0.54, 0.47, 0.52],
                "log2_ratio": [0.0, 0.20, -0.15, 0.10],
            }
        )
        summary = summarise_families(frame, baseline_fraction=0.50, baseline_log2=0.0)
        time_row = summary.set_index("family").loc["time"]
        self.assertAlmostEqual(float(time_row["max_abs_delta_loy_fraction"]), 0.04, places=6)
        singlet_row = summary.set_index("family").loc["singlet"]
        self.assertAlmostEqual(float(singlet_row["max_abs_delta_log2_ratio"]), 0.15, places=6)


class ScenarioDefaultsTests(unittest.TestCase):
    def test_scenario_defaults_match_primary_analysis(self) -> None:
        scenario = Scenario(family="time", label="baseline")
        self.assertTrue(scenario.time_qc)
        self.assertEqual(scenario.segment, "all")
        self.assertEqual(scenario.saturation_rule, "spectral")


class SvgTests(unittest.TestCase):
    def test_reporter_density_svg_marks_thresholds(self) -> None:
        rng = np.random.default_rng(3)
        values = np.r_[rng.normal(-2, 0.4, 500), rng.normal(4, 0.5, 500)]
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary, "reporter.svg")
            write_reporter_density_svg(values, threshold=1.0, output=output, title="eGFP-A")
            text = output.read_text(encoding="utf-8")
            self.assertIn("<svg", text)
            self.assertIn("</svg>", text)
            self.assertIn("stroke-dasharray", text)

    def test_time_qc_svg_marks_flagged_bins(self) -> None:
        bins = pd.DataFrame(
            {
                "segment": ["Run 1"] * 4,
                "bin": [1, 2, 3, 4],
                "time_start": [0.0, 1.0, 2.0, 3.0],
                "event_rate": [100.0, 102.0, 40.0, 99.0],
                "bad": [False, False, True, False],
            }
        )
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary, "time.svg")
            write_time_qc_svg(bins, output=output, title="F11 acquisition")
            self.assertIn("<svg", output.read_text(encoding="utf-8"))

    def test_saturation_and_delta_svgs_render(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            saturation = pd.DataFrame(
                {
                    "sample": ["F11", "G11"],
                    "any_saturated": [56.12, 47.03],
                    "spectral_saturated": [0.40, 0.93],
                    "imgb3_saturated": [56.09, 46.80],
                }
            )
            first = Path(temporary, "sat.svg")
            write_saturation_composition_svg(saturation, output=first, title="Saturation")
            self.assertIn("<svg", first.read_text(encoding="utf-8"))

            summary = pd.DataFrame(
                {
                    "family": ["time", "singlet"],
                    "max_abs_delta_loy_fraction": [0.004, 0.011],
                    "max_abs_delta_log2_ratio": [0.02, 0.05],
                }
            )
            second = Path(temporary, "delta.svg")
            write_sensitivity_delta_svg(summary, output=second, title="Sensitivity")
            self.assertIn("<svg", second.read_text(encoding="utf-8"))


class ReportTests(unittest.TestCase):
    def test_report_labels_counting_precision_and_viability_limit(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            data = {
                "title": "Synthetic A8 report",
                "interval_note": "Intervals show counting precision.",
                "limitations": ["There was no viability dye."],
                "tables": {},
                "figures": [],
            }
            paths = write_report(data, Path(temporary))
            text = paths.markdown.read_text(encoding="utf-8")
            self.assertIn("counting precision", text)
            self.assertIn("no viability dye", text)
            self.assertTrue(paths.html.exists())
            self.assertIn("<html", paths.html.read_text(encoding="utf-8").lower())

    def test_report_renders_tables_and_figure_links(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            data = {
                "title": "Synthetic",
                "interval_note": "counting precision",
                "limitations": ["no viability dye"],
                "tables": {
                    "Competition": pd.DataFrame({"sample": ["F11"], "loy_fraction": [0.318]})
                },
                "figures": [{"caption": "Scatter", "path": "figures/scatter.svg"}],
            }
            paths = write_report(data, Path(temporary))
            markdown = paths.markdown.read_text(encoding="utf-8")
            self.assertIn("Competition", markdown)
            self.assertIn("0.318", markdown)
            self.assertIn("figures/scatter.svg", markdown)


if __name__ == "__main__":
    unittest.main()
