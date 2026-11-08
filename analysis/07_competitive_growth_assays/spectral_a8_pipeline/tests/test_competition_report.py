from __future__ import annotations

import json
import math
import tempfile
import unittest
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np

from spectral.a8_pipeline import competition_report as report_module
from spectral.a8_pipeline.competition_report import (
    build_report_markdown,
    delta_log_odds,
    expected_gfp_fraction,
    fit_lda,
    load_cef_gates,
    pool_technical_counts,
    summarise_effects,
    write_competition_svg,
    write_outputs,
    write_ratio_ladder_svg,
)


class GateParsingTests(unittest.TestCase):
    def test_load_cef_gates_preserves_parameters_and_vertices(self) -> None:
        """A swapped parameter or vertex would apply the gate on the wrong axes."""
        payload = {
            "experiment": {
                "panels": [
                    {
                        "analysis": {
                            "gates": [
                                {
                                    "name": "Scatter",
                                    "gateKind": "Polygon",
                                    "gateId": "1",
                                    "parentPopulationId": "0-1",
                                    "parameters": [{"name": "X-A"}, {"name": "Y-A"}],
                                    "vertices": [
                                        {"x": 1.0, "y": 2.0},
                                        {"x": 3.0, "y": 4.0},
                                        {"x": 5.0, "y": 6.0},
                                    ],
                                }
                            ]
                        }
                    }
                ]
            }
        }
        with tempfile.TemporaryDirectory() as temporary:
            cef = Path(temporary, "synthetic.cef")
            with zipfile.ZipFile(cef, "w") as archive:
                archive.writestr("experiment.json", json.dumps(payload))
            gates = load_cef_gates(cef)

        gate = gates["Scatter"]
        self.assertEqual(gate.parameters, ("X-A", "Y-A"))
        np.testing.assert_array_equal(
            gate.vertices,
            np.array([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]]),
        )


class ReporterAndStatisticsTests(unittest.TestCase):
    def test_background_gate_uses_unstained_channel_quantiles(self) -> None:
        """Per-channel thresholds must be derived from the matching unstained cloud."""
        unstained = np.array([[0.0, 0.0], [1.0, 2.0], [2.0, 1.0], [3.0, 4.0]])
        thresholds = report_module.background_gate_thresholds(unstained, percentile=75)
        np.testing.assert_allclose(thresholds, [2.25, 2.5])
        sample = np.array([[2.0, 2.0], [2.3, 1.0], [1.0, 2.6]])
        calls = report_module.apply_background_gate(sample, thresholds)
        np.testing.assert_array_equal(calls, [0, 1, 2])

    def test_unstained_control_is_background_not_a_reporter_class(self) -> None:
        """Unstained cells must not be forced into GFP-like or mCherry-like calls."""
        classifier = report_module.LDAModel(
            weights=np.array([1.0, -1.0]), intercept=0.0
        )
        features = np.array([[2.0, 0.0], [0.0, 2.0], [1.0, 1.0]])
        payload = report_module.build_all_sample_plot_payload(
            classifier,
            {"A01-unstained_CTRL": {"category": "controls", "features": features}},
            maximum_points_per_sample=20,
        )
        panel = payload["groups"]["controls"][0]
        self.assertEqual(panel["gfp_count"], 0)
        self.assertEqual(panel["gfp_fraction"], 0.0)
        self.assertEqual(panel["mcherry_fraction"], 0.0)
        self.assertTrue(panel["unstained"])
        self.assertEqual(panel["series"][0]["reporter_class"], -1)
        self.assertEqual(len(panel["series"][0]["points"]), 3)

    def test_all_sample_payload_uses_full_counts_and_proportionate_display_points(self) -> None:
        """Plot subsampling must not replace full-event reporter percentages."""
        self.assertTrue(hasattr(report_module, "build_all_sample_plot_payload"))
        classifier = report_module.LDAModel(
            weights=np.array([1.0, -1.0]), intercept=0.0
        )
        features = np.array(
            [[2.0, 0.0], [3.0, 0.0], [4.0, 0.0], [0.0, 2.0]], dtype=float
        )
        samples = {
            "A01-control": {"category": "controls", "features": features},
            "B01-ratio": {"category": "ratio_ladder", "features": features},
            "D01-endpoint": {
                "category": "competition_endpoints",
                "features": features,
            },
        }
        payload = report_module.build_all_sample_plot_payload(
            classifier, samples, maximum_points_per_sample=20
        )

        self.assertEqual(
            {name: len(panels) for name, panels in payload["groups"].items()},
            {"controls": 1, "ratio_ladder": 1, "competition_endpoints": 1},
        )
        for panels in payload["groups"].values():
            panel = panels[0]
            self.assertEqual(panel["primary_singlets"], 4)
            self.assertEqual(panel["gfp_count"], 3)
            self.assertAlmostEqual(panel["gfp_fraction"], 0.75)
            counts = {
                series["reporter_class"]: len(series["points"])
                for series in panel["series"]
            }
            self.assertEqual(counts, {1: 3, 0: 1})
        self.assertLess(payload["x_limits"][0], 0.0)
        self.assertGreater(payload["x_limits"][1], 4.0)

    def test_plot_payload_uses_controls_and_preserves_mixture_class_counts(self) -> None:
        """Balancing plotted classes would visually misrepresent mixture composition."""
        self.assertTrue(hasattr(report_module, "build_classifier_plot_payload"))
        classifier = report_module.LDAModel(
            weights=np.array([1.0, -1.0]), intercept=0.0
        )
        controls = {
            "A03-CTRL_GFP": np.array([[2.0, 0.0], [3.0, 0.0]]),
            "A05-LOY_GFP": np.array([[2.0, 0.2], [3.0, 0.2]]),
            "A07-CTRL_cherry": np.array([[0.0, 2.0], [0.0, 3.0]]),
            "A09-LOY_cherry": np.array([[0.2, 2.0], [0.2, 3.0]]),
        }
        mixed = np.array([[2.0, 0.0], [3.0, 0.0], [4.0, 0.0], [0.0, 2.0]])
        payload = report_module.build_classifier_plot_payload(
            classifier,
            controls,
            {"LOY_GFP": mixed, "LOY_mCherry": mixed},
            {"LOY_GFP": mixed, "LOY_mCherry": mixed},
            maximum_points_per_panel=20,
        )

        self.assertEqual(len(payload["panels"]), 5)
        self.assertEqual(
            [series["label"] for series in payload["panels"][0]["series"]],
            ["CTRL-GFP", "LOY-GFP", "CTRL-mCherry", "LOY-mCherry"],
        )
        for panel in payload["panels"][1:]:
            counts = {
                series["reporter_class"]: len(series["points"])
                for series in panel["series"]
            }
            self.assertEqual(counts, {1: 3, 0: 1})
        self.assertEqual(payload["weights"], [1.0, -1.0])
        self.assertEqual(payload["intercept"], 0.0)

    def test_lda_separates_control_derived_reporter_clusters(self) -> None:
        """Reversing the discriminant sign would swap GFP and mCherry calls."""
        features = np.array([[-2.0, 0.0], [-1.5, 0.2], [1.5, -0.2], [2.0, 0.0]])
        labels = np.array([0, 0, 1, 1])
        model = fit_lda(features, labels)
        self.assertEqual(model.predict(np.array([[-3.0, 0.0], [3.0, 0.0]])).tolist(), [0, 1])

    def test_expected_gfp_fraction_respects_dye_orientation(self) -> None:
        """Treating the LOY:CTRL ratio as GFP:mCherry would reverse one ladder."""
        self.assertAlmostEqual(expected_gfp_fraction(0.25, "LOY_GFP"), 0.25)
        self.assertAlmostEqual(expected_gfp_fraction(0.25, "CTRL_GFP"), 0.75)

    def test_pooling_sums_counts_instead_of_averaging_percentages(self) -> None:
        """Unequal event yields require count pooling, not an unweighted mean."""
        positive, total, fraction = pool_technical_counts([(10, 100), (30, 50)])
        self.assertEqual((positive, total), (40, 150))
        self.assertAlmostEqual(fraction, 40 / 150)

    def test_delta_log_odds_uses_measured_d0(self) -> None:
        """Subtracting raw percentages would not estimate relative competition."""
        self.assertAlmostEqual(delta_log_odds(0.5, 0.25), -math.log(3.0))

    def test_effect_summary_balances_the_two_dye_orientations(self) -> None:
        """Pooling without an orientation term could let one dye arm dominate."""
        summary = summarise_effects(
            {"LOY_GFP": [-0.5, -0.7, -0.6], "LOY_mCherry": [-0.8, -0.6, -0.7]}
        )
        self.assertAlmostEqual(summary["common_delta"], -0.65)
        self.assertAlmostEqual(summary["interval_odds_ratio"], math.exp(-0.65))
        self.assertLess(summary["ci_low"], -0.65)
        self.assertGreater(summary["ci_high"], -0.65)


class ArtifactTests(unittest.TestCase):
    def test_all_sample_sheet_labels_each_sample_with_shared_axes_and_boundary(self) -> None:
        """Per-sample sheets must remain comparable and identify every plotted sample."""
        self.assertTrue(hasattr(report_module, "write_all_sample_sheet_svg"))
        panels = []
        for name in ("A01-control-one", "A02-control-two"):
            panels.append(
                {
                    "sample": name,
                    "primary_singlets": 4,
                    "gfp_count": 3,
                    "gfp_fraction": 0.75,
                    "mcherry_fraction": 0.25,
                    "series": [
                        {
                            "label": "GFP-like",
                            "reporter_class": 1,
                            "points": [[2.0, 0.0], [3.0, 0.0], [4.0, 0.0]],
                        },
                        {
                            "label": "mCherry-like",
                            "reporter_class": 0,
                            "points": [[0.0, 2.0]],
                        },
                    ],
                }
            )
        payload = {
            "weights": [1.0, -1.0],
            "intercept": 0.0,
            "x_limits": [-0.5, 4.5],
            "y_limits": [-0.5, 2.5],
            "groups": {
                "controls": panels,
                "ratio_ladder": [],
                "competition_endpoints": [],
            },
        }
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary, "controls.svg")
            report_module.write_all_sample_sheet_svg(payload, "controls", path)
            root = ET.parse(path).getroot()

        facets = [item for item in root.iter() if item.attrib.get("data-sample")]
        boundaries = [item for item in root.iter() if item.attrib.get("data-boundary")]
        labels = " ".join(item.text or "" for item in root.iter())
        self.assertEqual(len(facets), 2)
        self.assertEqual(len(boundaries), 2)
        self.assertEqual({item.attrib["data-x-min"] for item in facets}, {"-0.5"})
        self.assertEqual({item.attrib["data-x-max"] for item in facets}, {"4.5"})
        self.assertIn("A01-control-one", labels)
        self.assertIn("A02-control-two", labels)
        self.assertIn("n=4", labels)
        self.assertIn("GFP-like 75.0%", labels)
        font_sizes = [
            float(item.attrib["font-size"])
            for item in root.iter()
            if "font-size" in item.attrib
        ]
        self.assertGreaterEqual(min(font_sizes), 11.0)

    def test_workflow_svg_has_four_connected_analysis_stages(self) -> None:
        """Dropping a stage would make the classification procedure misleading."""
        self.assertTrue(hasattr(report_module, "write_classifier_workflow_svg"))
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary, "workflow.svg")
            report_module.write_classifier_workflow_svg(path)
            root = ET.parse(path).getroot()

        nodes = [item for item in root.iter() if item.attrib.get("data-stage")]
        arrows = [item for item in root.iter() if item.attrib.get("data-connector")]
        labels = " ".join(item.text or "" for item in root.iter())
        self.assertEqual(len(nodes), 4)
        self.assertEqual(len(arrows), 3)
        self.assertIn("Scatter", labels)
        self.assertIn("arcsinh", labels)
        self.assertIn("single-colour controls", labels)
        self.assertIn("score > 0: GFP-like", labels)
        self.assertIn("score ≤ 0: mCherry-like", labels)
        self.assertIn("Unstained controls: 0%/0% background", labels)

    def test_reporter_space_svg_shows_shared_boundary_in_all_five_panels(self) -> None:
        """A missing boundary or facet could hide how reporter calls transfer across days."""
        self.assertTrue(hasattr(report_module, "write_reporter_space_svg"))
        panels = []
        for index, title in enumerate(
            [
                "Training controls",
                "D0 | LOY-GFP / CTRL-mCherry",
                "D4 | LOY-GFP / CTRL-mCherry",
                "D0 | LOY-mCherry / CTRL-GFP",
                "D4 | LOY-mCherry / CTRL-GFP",
            ]
        ):
            panels.append(
                {
                    "title": title,
                    "series": [
                        {
                            "label": "GFP-like",
                            "reporter_class": 1,
                            "points": [[1.0 + index, 0.0], [1.2 + index, 0.1]],
                        },
                        {
                            "label": "mCherry-like",
                            "reporter_class": 0,
                            "points": [[0.0, 1.0 + index], [0.1, 1.2 + index]],
                        },
                    ],
                }
            )
        payload = {"weights": [1.0, -1.0], "intercept": 0.0, "panels": panels}
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary, "reporters.svg")
            report_module.write_reporter_space_svg(payload, path)
            root = ET.parse(path).getroot()

        facets = [item for item in root.iter() if item.attrib.get("data-panel")]
        boundaries = [item for item in root.iter() if item.attrib.get("data-boundary")]
        marks = [item for item in root.iter() if item.attrib.get("data-reporter-class")]
        ticks = [item for item in root.iter() if item.attrib.get("data-axis-tick")]
        labels = " ".join(item.text or "" for item in root.iter())
        self.assertEqual(len(facets), 5)
        self.assertEqual(len(boundaries), 5)
        self.assertEqual(len(marks), 20)
        self.assertEqual(len(ticks), 40)
        self.assertIn("arcsinh(eGFP-A / 1000)", labels)
        self.assertIn("arcsinh(mCherry-A / 1000)", labels)
        self.assertIn("CTRL control", labels)
        self.assertIn("LOY control", labels)
        legend = next(item for item in root.iter() if item.attrib.get("aria-label") == "legend")
        legend_bottom = max(float(item.attrib["y"]) for item in legend.iter() if "y" in item.attrib)
        training_title = next(item for item in root.iter() if item.text == "Training controls")
        self.assertGreaterEqual(float(training_title.attrib["y"]) - legend_bottom, 30.0)
        x_axis_labels = [
            item for item in root.iter() if item.text == "arcsinh(eGFP-A / 1000)"
        ]
        next_row_titles = [
            item
            for item in root.iter()
            if item.text and (item.text.startswith("D0 |") or item.text.startswith("D4 |"))
        ]
        self.assertGreaterEqual(
            min(float(item.attrib["y"]) for item in next_row_titles)
            - min(float(item.attrib["y"]) for item in x_axis_labels),
            30.0,
        )

    def test_all_sample_sheet_marks_unstained_background_as_zero_zero(self) -> None:
        """The controls sheet must show unstained as 0%/0%, not a forced class split."""
        panel = {
            "sample": "A01-unstained_CTRL",
            "primary_singlets": 3,
            "gfp_count": 0,
            "gfp_fraction": 0.0,
            "mcherry_fraction": 0.0,
            "unstained": True,
            "series": [
                {
                    "label": "Unstained background",
                    "reporter_class": -1,
                    "points": [[0.0, 0.0], [0.1, 0.1], [0.2, 0.2]],
                }
            ],
        }
        payload = {
            "weights": [1.0, -1.0],
            "intercept": 0.0,
            "x_limits": [-0.5, 1.0],
            "y_limits": [-0.5, 1.0],
            "groups": {"controls": [panel], "ratio_ladder": [], "competition_endpoints": []},
        }
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary, "controls.svg")
            report_module.write_all_sample_sheet_svg(payload, "controls", path)
            root = ET.parse(path).getroot()
        labels = " ".join(item.text or "" for item in root.iter())
        background = [
            item for item in root.iter() if item.attrib.get("data-reporter-class") == "-1"
        ]
        self.assertEqual(len(background), 3)
        self.assertIn("GFP+ 0.0%", labels)
        self.assertIn("mCherry+ 0.0%", labels)
        self.assertIn("unstained background", labels)

    def test_write_outputs_accepts_preexisting_empty_directory(self) -> None:
        """Dropbox may briefly lock an emptied output directory during --force."""
        result = {
            "summary": {
                "common_delta": -0.697,
                "ci_low": -0.863,
                "ci_high": -0.530,
                "interval_odds_ratio": 0.498,
            },
            "orientations": {},
            "ratio_ladder": [],
            "classifier": {"name": "LDA", "maximum_control_error": 0.01},
            "classifier_plot": {
                "weights": [1.0, -1.0],
                "intercept": 0.0,
                "panels": [],
            },
            "all_sample_plots": {
                "weights": [1.0, -1.0],
                "intercept": 0.0,
                "x_limits": [0.0, 1.0],
                "y_limits": [0.0, 1.0],
                "groups": {
                    "controls": [],
                    "ratio_ladder": [],
                    "competition_endpoints": [],
                },
            },
            "technical_max_difference_pp": 0.8,
        }
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary, "report")
            output.mkdir()
            (output / "figures").mkdir()
            write_outputs(result, output)
            self.assertTrue((output / "report.md").is_file())
            self.assertTrue((output / "results.json").is_file())
            self.assertTrue((output / "figures" / "classification_workflow.svg").is_file())
            self.assertTrue((output / "figures" / "reporter_classification.svg").is_file())
            self.assertTrue((output / "figures" / "all_samples" / "controls.svg").is_file())
            self.assertTrue((output / "figures" / "all_samples" / "ratio_ladder.svg").is_file())
            self.assertTrue(
                (output / "figures" / "all_samples" / "competition_endpoints.svg").is_file()
            )

    def test_clear_generated_output_reuses_directory_tree(self) -> None:
        """Force regeneration must not rmdir folders that Dropbox can transiently lock."""
        self.assertTrue(hasattr(report_module, "clear_generated_output"))
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary, "report")
            figures = output / "figures"
            figures.mkdir(parents=True)
            (output / "report.md").write_text("old", encoding="utf-8")
            (figures / "old.svg").write_text("old", encoding="utf-8")
            report_module.clear_generated_output(output)
            self.assertTrue(output.is_dir())
            self.assertTrue(figures.is_dir())
            self.assertEqual(list(output.rglob("*")), [figures])

    def test_svg_writers_include_axes_and_expected_reference(self) -> None:
        """A blank or unlabeled SVG would not communicate the QC result."""
        ladder = [
            {"orientation": "LOY_GFP", "nominal_loy": 0.25, "observed_loy": 0.26},
            {"orientation": "CTRL_GFP", "nominal_loy": 0.25, "observed_loy": 0.18},
        ]
        competition = {
            "orientations": {
                "LOY_GFP": {"d0_loy": 0.52, "plates": [{"replicate": 1, "d4_loy": 0.37}]},
                "LOY_mCherry": {"d0_loy": 0.42, "plates": [{"replicate": 1, "d4_loy": 0.22}]},
            }
        }
        with tempfile.TemporaryDirectory() as temporary:
            ratio_path = Path(temporary, "ratio.svg")
            competition_path = Path(temporary, "competition.svg")
            write_ratio_ladder_svg(ladder, ratio_path)
            write_competition_svg(competition, competition_path)
            ratio = ratio_path.read_text(encoding="utf-8")
            comp = competition_path.read_text(encoding="utf-8")
        self.assertIn("Nominal LOY", ratio)
        self.assertIn("Observed LOY", ratio)
        self.assertIn("expected", ratio)
        self.assertIn("Day 0", comp)
        self.assertIn("Day 4", comp)

    def test_report_contains_result_method_qc_and_limitation(self) -> None:
        """The report is incomplete if it omits interpretation boundaries."""
        result = {
            "summary": {
                "common_delta": -0.697,
                "ci_low": -0.863,
                "ci_high": -0.530,
                "interval_odds_ratio": 0.498,
            },
            "orientations": {},
            "ratio_ladder": [],
            "classifier": {"name": "LDA", "maximum_control_error": 0.01},
            "technical_max_difference_pp": 0.8,
        }
        text = build_report_markdown(result)
        self.assertIn("Methods", text)
        self.assertIn("Quality control", text)
        self.assertIn("shared D0", text)
        self.assertIn("0.498", text)
        self.assertIn("figures/classification_workflow.svg", text)
        self.assertIn("figures/reporter_classification.svg", text)
        self.assertIn("figures/all_samples/controls.svg", text)
        self.assertIn("figures/all_samples/ratio_ladder.svg", text)
        self.assertIn("figures/all_samples/competition_endpoints.svg", text)
        self.assertIn("Unstained controls", text)


class RepositoryIntegrationTests(unittest.TestCase):
    def test_real_analysis_carries_five_classifier_plot_panels(self) -> None:
        """Report regeneration must derive the diagnostic plots from the real FCS inputs."""
        repo = Path(__file__).resolve().parents[3]
        result = report_module.analyze_repository(repo)
        self.assertIn("classifier_plot", result)
        panels = result["classifier_plot"]["panels"]
        self.assertEqual(len(panels), 5)
        self.assertEqual(panels[0]["title"], "Training controls")
        self.assertTrue(all(sum(len(series["points"]) for series in panel["series"]) > 0 for panel in panels))
        self.assertTrue(
            all(
                sum(len(series["points"]) for series in panel["series"]) <= 800
                for panel in panels
            )
        )

    def test_real_analysis_groups_all_38_d4_samples_for_plotting(self) -> None:
        """No September acquisition may be omitted from the all-sample sheets."""
        repo = Path(__file__).resolve().parents[3]
        result = report_module.analyze_repository(repo)
        self.assertIn("all_sample_plots", result)
        groups = result["all_sample_plots"]["groups"]
        self.assertEqual(
            {name: len(panels) for name, panels in groups.items()},
            {"controls": 6, "ratio_ladder": 20, "competition_endpoints": 12},
        )
        panels = [panel for group in groups.values() for panel in group]
        self.assertEqual(len({panel["sample"] for panel in panels}), 38)
        self.assertTrue(
            all(sum(len(series["points"]) for series in panel["series"]) <= 300 for panel in panels)
        )
        unstained = [panel for panel in groups["controls"] if "unstained" in panel["sample"]]
        self.assertEqual(len(unstained), 2)
        self.assertTrue(all(panel["gfp_fraction"] == 0.0 for panel in unstained))
        self.assertTrue(all(panel["mcherry_fraction"] == 0.0 for panel in unstained))
        self.assertTrue(
            all(panel["series"][0]["reporter_class"] == -1 for panel in unstained)
        )


if __name__ == "__main__":
    unittest.main()
