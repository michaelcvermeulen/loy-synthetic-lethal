from __future__ import annotations

import csv
import json
import tempfile
import unittest
from pathlib import Path

from spectral.a8_pipeline.longitudinal_report import (
    analyze_new_acquisition,
    build_longitudinal_result,
    build_new_run_design,
    parse_new_run_sample,
    pool_new_run_measurements,
    verify_unmixing_compatibility,
    write_longitudinal_outputs,
)


class NewRunDesignTests(unittest.TestCase):
    def test_assigns_run1_wells_to_d7_with_two_technical_replicates(self) -> None:
        """Confusing acquisition-run labels with timepoints would mislabel A/B wells."""
        first = parse_new_run_sample(
            Path("A02-LOY_GFP_CTRL_CHERRY_rep1_run1_.fcs")
        )
        second = parse_new_run_sample(
            Path("B02-LOY_GFP_CTRL_CHERRY_rep1_run1_001.fcs")
        )

        self.assertEqual(
            (first.experiment_run, first.timepoint, first.orientation),
            ("run1", "D7", "LOY_GFP"),
        )
        self.assertEqual((first.biological_replicate, first.technical_replicate), (1, 1))
        self.assertEqual((second.biological_replicate, second.technical_replicate), (1, 2))

    def test_assigns_run2_wells_to_d0_with_two_technical_replicates(self) -> None:
        """Run 2 baseline wells must not be added to Run 1's D7 endpoint."""
        first = parse_new_run_sample(
            Path("D05-LOY_CHERRY_CTRL_GFP_rep1_run2_002.fcs")
        )
        second = parse_new_run_sample(
            Path("E05-LOY_CHERRY_CTRL_GFP_rep1_run2_003.fcs")
        )

        self.assertEqual(
            (first.experiment_run, first.timepoint, first.orientation),
            ("run2", "D0", "LOY_mCherry"),
        )
        self.assertEqual((first.biological_replicate, first.technical_replicate), (1, 1))
        self.assertEqual((second.biological_replicate, second.technical_replicate), (1, 2))

    def test_real_bundle_has_three_biological_replicates_and_two_technical_wells(self) -> None:
        """Missing or duplicated wells must be detected before quantitative analysis."""
        root = (
            Path(__file__).resolve().parents[3]
            / "spectral"
            / "Mike_Run 2"
            / "run2"
        )
        design = build_new_run_design(root.glob("*.fcs"))

        self.assertEqual(len(design), 24)
        groups: dict[tuple[str, str, str, int], list[int]] = {}
        for sample in design:
            key = (
                sample.experiment_run,
                sample.timepoint,
                sample.orientation,
                sample.biological_replicate,
            )
            groups.setdefault(key, []).append(sample.technical_replicate)
        self.assertEqual(len(groups), 12)
        self.assertTrue(all(sorted(replicates) == [1, 2] for replicates in groups.values()))

    def test_rejects_incomplete_technical_pairs(self) -> None:
        """A single acquired well must not be treated as a complete biological replicate."""
        with self.assertRaisesRegex(ValueError, "expected technical replicates 1 and 2"):
            build_new_run_design(
                [Path("A02-LOY_GFP_CTRL_CHERRY_rep1_run1_.fcs")]
            )

    def test_pools_technical_wells_by_event_counts(self) -> None:
        """Unequal well yields must be count-pooled rather than averaged as proportions."""
        samples = build_new_run_design(
            [
                Path("A02-LOY_GFP_CTRL_CHERRY_rep1_run1_.fcs"),
                Path("B02-LOY_GFP_CTRL_CHERRY_rep1_run1_001.fcs"),
            ]
        )
        measurements = {
            samples[0].path: {"loy_count": 10, "primary_singlets": 20, "saturated_fraction": 0.1},
            samples[1].path: {"loy_count": 30, "primary_singlets": 40, "saturated_fraction": 0.2},
        }

        pooled = pool_new_run_measurements(samples, measurements)

        self.assertEqual(len(pooled), 1)
        self.assertEqual(pooled[0]["loy_count"], 40)
        self.assertEqual(pooled[0]["primary_singlets"], 60)
        self.assertAlmostEqual(pooled[0]["loy_fraction"], 2 / 3)
        self.assertAlmostEqual(pooled[0]["technical_difference_pp"], 25.0)

    def test_builds_matched_run1_trajectory_and_separate_run2_baseline(self) -> None:
        """Run 2 D0 biological replicates must remain separate from Run 1 trajectories."""
        old = {
            "source_runs": {"d0": "old_d0", "d4": "old_d4"},
            "summary": {"common_delta": -0.4},
            "orientations": {
                "LOY_GFP": {
                    "d0_loy": 0.50,
                    "plates": [{"replicate": 1, "d4_loy": 0.40}],
                },
                "LOY_mCherry": {
                    "d0_loy": 0.45,
                    "plates": [{"replicate": 1, "d4_loy": 0.35}],
                },
            },
        }
        new = [
            {
                "experiment_run": "run1",
                "timepoint": "D7",
                "orientation": "LOY_GFP",
                "biological_replicate": 1,
                "loy_fraction": 0.30,
            },
            {
                "experiment_run": "run1",
                "timepoint": "D7",
                "orientation": "LOY_mCherry",
                "biological_replicate": 1,
                "loy_fraction": 0.25,
            },
            {
                "experiment_run": "run2",
                "timepoint": "D0",
                "orientation": "LOY_GFP",
                "biological_replicate": 1,
                "loy_fraction": 0.55,
            },
            {
                "experiment_run": "run2",
                "timepoint": "D0",
                "orientation": "LOY_mCherry",
                "biological_replicate": 1,
                "loy_fraction": 0.48,
            },
        ]

        result = build_longitudinal_result(old, new)

        self.assertEqual(result["run1"][0]["timepoints"], {"D0": 0.5, "D4": 0.4, "D7": 0.3})
        self.assertEqual([row["loy_fraction"] for row in result["run2_d0"]], [0.55, 0.48])
        self.assertIn("d0_to_d7", result["effects"])
        self.assertIn("d4_to_d7", result["effects"])

    def test_rejects_changed_spectral_unmixing(self) -> None:
        """A changed unmixing matrix must prevent cross-acquisition classification."""
        reference = {"bdspectral matrix name": "matrix", "spill": "same-values"}
        verify_unmixing_compatibility(reference, dict(reference))
        with self.assertRaisesRegex(ValueError, "spectral unmixing differs"):
            verify_unmixing_compatibility(
                reference,
                {"bdspectral matrix name": "matrix", "spill": "changed-values"},
            )

    def test_analyzes_all_real_new_acquisition_pairs(self) -> None:
        """The on-disk acquisition must yield one pooled row per biological replicate."""
        repo = Path(__file__).resolve().parents[3]

        result = analyze_new_acquisition(repo)

        self.assertEqual(len(result["pooled_rows"]), 12)
        self.assertEqual(result["sample_count"], 24)
        self.assertEqual(len(result["source_manifest"]), 24)
        self.assertTrue(
            all(len(row["sha256"]) == 64 and row["size_bytes"] > 0 for row in result["source_manifest"])
        )
        self.assertTrue(result["unmixing_compatible"])
        self.assertTrue(
            all(
                row["primary_singlets"] > 0 and 0 < row["loy_fraction"] < 1
                for row in result["pooled_rows"]
            )
        )

    def test_writes_machine_readable_longitudinal_outputs(self) -> None:
        """Generated artifacts must preserve trajectories, pooled rows, and source hashes."""
        result = {
            "title": "test",
            "source_runs": {},
            "run1": [
                {
                    "orientation": "LOY_GFP",
                    "biological_replicate": 1,
                    "timepoints": {"D0": 0.5, "D4": 0.4, "D7": 0.3},
                    "delta_log_odds_d0_d4": -0.4,
                    "delta_log_odds_d0_d7": -0.8,
                    "delta_log_odds_d4_d7": -0.4,
                    "d7_measurement": {},
                }
            ],
            "run2_d0": [
                {
                    "experiment_run": "run2",
                    "timepoint": "D0",
                    "orientation": "LOY_GFP",
                    "biological_replicate": 1,
                    "loy_count": 50,
                    "primary_singlets": 100,
                    "loy_fraction": 0.5,
                    "technical_difference_pp": 1.0,
                    "saturated_fraction": 0.1,
                    "technical_wells": [],
                }
            ],
            "effects": {
                "d0_to_d4": {"common_delta": -0.4, "interval_odds_ratio": 0.67},
                "d0_to_d7": {"common_delta": -0.8, "interval_odds_ratio": 0.45},
                "d4_to_d7": {"common_delta": -0.4, "interval_odds_ratio": 0.67},
            },
            "new_acquisition": {
                "source_manifest": [{"path": "sample.fcs", "size_bytes": 10, "sha256": "a" * 64}],
                "pooled_rows": [],
                "sample_count": 1,
                "unmixing_compatible": True,
                "unmixing_matrix": "matrix",
            },
        }
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "report"
            write_longitudinal_outputs(result, output)
            payload = json.loads((output / "results.json").read_text(encoding="utf-8"))
            with (output / "run1_trajectories.csv").open(newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))

            self.assertEqual(payload["title"], "test")
            self.assertEqual(rows[0]["D7"], "0.3")
            self.assertTrue((output / "source_manifest.csv").is_file())
            figure = output / "figures" / "longitudinal_loy_fraction.svg"
            self.assertIn("<polyline", figure.read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
