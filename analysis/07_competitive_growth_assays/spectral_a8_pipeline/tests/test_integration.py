from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from spectral.a8_pipeline.manifest import build_manifest, verify_manifest
from spectral.a8_pipeline.run_analysis import (
    OutputLocationError,
    resolve_output,
    run_pipeline,
)

INPUT_ROOT = Path(r"A:\flow_jo_AUG28")
CONFIG = Path("spectral/a8_pipeline/default_config.json")


class OutputValidationTests(unittest.TestCase):
    """These must not touch the data, so they stay outside the shared dry run."""

    def test_rejects_output_equal_to_input_root(self) -> None:
        with self.assertRaisesRegex(OutputLocationError, "must not be"):
            resolve_output(INPUT_ROOT, INPUT_ROOT)

    def test_rejects_output_above_input_root(self) -> None:
        with self.assertRaisesRegex(OutputLocationError, "must not be"):
            resolve_output(INPUT_ROOT, INPUT_ROOT.parent)

    def test_accepts_output_beneath_input_root(self) -> None:
        output, staging = resolve_output(INPUT_ROOT, INPUT_ROOT / "programmatic_analysis")
        self.assertEqual(output.name, "programmatic_analysis")
        self.assertEqual(staging.name, ".programmatic_analysis_staging")

    def test_refuses_to_overwrite_an_existing_output(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary, "programmatic_analysis")
            output.mkdir()
            (output / "keep.txt").write_text("existing", encoding="utf-8")
            with self.assertRaisesRegex(OutputLocationError, "already exists"):
                run_pipeline(INPUT_ROOT, output, CONFIG, dry_run=False)
            self.assertEqual((output / "keep.txt").read_text(encoding="utf-8"), "existing")


class IntegrationTests(unittest.TestCase):
    """One dry run shared by every assertion.

    A full run hashes 9.55 GB across 28,784 source files twice, so running the
    pipeline once per test method would dominate the suite runtime.
    """

    @classmethod
    def setUpClass(cls) -> None:
        cls._manifest_before = build_manifest(INPUT_ROOT)
        cls._temporary = tempfile.TemporaryDirectory()
        cls.output = Path(cls._temporary.name, "programmatic_analysis")
        cls.result = run_pipeline(INPUT_ROOT, cls.output, CONFIG, dry_run=True)

    @classmethod
    def tearDownClass(cls) -> None:
        cls._temporary.cleanup()

    def test_dry_run_uses_expected_files_without_mutation(self) -> None:
        self.assertEqual(self.result.raw_counts["F11"], 15000)
        self.assertEqual(self.result.raw_counts["G11"], 20000)
        self.assertEqual(self.result.segment_counts["G11"], [10000, 10000])
        self.assertEqual(self.result.segment_counts["F11"], [15000])
        verify_manifest(self._manifest_before)

    def test_dry_run_writes_nothing(self) -> None:
        self.assertFalse(self.output.exists())
        self.assertFalse(self.result.staging.exists())

    def test_dry_run_reports_expected_image_coverage(self) -> None:
        self.assertEqual(self.result.image_coverage["F11"], (8776, 15000))
        self.assertEqual(self.result.image_coverage["G11"], (20000, 20000))

    def test_dry_run_confirms_saturation_decoding_and_reporters(self) -> None:
        self.assertTrue(all(self.result.saturation_agreement.values()))
        self.assertEqual(self.result.matrix_names["F11"], "mcherry_GFP_matrix")
        self.assertEqual(self.result.matrix_names["G11"], "mcherry_GFP_matrix")

    def test_gate_counts_are_nested_and_reconcile(self) -> None:
        counts = self.result.gate_counts.set_index("sample")
        for sample in ("F11", "G11"):
            row = counts.loc[sample]
            self.assertEqual(row["raw_events"], self.result.raw_counts[sample])
            self.assertGreaterEqual(row["raw_events"], row["after_time_qc"])
            self.assertGreaterEqual(row["after_time_qc"], row["after_spectral_saturation"])
            self.assertGreaterEqual(row["after_spectral_saturation"], row["after_cell_polygon"])
            self.assertGreaterEqual(row["after_cell_polygon"], row["primary_singlets"])

    def test_competition_covers_both_orientations_and_segments(self) -> None:
        competition = self.result.competition
        overall = competition[competition["scope"] == "all"].set_index("sample")
        self.assertEqual(overall.loc["F11", "loy_reporter"], "GFP")
        self.assertEqual(overall.loc["G11", "loy_reporter"], "mCherry")
        g11_scopes = set(competition[competition["sample"] == "G11"]["scope"])
        self.assertEqual(g11_scopes, {"all", "Run 1", "Run 2"})

    def test_sensitivity_covers_every_family(self) -> None:
        families = set(self.result.family_summary["family"])
        self.assertEqual(
            families,
            {"time", "segment", "cell", "singlet", "saturation", "threshold", "cofactor"},
        )


if __name__ == "__main__":
    unittest.main()
