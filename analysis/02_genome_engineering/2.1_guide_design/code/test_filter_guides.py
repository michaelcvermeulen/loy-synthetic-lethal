"""Focused checks for consequential filtering and missing-data behaviour."""
import csv
import json
from pathlib import Path
import tempfile
import unittest
from contextlib import redirect_stdout, redirect_stderr
from io import StringIO

import filter_guides as pipeline


DEFAULTS = {
    "minimum_exact_y_matches": 10,
    "maximum_exact_non_y_matches": 0,
    "maximum_non_y_1_mismatch_sites": 0,
    "maximum_non_y_2_mismatch_sites": 0,
}


def candidate(guide, y=10, non_y=0):
    row = {field: 0 for field in pipeline.COUNT_COLS}
    row.update(guide=guide, chrY=y, chr1=non_y, sumY=y, sumOther=non_y)
    return row


def annotation(mm1=0, mm2=0, cfd=0):
    return dict(non_y_mm1=mm1, non_y_mm2=mm2, non_y_mm3=4, non_y_mm4=6,
                legacy_non_y_cfd_mean=.1, legacy_non_y_cfd_max=.5,
                legacy_non_y_cfd_sum=cfd)


class FilterTests(unittest.TestCase):
    def test_exact_match_boundaries_and_missing_is_not_zero(self):
        candidates = [candidate("A" * 20, y=9), candidate("C" * 20, non_y=1),
                      candidate("G" * 20), candidate("T" * 20)]
        rows = pipeline.evaluate(candidates, {"G" * 20: annotation()}, {}, DEFAULTS)
        decisions = {r["guide"]: r for r in rows}
        self.assertEqual(decisions["A" * 20]["decision"], "rejected_exact_match_filters")
        self.assertEqual(decisions["C" * 20]["decision"], "rejected_exact_match_filters")
        self.assertEqual(decisions["G" * 20]["decision"], "retained")
        self.assertEqual(decisions["T" * 20]["decision"], "unassessed_missing_mismatch_counts")
        self.assertIsNone(decisions["T" * 20]["non_y_mm1"])

    def test_partial_annotation_and_known_failure(self):
        rows = pipeline.evaluate([candidate("A" * 20), candidate("C" * 20)],
                                 {"A" * 20: annotation(None, 0), "C" * 20: annotation(None, 1)}, {}, DEFAULTS)
        self.assertEqual(rows[0]["decision"], "unassessed_missing_mismatch_counts")
        self.assertEqual(rows[1]["decision"], "rejected_mismatch_filters")

    def test_thresholds_are_effective_and_labels_do_not_force_retention(self):
        c = [candidate("A" * 20, y=11)]
        a = {"A" * 20: annotation(mm1=1)}
        labels = {"A" * 20: {"guide_id": "historic_anchor"}}
        self.assertEqual(pipeline.evaluate(c, a, labels, DEFAULTS)[0]["decision"], "rejected_mismatch_filters")
        relaxed = {**DEFAULTS, "maximum_non_y_1_mismatch_sites": 1}
        self.assertEqual(pipeline.evaluate(c, a, labels, relaxed)[0]["decision"], "retained")
        stricter = {**relaxed, "minimum_exact_y_matches": 12}
        self.assertEqual(pipeline.evaluate(c, a, labels, stricter)[0]["decision"], "rejected_exact_match_filters")

    def test_ranking_is_deterministic_and_cfd_is_not_a_filter(self):
        c = [candidate("T" * 20, 11), candidate("G" * 20, 10), candidate("A" * 20, 10)]
        a = {r["guide"]: annotation(cfd=9999) for r in c}
        x = pipeline.evaluate(c, a, {}, DEFAULTS)
        y = pipeline.evaluate(c[::-1], a, {}, DEFAULTS)
        self.assertEqual(x, y)
        self.assertEqual([r["guide"] for r in x], ["T" * 20, "A" * 20, "G" * 20])
        self.assertEqual([r["rank"] for r in x], [1, 2, 3])
        self.assertTrue(all(r["legacy_non_y_cfd_sum"] == 9999 for r in x))

    def test_invalid_counts_and_scope_are_rejected(self):
        for value in (-1, 1.5, "inf"):
            with self.assertRaises(ValueError):
                pipeline.number(value, "count", integer=True)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "annotations.csv"
            columns = ["guide", "mismatch_1", "mismatch_2", "chrY_rows"]
            pipeline.write_csv(path, [dict(guide="A" * 20, mismatch_1=0, mismatch_2=0, chrY_rows=1)], columns)
            with self.assertRaisesRegex(ValueError, "chrY"):
                pipeline.load_annotations(path)
            row = dict(guide="A" * 20, mismatch_1=0, mismatch_2=0, chrY_rows=0)
            pipeline.write_csv(path, [row, row], columns)
            with self.assertRaisesRegex(ValueError, "duplicate"):
                pipeline.load_annotations(path)

    def test_end_to_end_partial_coverage_and_provenance(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            pipeline.write_csv(root / "candidates.csv", [candidate("A" * 20), candidate("C" * 20)],
                               ["guide"] + pipeline.COUNT_COLS + ["sumY", "sumOther"])
            pipeline.write_csv(root / "annotations.csv", [{"guide": "A" * 20, "mismatch_1": 0, "mismatch_2": 0}],
                               ["guide", "mismatch_1", "mismatch_2"])
            pipeline.write_csv(root / "labels.csv", [], ["guide", "guide_id", "region"])
            pipeline.write_csv(root / "historical.csv", [], ["guide"])
            config = {"pipeline_version": pipeline.VERSION,
                      "inputs": {"candidates": "candidates.csv", "non_y_mismatch_annotations": "annotations.csv",
                                 "guide_labels": "labels.csv", "historical_shortlist_for_comparison_only": "historical.csv"},
                      "filters": DEFAULTS, "ranking": pipeline.RANKING, "annotation_scope": "non_Y",
                      "annotation_provenance": "synthetic fixture", "score_interpretation": "annotation only",
                      "reference_description": "synthetic fixture", "output_root": "runs"}
            config_path = root / "config.json"
            config_path.write_text(json.dumps(config), encoding="utf-8")
            with redirect_stdout(StringIO()), redirect_stderr(StringIO()):
                code = pipeline.run(config_path, root / "strict")
                allowed = pipeline.run(config_path, root / "partial", allow_partial=True)
            self.assertEqual(code, 2)
            self.assertEqual(allowed, 0)
            summary = json.loads((root / "partial/run_manifest.json").read_text())
            self.assertEqual(summary["counts"]["retained"], 1)
            self.assertEqual(summary["counts"]["unassessed"], 1)
            self.assertEqual(summary["status"], "partial_mismatch_coverage")
            self.assertEqual(summary["inputs"]["candidates"]["sha256"], pipeline.sha256(root / "candidates.csv"))
            with (root / "partial/unassessed_guides.csv").open(newline="") as f:
                row = next(csv.DictReader(f))
            self.assertEqual(row["non_y_mm1"], "")
            with self.assertRaisesRegex(ValueError, "already exists"):
                pipeline.run(config_path, root / "partial", allow_partial=True)


if __name__ == "__main__":
    unittest.main()
