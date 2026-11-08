from __future__ import annotations

import csv
import tempfile
import unittest
from pathlib import Path

from PIL import Image

from spectral.a8_pipeline.reporter_contact_sheets import (
    format_panel_title,
    generate_contact_sheets,
    inventory_reporter_samples,
)


class ReporterContactSheetTests(unittest.TestCase):
    def test_long_sample_titles_wrap_inside_one_panel(self) -> None:
        """Long adjacent sample names must not collide across panel boundaries."""
        title = format_panel_title(
            "B02-LOY_GFP_CTRL_CHERRY_rep1_run1_001"
        )
        self.assertGreaterEqual(len(title.splitlines()), 2)
        self.assertTrue(all(len(line) <= 42 for line in title.splitlines()))

    def test_inventory_assigns_all_67_fcs_files_to_one_sheet(self) -> None:
        """An incorrect group rule would silently omit or duplicate a sample panel."""
        repo = Path(__file__).resolve().parents[3]

        samples = inventory_reporter_samples(repo)

        self.assertEqual(len(samples), 67)
        self.assertEqual(len({sample.path.resolve() for sample in samples}), 67)
        counts: dict[str, int] = {}
        for sample in samples:
            counts[sample.sheet] = counts.get(sample.sheet, 0) + 1
        self.assertEqual(
            counts,
            {
                "run1_d0": 5,
                "run1_d4_controls": 6,
                "run1_d4_ratio_standards": 20,
                "run1_d4_experimental": 12,
                "run1_d7": 12,
                "run2_d0": 12,
            },
        )

    def test_real_generation_writes_six_nonblank_sheets_and_complete_manifest(self) -> None:
        """A successful command must produce viewable panels and an auditable manifest."""
        repo = Path(__file__).resolve().parents[3]
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "plots"

            generated = generate_contact_sheets(repo, output)

            self.assertEqual(len(generated), 6)
            with (output / "plot_manifest.csv").open(newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 67)
            self.assertEqual(len({row["source_path"] for row in rows}), 67)
            self.assertEqual(sum(row["status"] == "reporter_channels_unavailable" for row in rows), 3)
            for path in generated:
                self.assertGreater(path.stat().st_size, 10_000)
                with Image.open(path) as image:
                    self.assertGreaterEqual(image.width, 1800)
                    self.assertGreaterEqual(image.height, 900)
                    self.assertNotEqual(image.getextrema(), ((255, 255), (255, 255), (255, 255), (255, 255)))


if __name__ == "__main__":
    unittest.main()
