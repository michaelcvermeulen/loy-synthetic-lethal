from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import numpy as np
from PIL import Image

from spectral.a8_pipeline.images import (
    build_image_manifest,
    parse_event_number,
    read_tiff_channels,
    render_montage,
)

ROOT = Path(r"A:\flow_jo_AUG28")
F11_IMAGES = ROOT / "F11-LOY_GFP___CTRL_mCherry_images"
G11_IMAGES = ROOT / "G11-CTRL_GFP__LOY_mCherry_images"


def _write_stack(path: Path, frames: list[np.ndarray]) -> None:
    images = [Image.fromarray(frame.astype(np.float32), mode="F") for frame in frames]
    images[0].save(path, save_all=True, append_images=images[1:])


class ImageMappingTests(unittest.TestCase):
    def test_parses_zero_based_suffix_literally(self) -> None:
        path = Path("F11-LOY_GFP___CTRL_mCherry_00000007.tiff")
        self.assertEqual(parse_event_number(path), 7)

    def test_rejects_duplicate_event(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "a").mkdir()
            Path(root, "sample_00000002.tiff").touch()
            Path(root, "a", "copy_00000002.tiff").touch()
            with self.assertRaisesRegex(ValueError, "duplicate image event"):
                build_image_manifest(root, event_count=10)

    def test_rejects_out_of_range_event(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            Path(root, "sample_00000010.tiff").touch()
            with self.assertRaisesRegex(ValueError, "out of range"):
                build_image_manifest(root, event_count=10)

    def test_manifest_is_one_based_and_marks_missing_images(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            Path(root, "sample_00000000.tiff").touch()
            Path(root, "sample_00000003.tiff").touch()
            manifest = build_image_manifest(root, event_count=5)
            self.assertEqual(len(manifest), 5)
            self.assertEqual(manifest["event_number"].tolist(), [1, 2, 3, 4, 5])
            self.assertEqual(manifest["available"].tolist(), [True, False, False, True, False])
            self.assertEqual(manifest.loc[0, "image_index"], 0)
            self.assertIsNone(manifest.loc[1, "path"])


class TiffReadingTests(unittest.TestCase):
    def test_reads_every_frame_distinctly(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary, "stack_00000000.tiff")
            frames = [np.full((4, 5), value, dtype=np.float32) for value in (1.0, 2.0, 3.0)]
            _write_stack(path, frames)
            stack = read_tiff_channels(path)
            self.assertEqual(stack.shape, (3, 4, 5))
            self.assertEqual([float(frame.mean()) for frame in stack], [1.0, 2.0, 3.0])

    def test_real_tiff_has_six_frames(self) -> None:
        path = G11_IMAGES / "00000000" / "G11-CTRL_GFP__LOY_mCherry_00000000.tiff"
        stack = read_tiff_channels(path)
        self.assertEqual(stack.shape[0], 6)
        self.assertGreater(len({float(frame.mean()) for frame in stack}), 1)


class MontageTests(unittest.TestCase):
    def test_montage_pads_variable_heights_and_writes_png(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = []
            for index, height in enumerate((4, 7)):
                path = root / f"cell_{index:08d}.tiff"
                _write_stack(path, [np.full((height, 5), 1.0, dtype=np.float32)] * 6)
                paths.append(path)
            output = root / "montage.png"
            render_montage(
                {"Example": [(1, paths[0]), (2, paths[1])]},
                output,
                channel_map={"gray": 0, "green": 3, "red": 5},
            )
            self.assertTrue(output.exists())
            self.assertGreater(output.stat().st_size, 0)


class RealCoverageTests(unittest.TestCase):
    def test_reported_coverage_matches_export(self) -> None:
        f11 = build_image_manifest(F11_IMAGES, event_count=15000)
        g11 = build_image_manifest(G11_IMAGES, event_count=20000)
        self.assertEqual(int(f11["available"].sum()), 8776)
        self.assertEqual(int(g11["available"].sum()), 20000)
        self.assertEqual(len(f11), 15000)
        self.assertEqual(len(g11), 20000)


if __name__ == "__main__":
    unittest.main()
