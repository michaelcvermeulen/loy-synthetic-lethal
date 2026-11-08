from __future__ import annotations

import unittest
from pathlib import Path

import numpy as np

from spectral.a8_pipeline.fcs_io import load_fcs
from spectral.a8_pipeline.saturation import (
    decode_mask_matrix,
    decode_saturation,
    reinterpret_u32,
)


INPUT_ROOT = Path(r"A:\flow_jo_AUG28")


class SaturationTests(unittest.TestCase):
    def test_reinterprets_float_storage_bytes(self) -> None:
        bits = np.array([0, 1, 2**31, 2**32 - 1], dtype=np.uint32)
        stored = bits.view(np.float32)
        np.testing.assert_array_equal(reinterpret_u32(stored), bits)

    def test_crosses_mask_boundaries(self) -> None:
        masks = np.zeros((1, 3), dtype=np.uint32)
        masks[0, 0] |= np.uint32(1 << 31)
        masks[0, 1] |= np.uint32(1 << 0)
        masks[0, 2] |= np.uint32(1 << 5)
        decoded = decode_mask_matrix(masks.view(np.float32))
        self.assertEqual(np.flatnonzero(decoded[0]).tolist(), [31, 32, 69])

    def test_decoded_or_matches_saturated_for_every_source_event(self) -> None:
        for path in sorted(INPUT_ROOT.glob("*.fcs")):
            with self.subTest(path=path.name):
                data = load_fcs(path)
                flags = decode_saturation(data.events)
                expected = data.events["Saturated"].to_numpy() != 0
                np.testing.assert_array_equal(flags["any_saturated"].to_numpy(), expected)

    def test_mixed_spectral_and_imgb3_rates_match_export(self) -> None:
        expected = {
            "F11-LOY_GFP___CTRL_mCherry.fcs": (0.0040, 0.56087),
            "G11-CTRL_GFP__LOY_mCherry.fcs": (0.0093, 0.46795),
        }
        for filename, (spectral, imgb3) in expected.items():
            with self.subTest(filename=filename):
                flags = decode_saturation(load_fcs(INPUT_ROOT / filename).events)
                self.assertAlmostEqual(float(flags["spectral_saturated"].mean()), spectral, places=5)
                self.assertAlmostEqual(float(flags["imgb3_saturated"].mean()), imgb3, places=5)


if __name__ == "__main__":
    unittest.main()
