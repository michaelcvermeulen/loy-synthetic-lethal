from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from spectral.a8_pipeline.gating import (
    points_in_polygon,
    pulse_relationship_mask,
    pulse_singlet_mask,
    scale_polygon,
)
from spectral.a8_pipeline.svg import write_scatter_density_svg


class GatingTests(unittest.TestCase):
    def test_polygon_inclusion_and_scaling(self) -> None:
        polygon = np.array([[0, 0], [2, 0], [2, 2], [0, 2]], dtype=float)
        points = np.array([[1, 1], [3, 1]], dtype=float)
        self.assertEqual(points_in_polygon(points, polygon).tolist(), [True, False])
        expanded = scale_polygon(polygon, 1.05)
        np.testing.assert_allclose(expanded.mean(axis=0), polygon.mean(axis=0))

    def test_pulse_gate_rejects_line_outlier_and_nonpositive(self) -> None:
        area = np.array([100, 200, 300, 400, 500, 600, 0], dtype=float)
        height = np.array([10, 20, 30, 40, 50, 300, 5], dtype=float)
        mask, fit = pulse_relationship_mask(area, height, mad_limit=3.5)
        self.assertEqual(mask.tolist(), [True, True, True, True, True, False, False])
        self.assertAlmostEqual(fit.slope, 1.0, places=2)

    def test_primary_singlet_gate_intersects_fsc_and_ssc(self) -> None:
        events = pd.DataFrame(
            {
                "FSC-A": [100, 200, 300, 400, 500, 600],
                "FSC-H": [10, 20, 30, 40, 50, 300],
                "SSC (Violet)-A": [100, 200, 300, 400, 500, 600],
                "SSC (Violet)-H": [10, 20, 30, 40, 200, 60],
            }
        )
        mask, fits = pulse_singlet_mask(events, mad_limit=3.5)
        self.assertEqual(mask.tolist(), [True, True, True, True, False, False])
        self.assertEqual(set(fits), {"fsc", "ssc"})

    def test_scatter_svg_contains_density_and_polygon(self) -> None:
        events = pd.DataFrame({"FSC-A": [10, 20, 30], "SSC (Violet)-A": [5, 10, 20]})
        polygon = np.array([[5, 2], [40, 2], [40, 30], [5, 30]], dtype=float)
        with tempfile.TemporaryDirectory() as td:
            output = Path(td, "scatter.svg")
            write_scatter_density_svg(events, polygon, output, title="Test")
            text = output.read_text(encoding="utf-8")
            self.assertIn("<svg", text)
            self.assertIn("<polygon", text)
            self.assertGreater(output.stat().st_size, 500)


if __name__ == "__main__":
    unittest.main()
