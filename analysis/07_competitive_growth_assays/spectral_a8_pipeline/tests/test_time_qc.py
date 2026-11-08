from __future__ import annotations

import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from spectral.a8_pipeline.fcs_io import load_fcs
from spectral.a8_pipeline.time_qc import Segment, detect_segments, robust_outliers, time_qc


INPUT_ROOT = Path(r"A:\flow_jo_AUG28")
F11 = INPUT_ROOT / "F11-LOY_GFP___CTRL_mCherry.fcs"
G11 = INPUT_ROOT / "G11-CTRL_GFP__LOY_mCherry.fcs"


class TimeQCTests(unittest.TestCase):
    def test_detects_large_gap(self) -> None:
        time = np.array([0.0, 0.1, 0.2, 10.2, 10.3])
        self.assertEqual([(s.start, s.stop) for s in detect_segments(time)], [(0, 3), (3, 5)])

    def test_detects_time_reset(self) -> None:
        time = np.array([0.0, 0.1, 0.2, 0.0, 0.1])
        self.assertEqual([(s.start, s.stop) for s in detect_segments(time)], [(0, 3), (3, 5)])

    def test_flags_zero_mad_rate_outlier(self) -> None:
        rates = np.array([200.0, 200.0, 30.0, 200.0])
        self.assertEqual(robust_outliers(rates, z_limit=3.5).tolist(), [False, False, True, False])

    def test_real_acquisitions_have_expected_segments(self) -> None:
        f11 = load_fcs(F11, required=("Time",))
        g11 = load_fcs(G11, required=("Time",))
        self.assertEqual([s.size for s in detect_segments(f11.events["Time"].to_numpy())], [15000])
        self.assertEqual([s.size for s in detect_segments(g11.events["Time"].to_numpy())], [10000, 10000])

    def test_time_qc_marks_slow_event_bin(self) -> None:
        increments = np.r_[np.full(200, 0.005), np.full(200, 0.005), np.full(200, 0.05), np.full(200, 0.005)]
        events = pd.DataFrame(
            {
                "Time": np.cumsum(increments),
                "FSC-A": np.full(800, 100.0),
                "SSC (Violet)-A": np.full(800, 50.0),
                "eGFP-A": np.full(800, 10.0),
                "mCherry-A": np.full(800, 20.0),
                "Autofluorescence_001-A": np.full(800, 5.0),
            }
        )
        result = time_qc(events, [Segment(0, 800, "Run 1")], min_events=200)
        self.assertEqual(result.bins["bad"].tolist(), [False, False, True, False])
        self.assertEqual(result.good_events.sum(), 600)
        self.assertFalse(result.good_events[400:600].any())


if __name__ == "__main__":
    unittest.main()
