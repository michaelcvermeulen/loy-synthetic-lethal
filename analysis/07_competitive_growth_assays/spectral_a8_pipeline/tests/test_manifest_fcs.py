from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from spectral.a8_pipeline.fcs_io import load_fcs
from spectral.a8_pipeline.manifest import SourceRecord, verify_manifest


INPUT_ROOT = Path(r"A:\flow_jo_AUG28")
F11 = INPUT_ROOT / "F11-LOY_GFP___CTRL_mCherry.fcs"


class ManifestAndFCSTests(unittest.TestCase):
    def test_hash_verification_detects_changed_file(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            path = Path(td, "x.bin")
            path.write_bytes(b"first")
            record = SourceRecord.from_path(path)
            path.write_bytes(b"second")
            with self.assertRaisesRegex(RuntimeError, "source changed"):
                verify_manifest([record])

    def test_load_fcs_requires_named_parameters(self) -> None:
        with self.assertRaisesRegex(KeyError, "missing required parameters"):
            load_fcs(F11, required=("not-a-channel",))

    def test_mixed_fcs_exposes_hard_coded_reporters(self) -> None:
        data = load_fcs(F11, required=("eGFP-A", "mCherry-A", "Time"))
        self.assertEqual(data.events.shape, (15000, 253))
        self.assertEqual(data.matrix_name, "mcherry_GFP_matrix")


if __name__ == "__main__":
    unittest.main()
