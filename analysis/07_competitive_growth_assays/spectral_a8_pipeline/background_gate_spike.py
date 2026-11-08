"""Throwaway comparison of 99th-percentile unstained gates with LDA calls."""
import json
from pathlib import Path
import numpy as np
from spectral.a8_pipeline import competition_report as c

repo = Path(__file__).parents[2]
d4 = repo / "spectral/20260901_125230_Mike_Run1/Mike_Run1"
gates = c.load_cef_gates(d4 / "20260901_125230_Mike_Run1.cef")
files = c._file_map(d4 / "1")
unstained = np.vstack([c._primary_features(files[n], gates)[0] for n in ("A01-unstained_CTRL", "A02-unstained_LOY")])
thresholds = c.background_gate_thresholds(unstained, 99)
results = json.loads((repo / "analysis/spectral_competition_20260901/results.json").read_text())
print("pooled thresholds", thresholds.tolist())
print("sample | LDA GFP-like | gate GFP+ | gate mCherry+ | double+")
for name in sorted(files):
    if not name.startswith(("B", "C", "D", "E")):
        continue
    features = c._primary_features(files[name], gates)[0]
    calls = c.apply_background_gate(features, thresholds)
    panel = next(p for values in results["all_sample_plots"]["groups"].values() for p in values if p["sample"] == name)
    print(name, f"{100*panel['gfp_fraction']:.1f}%", f"{100*np.mean(calls==1):.1f}%", f"{100*np.mean(calls==2):.1f}%", f"{100*np.mean(calls==3):.1f}%", sep=" | ")
