"""Coverage matrices around differential and reference regions (run under WSL).

Reads region_heatmap_regions.csv (from prepare_ht1080_atac_region_heatmaps.py) and
the SEPT_08_2026 fragment-CPM bigWigs; +/-2 kb around each region centre in 50 bp
bins, base-level mean with uncovered bases counted as zero. Set TRACKS_DIR to a
local copy of the tracks (random reads over /mnt are slow).
"""
from pathlib import Path
import csv, json, os
import numpy as np
import pyBigWig

ROOT = Path(__file__).resolve().parents[3]
RUN = ROOT.parent / "ATACseq/SEPT_08_2026_REPROCESS"
DATA = ROOT / "figures/chapter3/source_data/ht1080_atac_accessibility"
SAMPLES = ["CTRL_11", "CTRL_12", "CTRL_21", "LOY_11", "LOY_12", "LOY_21"]
HALF, BIN = 2000, 50

regions = list(csv.DictReader(open(DATA / "region_heatmap_regions.csv")))
track_dir = Path(os.environ.get("TRACKS_DIR", RUN / "tracks"))
nb = 2 * HALF // BIN
mats = {}
for s in SAMPLES:
    bw = pyBigWig.open(str(track_dir / f"{s}.fragment_cpm.bw"))
    m = np.zeros((len(regions), nb), dtype=np.float32)
    for i, r in enumerate(regions):
        c = int(r["centre"])
        v = np.nan_to_num(bw.values(r["chrom"], c - HALF, c + HALF, numpy=True))
        m[i] = v.reshape(nb, BIN).mean(axis=1)
    mats[s] = m
    bw.close()
np.savez_compressed(DATA / "region_heatmap_matrix.npz", **mats, peak_id=np.array([r["peak_id"] for r in regions]))
json.dump({"regions": len(regions), "half_window_bp": HALF, "bin_bp": BIN,
           "signal": "fragment coverage per million non-Y fragments; base-level mean per bin, uncovered = 0"},
          open(DATA / "region_heatmap_matrix.json", "w"), indent=1)
print(len(regions), "regions")
