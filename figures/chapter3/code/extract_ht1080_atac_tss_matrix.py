"""TSS-centred coverage matrices for the HT1080 ATAC QC heatmaps (run under WSL).

Uses the SEPT_08_2026 per-library fragment-CPM bigWigs and the same TSS set and
filters as the pipeline's TSS-enrichment score (02_prepare.py): knownGene TSSs,
no chrY, unique (chrom, position, strand), window inside the chromosome and not
overlapping the ENCODE hg38 blacklist. The window here is +/-3 kb (50 bp bins);
minus-strand windows are reversed so upstream is always on the left.
Requires pyBigWig: /home/micha/.local/share/atacseq-sept08-2026-venv/bin/python
"""
from pathlib import Path
import csv, json, os
import numpy as np
import pyBigWig

ROOT = Path(__file__).resolve().parents[3]
RUN = ROOT.parent / "ATACseq/SEPT_08_2026_REPROCESS"
OUT = ROOT / "figures/chapter3/source_data/ht1080_atac_accessibility"
SAMPLES = ["CTRL_11", "CTRL_12", "CTRL_21", "LOY_11", "LOY_12", "LOY_21"]
HALF, BIN = 3000, 50

blacklist = {}
for chrom, start, end, *_ in csv.reader(open(RUN / "config/hg38-blacklist.v2.bed"), delimiter="\t"):
    blacklist.setdefault(chrom, []).append((int(start), int(end)))

# TRACKS_DIR may point at a local copy of RUN/tracks (random reads over /mnt are slow).
track_dir = Path(os.environ.get("TRACKS_DIR", RUN / "tracks"))
tracks = {s: pyBigWig.open(str(track_dir / f"{s}.fragment_cpm.bw")) for s in SAMPLES}
sizes = tracks[SAMPLES[0]].chroms()
sites, used = [], set()
for chrom, start, end, gene, score, strand in csv.reader(open(RUN / "config/tss_knownGene.bed"), delimiter="\t"):
    centre, key = int(start), (chrom, int(start), strand)
    if chrom == "chrY" or chrom not in sizes or key in used:
        continue
    if centre < HALF or centre + HALF > sizes[chrom]:
        continue
    if any(s < centre + HALF and e > centre - HALF for s, e in blacklist.get(chrom, [])):
        continue
    used.add(key); sites.append((chrom, centre, strand, gene))

nb = 2 * HALF // BIN
mats = {}
for s in SAMPLES:
    m = np.zeros((len(sites), nb), dtype=np.float32)
    for i, (chrom, centre, strand, _) in enumerate(sites):
        # Base-level values averaged into bins; uncovered bases (NaN) count as zero coverage.
        v = np.nan_to_num(tracks[s].values(chrom, centre - HALF, centre + HALF, numpy=True))
        row = v.reshape(nb, BIN).mean(axis=1).astype(np.float32)
        m[i] = row if strand == "+" else row[::-1]
    mats[s] = m
    tracks[s].close()
np.savez_compressed(OUT / "tss_heatmap_matrix.npz", **mats,
                    sites=np.array([f"{c}:{p}:{st}:{g}" for c, p, st, g in sites]))
json.dump({"tss_sites": len(sites), "half_window_bp": HALF, "bin_bp": BIN,
           "signal": "fragment coverage per million non-Y fragments; base-level mean per 50 bp bin, uncovered bases = 0"},
          open(OUT / "tss_heatmap_matrix.json", "w"), indent=1)
print(len(sites), "TSS sites")
