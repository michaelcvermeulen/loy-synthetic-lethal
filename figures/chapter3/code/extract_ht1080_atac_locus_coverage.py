"""Extract fragment-CPM coverage at the selected HT1080 ATAC loci (run under WSL).

Uses the SEPT_08_2026 per-library bigWigs (fragment coverage per million non-Y
fragments) and averages them into 50 bp bins. Requires pyBigWig, e.g.
/home/micha/.local/share/atacseq-sept08-2026-venv/bin/python.
Optional argument: a source_data subfolder name (default ht1080_atac_accessibility)
containing locus_definitions.csv; output is written beside it.
"""
from pathlib import Path
import csv, sys
import numpy as np
import pyBigWig

ROOT = Path(__file__).resolve().parents[3]
TRACKS = ROOT.parent / "ATACseq/SEPT_08_2026_REPROCESS/tracks"
DATA = ROOT / "figures/chapter3/source_data" / (sys.argv[1] if len(sys.argv) > 1 else "ht1080_atac_accessibility")
SAMPLES = ["CTRL_11", "CTRL_12", "CTRL_21", "LOY_11", "LOY_12", "LOY_21"]
BIN = 50

loci = list(csv.DictReader(open(DATA / "locus_definitions.csv")))
with open(DATA / "locus_coverage_50bp.csv", "w", newline="") as handle:
    out = csv.writer(handle)
    out.writerow(["peak_id", "sample", "bin_start", "bin_end", "mean_fragment_cpm"])
    for sample in SAMPLES:
        bw = pyBigWig.open(str(TRACKS / f"{sample}.fragment_cpm.bw"))
        for locus in loci:
            start, end = int(locus["window_start"]), int(locus["window_end"])
            values = np.nan_to_num(bw.values(locus["chrom"], start, end, numpy=True))
            for i in range(0, end - start, BIN):
                out.writerow([locus["peak_id"], sample, start + i, start + i + BIN,
                              f"{values[i:i + BIN].mean():.5f}"])
        bw.close()
print("done")
