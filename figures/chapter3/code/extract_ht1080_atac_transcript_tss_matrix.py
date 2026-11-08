"""Coverage at transcript-level TSSs of the 13,300 integration-universe genes (run under WSL).

Reads tss_transcript_sites.csv (prepare_ht1080_atac_transcript_tss.R) and the SEPT_08_2026
fragment-CPM bigWigs. Filters as extract_ht1080_atac_tss_matrix.py: window inside the chromosome and
not overlapping the ENCODE hg38 blacklist. +/-3 kb, 50 bp bins, base-level mean with uncovered bases
counted as zero; minus-strand windows reversed so upstream is on the left.
Each gene keeps one TSS: the one with the highest central (+/-100 bp) coverage pooled over all six
libraries (condition-independent choice). Full profiles are saved for the chosen TSSs only; the
central coverage of every candidate TSS is saved as CSV.
Set TRACKS_DIR to a local copy of the tracks (random reads over /mnt are slow).
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
track_dir = Path(os.environ.get("TRACKS_DIR", RUN / "tracks"))
tracks = {s: pyBigWig.open(str(track_dir / f"{s}.fragment_cpm.bw")) for s in SAMPLES}
sizes = tracks[SAMPLES[0]].chroms()

sites, dropped = [], 0
for r in csv.DictReader(open(OUT / "tss_transcript_sites.csv")):
    chrom, centre = r["chrom"], int(r["tss"])
    if chrom not in sizes or centre < HALF or centre + HALF > sizes[chrom] or \
            any(s < centre + HALF and e > centre - HALF for s, e in blacklist.get(chrom, [])):
        dropped += 1
        continue
    sites.append((r["entrez_id"], chrom, centre, r["strand"]))

nb = 2 * HALF // BIN
c0, c1 = nb // 2 - 2, nb // 2 + 2
full = {s: np.zeros((len(sites), nb), dtype=np.float32) for s in SAMPLES}
for s in SAMPLES:
    for i, (_, chrom, centre, strand) in enumerate(sites):
        v = np.nan_to_num(tracks[s].values(chrom, centre - HALF, centre + HALF, numpy=True))
        row = v.reshape(nb, BIN).mean(axis=1)
        full[s][i] = row if strand == "+" else row[::-1]
    tracks[s].close()
    print("done", s, flush=True)

central = np.stack([full[s][:, c0:c1].mean(1) for s in SAMPLES], axis=1)
pooled = central.mean(1)
best = {}
for i, (gene, *_rest) in enumerate(sites):
    if gene not in best or pooled[i] > pooled[best[gene]]:
        best[gene] = i
chosen = np.array(sorted(best.values()))
with open(OUT / "tss_transcript_central.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["entrez_id", "chrom", "tss", "strand", *SAMPLES, "pooled_centre", "chosen"])
    chosen_set = set(chosen.tolist())
    for i, (gene, chrom, centre, strand) in enumerate(sites):
        w.writerow([gene, chrom, centre, strand, *np.round(central[i], 5), round(float(pooled[i]), 5), int(i in chosen_set)])
np.savez_compressed(OUT / "tss_transcript_matrix.npz", **{s: full[s][chosen] for s in SAMPLES},
                    sites=np.array([f"{sites[i][1]}:{sites[i][2]}:{sites[i][3]}:{sites[i][0]}" for i in chosen]))
json.dump({"candidate_tss": len(sites), "dropped_edge_or_blacklist": dropped, "genes": len(chosen),
           "half_window_bp": HALF, "bin_bp": BIN, "choice": "max pooled central (+/-100 bp) coverage over six libraries",
           "signal": "fragment coverage per million non-Y fragments; base-level mean per 50 bp bin, uncovered = 0"},
          open(OUT / "tss_transcript_matrix.json", "w"), indent=1)
print(len(sites), "candidate TSSs;", len(chosen), "genes;", dropped, "dropped")
