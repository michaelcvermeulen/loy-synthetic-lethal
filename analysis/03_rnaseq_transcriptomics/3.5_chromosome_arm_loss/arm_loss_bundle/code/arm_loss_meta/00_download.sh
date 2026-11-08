#!/usr/bin/env bash
# Arm-loss meta-analysis: fetch the deposited count sets + hg38 cytoband table.
# Everything is public; nothing is realigned. Run from the arm_loss_meta/ root.
set -euo pipefail
W="${1:-$(pwd)}/data"
mkdir -p "$W"/{GSE222379,GSE336576,GSE150686}
GEO=https://ftp.ncbi.nlm.nih.gov/geo/series

# GSE222379  ReDACT engineered 1q/7p/8q loss (Girish et al., Science 2023) - 62 per-sample counts
( cd "$W/GSE222379" && curl -sSfO "$GEO/GSE222nnn/GSE222379/suppl/GSE222379_RAW.tar" \
    && tar xf GSE222379_RAW.tar && rm GSE222379_RAW.tar )
# GSE336576  RPTEC 3p-loss clones, day 18 + day 90 (Dahiya, van Belzen et al., bioRxiv 2026) - one matrix
( cd "$W/GSE336576" && curl -sSfO "$GEO/GSE336nnn/GSE336576/suppl/GSE336576_raw_counts_GEO.tsv.gz" )
# GSE150686  RPE1 monosomies + inducible p53 (Chunduri et al., Nat Commun 2021) - 33 per-sample counts
( cd "$W/GSE150686" && curl -sSfO "$GEO/GSE150nnn/GSE150686/suppl/GSE150686_RAW.tar" \
    && tar xf GSE150686_RAW.tar && rm GSE150686_RAW.tar \
    && curl -sSfO "$GEO/GSE150nnn/GSE150686/matrix/GSE150686_series_matrix.txt.gz" )
# hg38 cytoband (centromere = acen bands) for arm assignment
( cd "$W" && curl -sSfO https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz )

for g in GSE222379 GSE336576 GSE150686; do echo "$g: $(ls "$W/$g" | grep -c '\.gz$') gz files"; done
# expected: GSE222379 62, GSE336576 1, GSE150686 34 (33 counts + series matrix)
