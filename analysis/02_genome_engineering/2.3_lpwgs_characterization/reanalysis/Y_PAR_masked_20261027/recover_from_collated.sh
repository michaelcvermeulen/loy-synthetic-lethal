#!/usr/bin/env bash
# Resume after a completed name-collated BAM was produced during the pilot.
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
export PATH="$ROOT/runtime/env/bin:$PATH"
SAMPLE=A549-CTRL2
BAM="$ROOT/tmp/$SAMPLE.collate.bam"
samtools quickcheck -v "$BAM"
test -s "$ROOT/qc/$SAMPLE.original_read_qc.json"
test ! -e "$ROOT/fastq/$SAMPLE.R1.fastq.gz"
samtools fastq -@ 4 -c 1 -n -F 0x900 \
    -1 "$ROOT/fastq/$SAMPLE.R1.fastq.gz" -2 "$ROOT/fastq/$SAMPLE.R2.fastq.gz" \
    -s "$ROOT/fastq/$SAMPLE.singletons.fastq.gz" -0 "$ROOT/fastq/$SAMPLE.other.fastq.gz" "$BAM" \
    2> "$ROOT/logs/$SAMPLE.fastq_resumed.log"
test "$(gzip -cd "$ROOT/fastq/$SAMPLE.singletons.fastq.gz" | wc -c)" -eq 0
test "$(gzip -cd "$ROOT/fastq/$SAMPLE.other.fastq.gz" | wc -c)" -eq 0
python "$ROOT/pilot_qc.py" fastq "$ROOT/fastq/$SAMPLE.R1.fastq.gz" "$ROOT/fastq/$SAMPLE.R2.fastq.gz" \
    "$ROOT/qc/$SAMPLE.original_read_qc.json" "$ROOT/qc/$SAMPLE.fastq_validation.json"
date -Is > "$ROOT/qc/FASTQ_RECOVERY_COMPLETE.txt"
