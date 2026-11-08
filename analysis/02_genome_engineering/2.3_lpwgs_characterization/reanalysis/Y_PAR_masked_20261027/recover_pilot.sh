#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
export PATH="$ROOT/runtime/env/bin:$PATH"
SAMPLE=A549-CTRL2
BAM=/mnt/a/Dropbox/LOY_paralogs/projects/WGS/t2t_bam/$SAMPLE.bam
samtools quickcheck -v "$BAM"
samtools view -H "$BAM" > "$ROOT/qc/$SAMPLE.original_header.sam"
samtools flagstat -@ 4 "$BAM" > "$ROOT/qc/$SAMPLE.original_flagstat.txt"
python "$ROOT/pilot_qc.py" bam "$BAM" "$ROOT/qc/$SAMPLE.original_read_qc.json"
if [[ -e "$ROOT/fastq/$SAMPLE.R1.fastq.gz" ]]; then echo 'FASTQ output already exists; refusing overwrite' >&2; exit 1; fi
samtools collate -O -u -@ 4 "$BAM" "$ROOT/tmp/$SAMPLE.collate" |
    samtools fastq -@ 4 -c 1 -n -F 0x900 \
        -1 "$ROOT/fastq/$SAMPLE.R1.fastq.gz" -2 "$ROOT/fastq/$SAMPLE.R2.fastq.gz" \
        -s "$ROOT/fastq/$SAMPLE.singletons.fastq.gz" -0 "$ROOT/fastq/$SAMPLE.other.fastq.gz" - \
        2> "$ROOT/logs/$SAMPLE.fastq.log"
test "$(gzip -cd "$ROOT/fastq/$SAMPLE.singletons.fastq.gz" | wc -c)" -eq 0
test "$(gzip -cd "$ROOT/fastq/$SAMPLE.other.fastq.gz" | wc -c)" -eq 0
python "$ROOT/pilot_qc.py" fastq "$ROOT/fastq/$SAMPLE.R1.fastq.gz" "$ROOT/fastq/$SAMPLE.R2.fastq.gz" \
    "$ROOT/qc/$SAMPLE.original_read_qc.json" "$ROOT/qc/$SAMPLE.fastq_validation.json"
date -Is > "$ROOT/qc/FASTQ_RECOVERY_COMPLETE.txt"
