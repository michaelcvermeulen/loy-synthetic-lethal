#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
NATIVE=/home/micha/codex-work/wgs-ypar-pilot-20261027
export PATH="$ROOT/runtime/env/bin:$PATH"
SAMPLE=A549-CTRL2
test -s "$ROOT/reference/INDEX_COMPLETE.txt"
test -s "$ROOT/qc/FASTQ_RECOVERY_COMPLETE.txt"
mkdir -p "$NATIVE/fastq" "$NATIVE/alignment" "$NATIVE/tmp"
cp -- "$ROOT/fastq/$SAMPLE.R1.fastq.gz" "$ROOT/fastq/$SAMPLE.R2.fastq.gz" "$NATIVE/fastq/"
test ! -e "$ROOT/alignment/$SAMPLE.YPARmasked.bam"
# Fixed batch size records original 16-thread default of 160 million input bases.
# Use 8 mapping workers plus sorting workers to keep this pilot within 12 CPUs.
bwa-mem2 mem -t 8 -K 160000000 -R '@RG\tID:A549-CTRL2\tSM:A549-CTRL2\tPL:ILLUMINA' \
    "$NATIVE/reference/chm13v2.0.YPARmasked.fa" \
    "$NATIVE/fastq/$SAMPLE.R1.fastq.gz" "$NATIVE/fastq/$SAMPLE.R2.fastq.gz" \
    2> "$ROOT/logs/$SAMPLE.alignment.log" |
    samtools sort -@ 3 -m 1G -T "$NATIVE/tmp/$SAMPLE.sort" -o "$NATIVE/alignment/$SAMPLE.YPARmasked.bam" -
samtools quickcheck -v "$NATIVE/alignment/$SAMPLE.YPARmasked.bam"
samtools index -@ 4 "$NATIVE/alignment/$SAMPLE.YPARmasked.bam"
cp -- "$NATIVE/alignment/$SAMPLE.YPARmasked.bam" "$NATIVE/alignment/$SAMPLE.YPARmasked.bam.bai" "$ROOT/alignment/"
samtools flagstat -@ 4 "$ROOT/alignment/$SAMPLE.YPARmasked.bam" > "$ROOT/qc/$SAMPLE.masked_flagstat.txt"
samtools view -H "$ROOT/alignment/$SAMPLE.YPARmasked.bam" > "$ROOT/qc/$SAMPLE.masked_header.sam"
python "$ROOT/pilot_qc.py" bam "$ROOT/alignment/$SAMPLE.YPARmasked.bam" "$ROOT/qc/$SAMPLE.masked_read_qc.json"
python - "$ROOT" <<'PY'
import json,sys
from pathlib import Path
p=Path(sys.argv[1])/'qc'
a=json.loads((p/'A549-CTRL2.original_read_qc.json').read_text())
b=json.loads((p/'A549-CTRL2.masked_read_qc.json').read_text())
assert a['read_fingerprint']==b['read_fingerprint'],'Read content changed during alignment'
(p/'REALIGNMENT_COMPLETE.txt').write_text('All primary read names, sequences and qualities preserved.\n')
PY
