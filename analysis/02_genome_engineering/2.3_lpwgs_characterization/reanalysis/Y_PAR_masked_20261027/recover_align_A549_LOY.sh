#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
NATIVE=/home/micha/codex-work/wgs-ypar-pilot-20261027
export PATH="$ROOT/runtime/env/bin:$PATH"
SAMPLE=A549-LOY
BAM=/mnt/a/Dropbox/LOY_paralogs/projects/WGS/t2t_bam/$SAMPLE.bam
test -s "$ROOT/reference/INDEX_COMPLETE.txt"
test -s "$NATIVE/reference/chm13v2.0.YPARmasked.fa.bwt.2bit.64"
test ! -e "$ROOT/alignment/$SAMPLE.YPARmasked.bam"
test ! -e "$ROOT/fastq/$SAMPLE.R1.fastq.gz"
mkdir -p "$NATIVE/fastq" "$NATIVE/alignment" "$NATIVE/tmp"
sha256sum "$BAM" > "$ROOT/qc/$SAMPLE.original.sha256"
samtools quickcheck -v "$BAM"
samtools view -H "$BAM" > "$ROOT/qc/$SAMPLE.original_header.sam"
samtools flagstat -@ 4 "$BAM" > "$ROOT/qc/$SAMPLE.original_flagstat.txt"
echo "$(date -Is) Checking original primary reads"
python "$ROOT/pilot_qc.py" bam "$BAM" "$ROOT/qc/$SAMPLE.original_read_qc.json"
echo "$(date -Is) Recovering complete paired FASTQs"
samtools collate -O -u -@ 4 "$BAM" "$NATIVE/tmp/$SAMPLE.collate" |
    samtools fastq -@ 4 -c 1 -n -F 0x900 \
        -1 "$ROOT/fastq/$SAMPLE.R1.fastq.gz" -2 "$ROOT/fastq/$SAMPLE.R2.fastq.gz" \
        -s "$ROOT/fastq/$SAMPLE.singletons.fastq.gz" -0 "$ROOT/fastq/$SAMPLE.other.fastq.gz" - \
        2> "$ROOT/logs/$SAMPLE.fastq.log"
test "$(gzip -cd "$ROOT/fastq/$SAMPLE.singletons.fastq.gz" | wc -c)" -eq 0
test "$(gzip -cd "$ROOT/fastq/$SAMPLE.other.fastq.gz" | wc -c)" -eq 0
python "$ROOT/pilot_qc.py" fastq "$ROOT/fastq/$SAMPLE.R1.fastq.gz" "$ROOT/fastq/$SAMPLE.R2.fastq.gz" \
    "$ROOT/qc/$SAMPLE.original_read_qc.json" "$ROOT/qc/$SAMPLE.fastq_validation.json"
date -Is > "$ROOT/qc/$SAMPLE.FASTQ_RECOVERY_COMPLETE.txt"
cp -- "$ROOT/fastq/$SAMPLE.R1.fastq.gz" "$ROOT/fastq/$SAMPLE.R2.fastq.gz" "$NATIVE/fastq/"
echo "$(date -Is) Aligning against the same masked reference as A549-CTRL2"
bwa-mem2 mem -t 8 -K 160000000 -R '@RG\tID:A549-LOY\tSM:A549-LOY\tPL:ILLUMINA' \
    "$NATIVE/reference/chm13v2.0.YPARmasked.fa" \
    "$NATIVE/fastq/$SAMPLE.R1.fastq.gz" "$NATIVE/fastq/$SAMPLE.R2.fastq.gz" \
    2> "$ROOT/logs/$SAMPLE.alignment.log" |
    samtools sort -@ 3 -m 1G -T "$NATIVE/tmp/$SAMPLE.sort" -o "$NATIVE/alignment/$SAMPLE.YPARmasked.bam" -
samtools quickcheck -v "$NATIVE/alignment/$SAMPLE.YPARmasked.bam"
samtools index -@ 4 "$NATIVE/alignment/$SAMPLE.YPARmasked.bam"
cp -- "$NATIVE/alignment/$SAMPLE.YPARmasked.bam" "$NATIVE/alignment/$SAMPLE.YPARmasked.bam.bai" "$ROOT/alignment/"
samtools flagstat -@ 4 "$ROOT/alignment/$SAMPLE.YPARmasked.bam" > "$ROOT/qc/$SAMPLE.masked_flagstat.txt"
samtools view -H "$ROOT/alignment/$SAMPLE.YPARmasked.bam" > "$ROOT/qc/$SAMPLE.masked_header.sam"
echo "$(date -Is) Validating all realigned primary reads"
python "$ROOT/pilot_qc.py" bam "$ROOT/alignment/$SAMPLE.YPARmasked.bam" "$ROOT/qc/$SAMPLE.masked_read_qc.json"
python - "$ROOT" "$SAMPLE" <<'PY'
import json,sys
from pathlib import Path
p=Path(sys.argv[1])/'qc'; sample=sys.argv[2]
a=json.loads((p/f'{sample}.original_read_qc.json').read_text())
b=json.loads((p/f'{sample}.masked_read_qc.json').read_text())
assert a['read_fingerprint']==b['read_fingerprint'],'Read content changed during alignment'
(p/f'{sample}.REALIGNMENT_COMPLETE.txt').write_text('All primary read names, sequences and qualities preserved.\n')
PY
sha256sum -c "$ROOT/qc/$SAMPLE.original.sha256"
echo "$(date -Is) A549 LOY recovery and realignment complete"
