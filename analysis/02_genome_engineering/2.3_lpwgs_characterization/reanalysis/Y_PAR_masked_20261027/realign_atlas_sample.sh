#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
NATIVE=/home/micha/codex-work/wgs-ypar-pilot-20261027
export PATH="$ROOT/runtime/env/bin:$PATH"
SAMPLE=${1:?sample required}
case "$SAMPLE" in HT1080-sgCTRL-D0|HT1080-LOY-D0|LNCaP-sgCTRL|LNCaP-LOY-2X) ;; *) echo 'Sample outside authorised four-sample scope' >&2;exit 2;; esac
BAM=/mnt/a/Dropbox/LOY_paralogs/projects/WGS/t2t_bam/$SAMPLE.bam
test -s "$ROOT/reference/INDEX_COMPLETE.txt"
test -s "$NATIVE/reference/chm13v2.0.YPARmasked.fa.bwt.2bit.64"
mkdir -p "$NATIVE/fastq" "$NATIVE/alignment" "$NATIVE/tmp"
if [[ -s "$ROOT/qc/$SAMPLE.REALIGNMENT_COMPLETE.txt" ]]; then
    samtools quickcheck -v "$ROOT/alignment/$SAMPLE.YPARmasked.bam"
    sha256sum -c "$ROOT/qc/$SAMPLE.original.sha256"
    echo "$(date -Is) Already validated: $SAMPLE";exit 0
fi
if [[ ! -s "$ROOT/qc/$SAMPLE.ORIGINAL_QC_COMPLETE.txt" ]]; then
    sha256sum "$BAM" > "$ROOT/qc/$SAMPLE.original.sha256"
    samtools quickcheck -v "$BAM"
    samtools view -H "$BAM" > "$ROOT/qc/$SAMPLE.original_header.sam"
    samtools flagstat -@ 2 "$BAM" > "$ROOT/qc/$SAMPLE.original_flagstat.txt"
    echo "$(date -Is) Original read QC: $SAMPLE"
    python "$ROOT/pilot_qc.py" bam "$BAM" "$ROOT/qc/$SAMPLE.original_read_qc.json"
    date -Is > "$ROOT/qc/$SAMPLE.ORIGINAL_QC_COMPLETE.txt"
fi
if [[ ! -s "$ROOT/qc/$SAMPLE.FASTQ_RECOVERY_COMPLETE.txt" ]]; then
    test ! -e "$ROOT/fastq/$SAMPLE.R1.fastq.gz"
    test ! -e "$NATIVE/fastq/$SAMPLE.R1.fastq.gz"
    echo "$(date -Is) Recovering paired FASTQs: $SAMPLE"
    samtools collate -O -u -@ 2 "$BAM" "$NATIVE/tmp/$SAMPLE.collate" |
        samtools fastq -@ 2 -c 1 -n -F 0x900 \
            -1 "$NATIVE/fastq/$SAMPLE.R1.fastq.gz" -2 "$NATIVE/fastq/$SAMPLE.R2.fastq.gz" \
            -s "$ROOT/fastq/$SAMPLE.singletons.fastq.gz" -0 "$ROOT/fastq/$SAMPLE.other.fastq.gz" - \
            2> "$ROOT/logs/$SAMPLE.fastq.log"
    test "$(gzip -cd "$ROOT/fastq/$SAMPLE.singletons.fastq.gz" | wc -c)" -eq 0
    test "$(gzip -cd "$ROOT/fastq/$SAMPLE.other.fastq.gz" | wc -c)" -eq 0
    python "$ROOT/pilot_qc.py" fastq "$NATIVE/fastq/$SAMPLE.R1.fastq.gz" "$NATIVE/fastq/$SAMPLE.R2.fastq.gz" \
        "$ROOT/qc/$SAMPLE.original_read_qc.json" "$ROOT/qc/$SAMPLE.fastq_validation.json"
    cp -- "$NATIVE/fastq/$SAMPLE.R1.fastq.gz" "$NATIVE/fastq/$SAMPLE.R2.fastq.gz" "$ROOT/fastq/"
    date -Is > "$ROOT/qc/$SAMPLE.FASTQ_RECOVERY_COMPLETE.txt"
fi
if [[ ! -s "$ROOT/qc/$SAMPLE.ALIGNMENT_FILE_COMPLETE.txt" ]]; then
    test ! -e "$ROOT/alignment/$SAMPLE.YPARmasked.bam"
    test ! -e "$NATIVE/alignment/$SAMPLE.YPARmasked.bam"
    test -s "$NATIVE/fastq/$SAMPLE.R1.fastq.gz"
    echo "$(date -Is) Alignment started: $SAMPLE"
    # The fixed -K matches A549's 160-million-base batch size. Five mapping
    # workers per sample permit two jobs within the available 12 CPUs.
    RG="@RG\tID:$SAMPLE\tSM:$SAMPLE\tPL:ILLUMINA"
    bwa-mem2 mem -t 5 -K 160000000 -R "$RG" \
        "$NATIVE/reference/chm13v2.0.YPARmasked.fa" \
        "$NATIVE/fastq/$SAMPLE.R1.fastq.gz" "$NATIVE/fastq/$SAMPLE.R2.fastq.gz" \
        2> "$ROOT/logs/$SAMPLE.alignment.log" |
        samtools sort -@ 1 -m 1G -T "$NATIVE/tmp/$SAMPLE.sort" -o "$NATIVE/alignment/$SAMPLE.YPARmasked.bam" -
    samtools quickcheck -v "$NATIVE/alignment/$SAMPLE.YPARmasked.bam"
    samtools index -@ 2 "$NATIVE/alignment/$SAMPLE.YPARmasked.bam"
    cp -- "$NATIVE/alignment/$SAMPLE.YPARmasked.bam" "$NATIVE/alignment/$SAMPLE.YPARmasked.bam.bai" "$ROOT/alignment/"
    samtools flagstat -@ 2 "$ROOT/alignment/$SAMPLE.YPARmasked.bam" > "$ROOT/qc/$SAMPLE.masked_flagstat.txt"
    samtools view -H "$ROOT/alignment/$SAMPLE.YPARmasked.bam" > "$ROOT/qc/$SAMPLE.masked_header.sam"
    date -Is > "$ROOT/qc/$SAMPLE.ALIGNMENT_FILE_COMPLETE.txt"
fi
echo "$(date -Is) Validating realigned read content: $SAMPLE"
python "$ROOT/pilot_qc.py" bam "$ROOT/alignment/$SAMPLE.YPARmasked.bam" "$ROOT/qc/$SAMPLE.masked_read_qc.json"
python - "$ROOT" "$SAMPLE" <<'PY'
import json,sys
from pathlib import Path
p=Path(sys.argv[1])/'qc';s=sys.argv[2]
a=json.loads((p/f'{s}.original_read_qc.json').read_text())
b=json.loads((p/f'{s}.masked_read_qc.json').read_text())
assert a['read_fingerprint']==b['read_fingerprint'],'Read content changed'
assert [(r['SN'],r['LN']) for r in a['header']['SQ']]==[(r['SN'],r['LN']) for r in b['header']['SQ']]
PY
sha256sum -c "$ROOT/qc/$SAMPLE.original.sha256"
date -Is > "$ROOT/qc/$SAMPLE.REALIGNMENT_COMPLETE.txt"
echo "$(date -Is) COMPLETE: $SAMPLE; original unchanged and all primary read content preserved"
