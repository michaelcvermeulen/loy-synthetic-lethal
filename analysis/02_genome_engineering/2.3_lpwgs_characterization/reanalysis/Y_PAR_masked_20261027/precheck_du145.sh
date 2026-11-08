#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
export PATH="$ROOT/runtime/env/bin:$PATH"
precheck() {
    local sample=$1
    local bam=/mnt/a/Dropbox/LOY_paralogs/projects/WGS/t2t_bam/$sample.bam
    if [[ -s "$ROOT/qc/$sample.ORIGINAL_QC_COMPLETE.txt" ]]; then return; fi
    sha256sum "$bam" > "$ROOT/qc/$sample.original.sha256"
    samtools quickcheck -v "$bam"
    samtools view -H "$bam" > "$ROOT/qc/$sample.original_header.sam"
    samtools flagstat "$bam" > "$ROOT/qc/$sample.original_flagstat.txt"
    python "$ROOT/pilot_qc.py" bam "$bam" "$ROOT/qc/$sample.original_read_qc.json"
    date -Is > "$ROOT/qc/$sample.ORIGINAL_QC_COMPLETE.txt"
    echo "$(date -Is) Original QC ready: $sample"
}
precheck LNCaP-sgCTRL &
control_pid=$!
precheck LNCaP-LOY-2X &
loy_pid=$!
wait "$control_pid"
wait "$loy_pid"
