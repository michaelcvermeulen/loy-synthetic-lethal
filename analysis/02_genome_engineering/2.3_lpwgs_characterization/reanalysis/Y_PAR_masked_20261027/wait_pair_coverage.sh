#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
cell=${1:?cell required}
case "$cell" in
    HT1080) samples=(HT1080-sgCTRL-D0 HT1080-LOY-D0);;
    DU145) samples=(LNCaP-sgCTRL LNCaP-LOY-2X);;
    *) exit 2;;
esac
deadline=$((SECONDS+14400))
while [[ ! -s "$ROOT/qc/${samples[0]}.ALIGNMENT_FILE_COMPLETE.txt" || ! -s "$ROOT/qc/${samples[1]}.ALIGNMENT_FILE_COMPLETE.txt" ]]; do
    if (( SECONDS >= deadline )); then echo "Timed out: $cell BAM files" >&2;exit 1;fi
    sleep 5
done
"$ROOT/runtime/env/bin/python" "$ROOT/atlas_pair_coverage.py" "$cell" masked
