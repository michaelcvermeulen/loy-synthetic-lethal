#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
FAILED=0
for CELL in HT1080 DU145; do
    if [[ "$CELL" == HT1080 ]]; then S1=HT1080-sgCTRL-D0;S2=HT1080-LOY-D0;else S1=LNCaP-sgCTRL;S2=LNCaP-LOY-2X;fi
    echo "$(date -Is) Starting pair: $CELL"
    bash "$ROOT/realign_atlas_sample.sh" "$S1" > "$ROOT/logs/$S1.pipeline.log" 2>&1 & P1=$!
    bash "$ROOT/realign_atlas_sample.sh" "$S2" > "$ROOT/logs/$S2.pipeline.log" 2>&1 & P2=$!
    if wait "$P1";then echo "Validated: $S1";else echo "FAILED: $S1; see its pipeline log";FAILED=1;fi
    if wait "$P2";then echo "Validated: $S2";else echo "FAILED: $S2; see its pipeline log";FAILED=1;fi
done
exit "$FAILED"
