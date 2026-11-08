#!/usr/bin/env bash
# Run under WSL root solely to add/remove a task-specific temporary swap file.
# BWA itself runs as the normal WSL user. No WSL configuration is edited.
set -euo pipefail
ROOT=/mnt/a/Dropbox/LOY_paralogs/projects/WGS/reanalysis/Y_PAR_masked_20261027
NATIVE=/home/micha/codex-work/wgs-ypar-pilot-20261027
SWAP="$NATIVE/index.swap"
test -s "$ROOT/reference/official_reference_verification.json"
mkdir -p "$NATIVE/reference" "$NATIVE/tmp"
chown micha:micha "$NATIVE" "$NATIVE/reference" "$NATIVE/tmp"
test ! -e "$SWAP"
cleanup() {
    if grep -qF "$SWAP" /proc/swaps; then swapoff "$SWAP"; fi
    if [[ -f "$SWAP" ]] && ! grep -qF "$SWAP" /proc/swaps; then rm -- "$SWAP"; fi
}
trap cleanup EXIT
fallocate -l 16G "$SWAP"
chmod 600 "$SWAP"
mkswap "$SWAP"
swapon "$SWAP"
cp -- "$ROOT/reference/chm13v2.0.YPARmasked.fa" "$NATIVE/reference/chm13v2.0.YPARmasked.fa"
chown micha:micha "$NATIVE/reference/chm13v2.0.YPARmasked.fa"
runuser -u micha -- "$ROOT/runtime/env/bin/samtools" faidx "$NATIVE/reference/chm13v2.0.YPARmasked.fa"
runuser -u micha -- /usr/bin/time -v "$ROOT/runtime/env/bin/bwa-mem2" index "$NATIVE/reference/chm13v2.0.YPARmasked.fa"
date -Is > "$ROOT/reference/INDEX_COMPLETE.txt"
