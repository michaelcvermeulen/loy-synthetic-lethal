#!/bin/bash
# Copy the DepMap NextGen inputs of stages 31 to 33 from a folder of downloads
# into the collection and check each against the checksum recorded under
# "nextgen_inputs" in extension_parameters.json. The DepMap portal answers
# scripts with a verification page, so these files are downloaded by browser.
# A file already in place is verified, never overwritten; the source folder is
# left as it is. SHA-256 sums go to docs/SHA256SUMS_nextgen_inputs.
# Usage: stage_nextgen_inputs.sh <folder with the downloads>
set -euo pipefail
source_dir="${1:?usage: stage_nextgen_inputs.sh <folder with the downloads>}"
root=/global/home/hpc5173/projects/PRISM/data/CCLE/26Q1
python=/global/home/hpc5173/miniconda3/bin/python
cd "$root"
checksum_of() {
  if [ "$1" = "md5" ]; then md5sum "$2" | cut -d' ' -f1; else sha256sum "$2" | cut -d' ' -f1; fi
}
list=$(mktemp)
"$python" - <<'PY' > "$list"
import json
inputs = json.load(open("scripts/loy_reanalysis_v1/extension_parameters.json"))["nextgen_inputs"]
for row in zip(inputs["file"], inputs["origin"], inputs["algorithm"], inputs["checksum"]):
    print("\t".join(row))
PY
while IFS=$'\t' read -r file origin algorithm checksum; do
  mkdir -p "$(dirname "$file")"
  if [ ! -e "$file" ]; then
    partial="$file.partial"
    case "$origin" in
      *.zip:*) unzip -p "$source_dir/${origin%%:*}" "${origin#*:}" > "$partial" ;;
      *) cp "$source_dir/$origin" "$partial" ;;
    esac
    found=$(checksum_of "$algorithm" "$partial")
    if [ "$found" != "$checksum" ]; then
      rm -f "$partial"
      echo "checksum mismatch for $file: $found, expected $checksum" >&2
      exit 1
    fi
    mv "$partial" "$file"
    chmod a-w "$file"
    echo "copied     $file"
  fi
  found=$(checksum_of "$algorithm" "$file")
  if [ "$found" != "$checksum" ]; then
    echo "checksum mismatch for $file: $found, expected $checksum" >&2
    exit 1
  fi
  echo "verified   $file ($algorithm)"
done < "$list"
cut -f1 "$list" | xargs sha256sum > docs/SHA256SUMS_nextgen_inputs
rm -f "$list"
echo "wrote docs/SHA256SUMS_nextgen_inputs"
