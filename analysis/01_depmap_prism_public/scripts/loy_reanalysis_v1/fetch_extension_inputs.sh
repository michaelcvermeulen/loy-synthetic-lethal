#!/bin/bash
# Download the public inputs of the extension stages (25 to 30) and check each
# against the checksum recorded in extension_parameters.json. Run on a node with
# internet access, from anywhere. A file already in place is verified, never
# overwritten. SHA-256 sums of all inputs go to docs/SHA256SUMS_extension_inputs.
set -euo pipefail
root=/global/home/hpc5173/projects/PRISM/data/CCLE/26Q1
python=/global/home/hpc5173/miniconda3/bin/python
cd "$root"
checksum_of() {
  if [ "$1" = "md5" ]; then md5sum "$2" | cut -d' ' -f1; else sha256sum "$2" | cut -d' ' -f1; fi
}
"$python" - <<'EOF' > /tmp/extension_inputs_$$.tsv
import json
inputs = json.load(open("scripts/loy_reanalysis_v1/extension_parameters.json"))["inputs"]
for row in zip(inputs["file"], inputs["url"], inputs["algorithm"], inputs["checksum"]):
    print("\t".join(row))
EOF
while IFS=$'\t' read -r file url algorithm checksum; do
  mkdir -p "$(dirname "$file")"
  if [ ! -e "$file" ]; then
    partial="$file.partial"
    curl --fail --silent --show-error --location --max-time 1800 --output "$partial" "$url"
    found=$(checksum_of "$algorithm" "$partial")
    if [ "$found" != "$checksum" ]; then
      echo "checksum mismatch for $file: $found, expected $checksum" >&2
      exit 1
    fi
    mv "$partial" "$file"
    chmod a-w "$file"
    echo "downloaded $file"
  fi
  found=$(checksum_of "$algorithm" "$file")
  if [ "$found" != "$checksum" ]; then
    echo "checksum mismatch for $file: $found, expected $checksum" >&2
    exit 1
  fi
  echo "verified   $file ($algorithm)"
done < /tmp/extension_inputs_$$.tsv
cut -f1 /tmp/extension_inputs_$$.tsv | xargs sha256sum > docs/SHA256SUMS_extension_inputs
rm -f /tmp/extension_inputs_$$.tsv
echo "wrote docs/SHA256SUMS_extension_inputs"
