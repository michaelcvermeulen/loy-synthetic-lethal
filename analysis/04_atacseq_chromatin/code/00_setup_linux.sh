#!/usr/bin/env bash
set -euo pipefail
run_root="$(cd "$(dirname "$0")/.." && pwd)"
atac_env=/home/micha/.local/share/atacseq-sept08-2026-venv
bootstrap=/tmp/atacseq-sept08-2026-virtualenv.pyz
if ! command -v x86_64-linux-gnu-gcc >/dev/null || [ ! -f /usr/include/python3.12/Python.h ]; then
    sudo -n apt-get update
    sudo -n apt-get install -y --no-install-recommends build-essential python3-dev
fi
if [ ! -x "$atac_env/bin/python" ]; then
    python3 - "$bootstrap" <<'PY'
import sys, urllib.request
urllib.request.urlretrieve('https://bootstrap.pypa.io/virtualenv.pyz', sys.argv[1])
PY
    python3 "$bootstrap" --no-periodic-update "$atac_env"
fi
requirements="$run_root/config/requirements.in"
if [ -f "$run_root/config/requirements.lock.txt" ]; then requirements="$run_root/config/requirements.lock.txt"; fi
"$atac_env/bin/python" -m pip install -r "$requirements"
"$atac_env/bin/python" -m pip check
"$atac_env/bin/python" -m pip freeze > "$run_root/config/requirements.lock.txt"
"$atac_env/bin/macs3" --version
"$atac_env/bin/macs3" callpeak --help > "$run_root/logs/macs3_callpeak_help.txt"
"$atac_env/bin/python" - <<'PY'
import pysam
print('pysam', pysam.__version__, 'samtools', pysam.__samtools_version__)
PY
