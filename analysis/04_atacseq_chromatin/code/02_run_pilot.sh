#!/usr/bin/env bash
set -euo pipefail
run_root="$(cd "$(dirname "$0")/.." && pwd)"
atac_python=/home/micha/.local/share/atacseq-sept08-2026-venv/bin/python
"$atac_python" "$run_root/code/test_preparation.py" > "$run_root/logs/preparation_checks.log" 2>&1
"$atac_python" -u "$run_root/code/02_prepare.py" CTRL_11 > "$run_root/logs/CTRL_11_preparation.log" 2>&1
"$atac_python" -u "$run_root/code/03_call_peaks.py" CTRL_11
