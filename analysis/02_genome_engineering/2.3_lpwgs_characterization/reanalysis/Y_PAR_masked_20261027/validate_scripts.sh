#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
for script in "$ROOT"/*.sh; do bash -n "$script"; done
"$ROOT/runtime/env/bin/python" -m py_compile "$ROOT"/*.py
echo 'Shell and Python syntax validated.'
free -h
ps -eo pid,pcpu,rss,etime,comm --sort=-rss | head -8
