#!/usr/bin/env bash
set -euo pipefail
printf 'Operating system\n'
cat /etc/os-release
printf '\nTool paths\n'
for tool in samtools macs3 macs2 bedtools python3 Rscript conda micromamba pip3; do
    command -v "$tool" || true
done
printf '\nResources\n'
free -h
df -h / /mnt/a /mnt/h 2>/dev/null || true
printf '\nProject access\n'
ls -ld /mnt/a/Dropbox/LOY_paralogs/projects/ATACseq/filtered 2>/dev/null || true
printf '\nPython modules\n'
python3 - <<'PY'
import importlib.util
import sys
print(sys.version)
for name in ['pysam', 'MACS3', 'numpy', 'scipy', 'pandas', 'pyBigWig', 'matplotlib']:
    spec = importlib.util.find_spec(name)
    print(name, spec.origin if spec else 'not installed')
PY
printf '\nEnvironment managers\n'
for path in /home/micha/miniconda3/bin/conda /home/micha/miniforge3/bin/conda /opt/conda/bin/conda /home/micha/.local/bin/uv; do
    if [ -x "$path" ]; then printf '%s\n' "$path"; fi
done
python3 -c 'import importlib.util; print("ensurepip:", importlib.util.find_spec("ensurepip"))'
