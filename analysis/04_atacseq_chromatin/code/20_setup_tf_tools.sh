#!/usr/bin/env bash
set -euo pipefail
root=/mnt/a/Dropbox/LOY_paralogs/projects/ATACseq/SEPT_08_2026_REPROCESS
tf_tools=/tmp/atacseq-sept08-2026-tf-tools
mkdir -p "$tf_tools/homer"
bootstrap="$root/tools/TF_extension/virtualenv.pyz"
curl --fail --location --retry 2 https://bootstrap.pypa.io/virtualenv.pyz -o "$bootstrap"
python3 "$bootstrap" --no-periodic-update "$tf_tools/tobias_venv"
if [ -f "$root/config/TF_extension/tobias_requirements.lock.txt" ]; then
    "$tf_tools/tobias_venv/bin/python" -m pip install -r "$root/config/TF_extension/tobias_requirements.lock.txt" > "$root/logs/20_tobias_linux_install.log" 2>&1
else
    "$tf_tools/tobias_venv/bin/python" -m pip install tobias > "$root/logs/20_tobias_linux_install.log" 2>&1
fi
"$tf_tools/tobias_venv/bin/python" -m pip freeze > "$root/config/TF_extension/tobias_requirements.lock.txt"
if [ ! -f "$root/tools/TF_extension/homer-5.1.tar.bz2" ]; then
    curl --fail --location --retry 2 https://conda.anaconda.org/bioconda/linux-64/homer-5.1-pl5262h9948957_0.tar.bz2 -o "$root/tools/TF_extension/homer-5.1.tar.bz2"
fi
"$tf_tools/tobias_venv/bin/python" "$root/code/20b_unpack_homer.py" > "$root/logs/20_homer_install.log" 2>&1
