#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
mkdir -p "$ROOT/runtime" "$ROOT/logs" "$ROOT/reference" "$ROOT/qc" "$ROOT/fastq" "$ROOT/alignment" "$ROOT/tmp"
cd "$ROOT/runtime"
if [[ ! -x bin/micromamba ]]; then
    curl --fail --location --retry 3 https://micro.mamba.pm/api/micromamba/linux-64/latest -o micromamba.tar.bz2
    tar -xjf micromamba.tar.bz2 bin/micromamba
fi
export MAMBA_ROOT_PREFIX="$ROOT/runtime/mamba"
if [[ ! -x "$ROOT/runtime/env/bin/bwa-mem2" ]]; then
    "$ROOT/runtime/bin/micromamba" create -y -p "$ROOT/runtime/env" -c conda-forge -c bioconda --strict-channel-priority \
        bwa-mem2=2.2.1 samtools=1.23.1 pysam=0.23.3 python=3.12
fi
"$ROOT/runtime/bin/micromamba" list -p "$ROOT/runtime/env" --explicit > "$ROOT/runtime/environment_explicit.txt"
"$ROOT/runtime/env/bin/samtools" --version > "$ROOT/runtime/samtools_version.txt"
"$ROOT/runtime/env/bin/bwa-mem2" version > "$ROOT/runtime/bwa_mem2_version.txt" 2>&1
curl --fail --location --retry 3 \
    https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0_PAR.bed \
    -o "$ROOT/reference/chm13v2.0_PAR.bed"
echo "Runtime ready"
