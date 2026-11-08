#!/bin/bash
# Run one stage from an immutable code snapshot so later edits cannot alter a running job.
# Usage: submit_stage.sh <stage_name> <script> [args...]
# LR_TIME, LR_MEMORY and LR_CPUS override resources; LR_WAIT=1 blocks until the job ends;
# LR_ARRAY=1-N runs the script as an array (it reads SLURM_ARRAY_TASK_ID).
set -euo pipefail
cd /global/home/hpc5173/projects/PRISM/data/CCLE/26Q1
stage_name="$1"
program_name="$2"
shift 2
scripts=scripts/loy_reanalysis_v1
results=results/loy_reanalysis_v1
mkdir -p "$results/code_snapshots" "$results/logs"
snapshot_path=$(mktemp -d "$results/code_snapshots/${stage_name}_XXXXXXXX")
for source_path in "$scripts"/*; do
  if [ -f "$source_path" ]; then
    cp "$source_path" "$snapshot_path/"
  fi
done
find "$snapshot_path" -maxdepth 1 -type f -print0 | sort -z | xargs -0 sha256sum > "$snapshot_path.code.sha256"
# Forest stages need ranger, which only the main environment has.
runner=/global/home/hpc5173/miniconda3/envs/scRNA/bin/Rscript
if [[ "$program_name" == *_forest_*.R ]]; then
  runner=/global/home/hpc5173/miniconda3/envs/main/bin/Rscript
fi
wait_flag=""
if [ "${LR_WAIT:-0}" = "1" ]; then
  wait_flag="--wait"
fi
array_flag=""
if [ -n "${LR_ARRAY:-}" ]; then
  array_flag="--array=${LR_ARRAY}"
fi
status=0
job_id=$(sbatch --parsable $wait_flag $array_flag --job-name="lr_${stage_name}" --time="${LR_TIME:-02:00:00}" --mem="${LR_MEMORY:-16G}" --cpus-per-task="${LR_CPUS:-4}" "$scripts/job.sbatch" "$runner" --vanilla "$snapshot_path/$program_name" "$@") || status=$?
printf '%s\t%s\t%s\t%s\n' "$(date -u +%FT%TZ)" "$stage_name" "$job_id" "$snapshot_path" >> "$results/submission_ledger.tsv"
printf '%s\n' "$job_id"
exit "$status"
