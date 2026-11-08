#!/bin/bash
set -euo pipefail
R=/global/home/hpc5173/projects/LOY/RNA
B=$R/updated/PRJNA930361
RS=/global/home/hpc5173/miniconda3/envs/rnaseq_analysis/bin/Rscript
echo "=================== merge_counts ==================="
$RS $R/updated/code/merge_counts.R $B
echo; echo "=================== qc_report ==================="
bash $R/updated/code/qc_report.sh $B || echo "(qc_report returned non-zero)"
echo; echo "=================== chrY / PAR check ==================="
$RS $B/check_chrY_par.R $B
