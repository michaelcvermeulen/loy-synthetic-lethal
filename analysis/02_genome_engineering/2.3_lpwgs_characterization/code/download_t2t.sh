#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=128GB
#SBATCH --cpus-per-task=16
#SBATCH --job-name=download_t2t
#SBATCH --time=2-0

# Download T2T-CHM13v2.0 reference (chr-style names, includes chrY)
cd /global/home/hpc5173/data

wget -O chm13v2.0.fa.gz \
  https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz

gunzip chm13v2.0.fa.gz

# Build samtools faidx index
/global/home/hpc5173/miniconda3/envs/sv_analysis/bin/samtools faidx chm13v2.0.fa

# Build bwa-mem2 index
/global/home/hpc5173/miniconda3/envs/sv_analysis/bin/bwa-mem2 index chm13v2.0.fa
