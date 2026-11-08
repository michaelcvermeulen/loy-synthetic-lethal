SAMTOOLS=/global/home/hpc5173/miniconda3/envs/sv_analysis/bin/samtools
MOSDEPTH=/global/home/hpc5173/miniconda3/envs/sv_analysis/bin/mosdepth
BAMDIR=/global/home/hpc5173/projects/LOY/WGS/alignment/t2t_bam
CHRYBAMDIR=/global/home/hpc5173/projects/LOY/WGS/alignment/t2t_chrY
COVDIR=/global/home/hpc5173/projects/LOY/WGS/coverage

SAMPLES="A549-CTRL2 A549-LOY A549-BC14"

for i in $SAMPLES; do

echo '#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=64GB
#SBATCH --cpus-per-task=16
#SBATCH --job-name=cov_'"${i}"'
#SBATCH --time=2-0
' > tmp.sh; echo "

# Extract chrY reads
${SAMTOOLS} view -b -@ 16 -o ${CHRYBAMDIR}/${i}.chrY.bam ${BAMDIR}/${i}.bam chrY
${SAMTOOLS} index ${CHRYBAMDIR}/${i}.chrY.bam

# chrY coverage in 10kb windows
${MOSDEPTH} --by 10000 --chrom chrY -t 16 ${COVDIR}/${i}.chrY ${BAMDIR}/${i}.bam

# Genome-wide summary (for normalization)
${MOSDEPTH} --no-per-base -t 16 ${COVDIR}/${i}.genome ${BAMDIR}/${i}.bam

" >> tmp.sh; sbatch tmp.sh; done
