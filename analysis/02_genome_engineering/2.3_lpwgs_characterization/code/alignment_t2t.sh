REF=/global/home/hpc5173/data/chm13v2.0.fa
SAMTOOLS=/global/home/hpc5173/miniconda3/envs/sv_analysis/bin/samtools
BWA=/global/home/hpc5173/miniconda3/envs/sv_analysis/bin/bwa-mem2
INDIR=/global/home/hpc5173/projects/LOY/WGS/fastp
OUTDIR=/global/home/hpc5173/projects/LOY/WGS/alignment/t2t_bam

# A549 samples only for now — add more samples to this list to expand
SAMPLES="A549-CTRL2 A549-LOY A549-BC14"

for i in $SAMPLES; do

echo '#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=64GB
#SBATCH --cpus-per-task=16
#SBATCH --job-name=t2t_'"${i}"'
#SBATCH --time=2-0
' > tmp.sh; echo "

${BWA} mem -t 16 \
  -R '@RG\tID:${i}\tSM:${i}\tPL:ILLUMINA' \
  ${REF} \
  ${INDIR}/${i}_R1_001.fastq.gz \
  ${INDIR}/${i}_R2_001.fastq.gz | \
  ${SAMTOOLS} sort -@ 16 -o ${OUTDIR}/${i}.bam

${SAMTOOLS} index ${OUTDIR}/${i}.bam

" >> tmp.sh; sbatch tmp.sh; done
