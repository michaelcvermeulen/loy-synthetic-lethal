cd /global/home/hpc5173/projects/LOY/WGS/fastp

for i in $(ls *.gz* | grep -oP '^[A-Za-z0-9\-]+(?=_R[12]_001)' | sort -u); do

echo '#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=64GB
#SBATCH --cpus-per-task=16
#SBATCH --job-name=test_${SAMPLE}
#SBATCH --time=2-0
' > tmp.sh; echo "

/global/home/hpc5173/miniconda3/envs/hisat2/bin/bwa-mem2 mem -t 16 \
/global/home/hpc5173/data/GRCh38.primary_assembly.genome.fa \
/global/home/hpc5173/projects/LOY/WGS/fastp/${i}_R1_001.fastq.gz \
/global/home/hpc5173/projects/LOY/WGS/fastp/${i}_R2_001.fastq.gz | /global/home/hpc5173/miniconda3/envs/main/bin/samtools sort -@ 16 -o \
/global/home/hpc5173/projects/LOY/WGS/alignment/bam/${i}.bam

" >> tmp.sh; sbatch tmp.sh; done
