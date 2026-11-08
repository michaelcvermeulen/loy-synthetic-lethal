#cd /global/home/hpc5173/data
#wget https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz
#gunzip /global/home/hpc5173/data/grch38_genome.tar.gz
#tar -xvf grch38_genome.tar


cd /global/home/hpc5173/projects/LOY/WGS/fastp

for i in $(ls *.gz* | grep -oP '^[A-Za-z0-9\-]+(?=_R[12]_001)' | sort -u); do

echo '#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=64GB
#SBATCH --cpus-per-task=16
#SBATCH --job-name=test_${SAMPLE}
#SBATCH --time=2-0
' > tmp.sh; echo "


/global/home/hpc5173/miniconda3/envs/hisat2/bin/hisat2 -p 16 -q --no-spliced-alignment --maxins 800 \
-x /global/home/hpc5173/data/grch38/genome \
-1 /global/home/hpc5173/projects/LOY/WGS/fastp/${i}_R1_001.fastq.gz \
-2 /global/home/hpc5173/projects/LOY/WGS/fastp/${i}_R2_001.fastq.gz \
-S /global/home/hpc5173/projects/LOY/WGS/alignment/${i}.sam


# convert sam to bam format
/global/home/hpc5173/miniconda3/envs/main/bin/samtools view -b -o \
/global/home/hpc5173/projects/LOY/WGS/alignment/bam/${i}.bam \
-@ 16 /global/home/hpc5173/projects/LOY/WGS/alignment/${i}.sam
" >> tmp.sh; sbatch tmp.sh; done
