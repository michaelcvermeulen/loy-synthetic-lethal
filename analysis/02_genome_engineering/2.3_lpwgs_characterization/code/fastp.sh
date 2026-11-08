cd /global/home/hpc5173/projects/LOY/WGS/00_fastq

for i in $(ls *.gz* | grep -oP '^[A-Za-z0-9\-]+(?=_R[12]_001)' | sort -u); do

echo '#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=64GB
#SBATCH --cpus-per-task=16
#SBATCH --job-name=test_${SAMPLE}
#SBATCH --time=2-0
' > tmp.sh; echo "

/global/home/hpc5173/miniconda3/envs/hisat2/bin/fastp \
 -i /global/home/hpc5173/projects/LOY/WGS/00_fastq/${i}_R1_001.fastq.gz \
 -I /global/home/hpc5173/projects/LOY/WGS/00_fastq/${i}_R2_001.fastq.gz \
 -o /global/home/hpc5173/projects/LOY/WGS/fastp/${i}_R1_001.fastq.gz \
 -O /global/home/hpc5173/projects/LOY/WGS/fastp/${i}_R2_001.fastq.gz \
  --trim_front1 15 \
  --trim_front2 15 \
  --length_required 50 \
  --detect_adapter_for_pe \
  --html /global/home/hpc5173/projects/LOY/WGS/00_fastq/${i}.html \
  --thread 16


" >> tmp.sh; sbatch tmp.sh; done

