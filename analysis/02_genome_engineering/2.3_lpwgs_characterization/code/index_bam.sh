cd /global/home/hpc5173/projects/LOY/WGS/alignment/bam

for i in $(ls); do

echo '#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=64GB
#SBATCH --cpus-per-task=16
#SBATCH --job-name=test_${SAMPLE}
#SBATCH --time=2-0
' > tmp.sh; echo "

/global/home/hpc5173/miniconda3/envs/main/bin/samtools index  /global/home/hpc5173/projects/LOY/WGS/alignment/bam/${i} 
" >> tmp.sh; sbatch tmp.sh; done

