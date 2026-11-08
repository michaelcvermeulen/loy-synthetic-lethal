cd /global/home/hpc5173/projects/LOY/WGS/t2t_bam

for i in A549-CTRL2 A549-LOY HT1080-sgCTRL-D0 HT1080-LOY-D0; do

echo '#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=16GB
#SBATCH --cpus-per-task=4
#SBATCH --job-name=mosdepth_${SAMPLE}
#SBATCH --time=1-0
' > tmp.sh; echo "

/global/home/hpc5173/miniconda3/envs/main/bin/mosdepth \
  --by 100000 \
  --no-per-base \
  -t 4 \
  /global/home/hpc5173/projects/LOY/WGS/coverage/${i}.genome.100kb \
  /global/home/hpc5173/projects/LOY/WGS/t2t_bam/${i}.bam

" >> tmp.sh; sbatch tmp.sh; done
