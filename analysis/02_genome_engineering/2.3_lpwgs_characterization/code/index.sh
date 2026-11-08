echo '#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=128GB
#SBATCH --cpus-per-task=16
#SBATCH --job-name=test_${SAMPLE}
#SBATCH --time=2-0
' > tmp.sh; echo "

cd /global/home/hpc5173/data
/global/home/hpc5173/miniconda3/envs/hisat2/bin/bwa-mem2 index /global/home/hpc5173/data/GRCh38.primary_assembly.genome.fa


" >> tmp.sh; sbatch tmp.sh



