REF=/global/home/hpc5173/data/chm13v2.0.fa
DELLY=/global/home/hpc5173/miniconda3/envs/sv_analysis/bin/delly
BCFTOOLS=/global/home/hpc5173/miniconda3/envs/sv_analysis/bin/bcftools
BAMDIR=/global/home/hpc5173/projects/LOY/WGS/alignment/t2t_bam
SVDIR=/global/home/hpc5173/projects/LOY/WGS/sv/delly

SAMPLES="A549-CTRL2 A549-LOY A549-BC14"

# --- Step 1: Per-sample SV discovery ---
for i in $SAMPLES; do

echo '#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=64GB
#SBATCH --cpus-per-task=16
#SBATCH --job-name=delly_'"${i}"'
#SBATCH --time=2-0
' > tmp.sh; echo "

${DELLY} call -g ${REF} -o ${SVDIR}/${i}.bcf ${BAMDIR}/${i}.bam

" >> tmp.sh; sbatch tmp.sh; done

echo ""
echo "=== WAIT for all per-sample jobs to finish before running the merge/genotype steps below ==="
echo "=== Then run: bash code/delly_merge_genotype.sh ==="
