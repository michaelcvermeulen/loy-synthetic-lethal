REF=/global/home/hpc5173/data/chm13v2.0.fa
MANTA_CONFIG=/global/home/hpc5173/miniconda3/envs/sv_analysis/bin/configManta.py
BAMDIR=/global/home/hpc5173/projects/LOY/WGS/alignment/t2t_bam
SVDIR=/global/home/hpc5173/projects/LOY/WGS/sv/manta

SAMPLES="A549-CTRL2 A549-LOY A549-BC14"

# --- Per-sample Manta calling ---
for i in $SAMPLES; do

echo '#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=64GB
#SBATCH --cpus-per-task=16
#SBATCH --job-name=manta_'"${i}"'
#SBATCH --time=2-0
' > tmp.sh; echo "

${MANTA_CONFIG} \
  --bam ${BAMDIR}/${i}.bam \
  --referenceFasta ${REF} \
  --runDir ${SVDIR}/${i}

${SVDIR}/${i}/runWorkflow.py -j 16

" >> tmp.sh; sbatch tmp.sh; done

# --- Paired analysis: CTRL2 (normal) vs LOY (tumor) ---
echo '#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=64GB
#SBATCH --cpus-per-task=16
#SBATCH --job-name=manta_paired
#SBATCH --time=2-0
' > tmp.sh; echo "

${MANTA_CONFIG} \
  --normalBam ${BAMDIR}/A549-CTRL2.bam \
  --tumorBam ${BAMDIR}/A549-LOY.bam \
  --referenceFasta ${REF} \
  --runDir ${SVDIR}/A549_paired

${SVDIR}/A549_paired/runWorkflow.py -j 16

" >> tmp.sh; sbatch tmp.sh
