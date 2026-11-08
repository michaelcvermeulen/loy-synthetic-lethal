REF=/global/home/hpc5173/data/chm13v2.0.fa
DELLY=/global/home/hpc5173/miniconda3/envs/sv_analysis/bin/delly
BCFTOOLS=/global/home/hpc5173/miniconda3/envs/sv_analysis/bin/bcftools
BAMDIR=/global/home/hpc5173/projects/LOY/WGS/alignment/t2t_bam
SVDIR=/global/home/hpc5173/projects/LOY/WGS/sv/delly

SAMPLES="A549-CTRL2 A549-LOY A549-BC14"

# --- Step 2: Merge SV sites across samples ---
echo '#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=64GB
#SBATCH --cpus-per-task=16
#SBATCH --job-name=delly_merge
#SBATCH --time=2-0
' > tmp.sh; echo "

${DELLY} merge -o ${SVDIR}/sites.bcf ${SVDIR}/A549-CTRL2.bcf ${SVDIR}/A549-LOY.bcf ${SVDIR}/A549-BC14.bcf

# --- Step 3: Joint genotyping at merged sites ---
for i in ${SAMPLES}; do
  ${DELLY} call -g ${REF} -v ${SVDIR}/sites.bcf -o ${SVDIR}/\${i}.geno.bcf ${BAMDIR}/\${i}.bam
done

# --- Step 4: Merge genotyped BCFs and filter ---
${BCFTOOLS} merge ${SVDIR}/A549-CTRL2.geno.bcf ${SVDIR}/A549-LOY.geno.bcf ${SVDIR}/A549-BC14.geno.bcf -O b -o ${SVDIR}/all_samples.bcf
${BCFTOOLS} index ${SVDIR}/all_samples.bcf

# Filter by SV type
for svtype in DEL DUP INV BND; do
  ${DELLY} filter -t \${svtype} -o ${SVDIR}/\${svtype}_filtered.bcf ${SVDIR}/all_samples.bcf
done

" >> tmp.sh; sbatch tmp.sh
