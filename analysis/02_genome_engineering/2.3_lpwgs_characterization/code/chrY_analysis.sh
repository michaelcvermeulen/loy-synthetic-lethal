BCFTOOLS=/global/home/hpc5173/miniconda3/envs/sv_analysis/bin/bcftools
MERGEDIR=/global/home/hpc5173/projects/LOY/WGS/sv/merged
CHRYDIR=/global/home/hpc5173/projects/LOY/WGS/sv/chrY

echo '#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=64GB
#SBATCH --cpus-per-task=16
#SBATCH --job-name=chrY_analysis
#SBATCH --time=2-0
' > tmp.sh; echo "

# Extract chrY SVs from per-sample merged callsets
for i in A549-CTRL2 A549-LOY A549-BC14; do
  ${BCFTOOLS} view -r chrY ${MERGEDIR}/\${i}.merged.vcf > ${CHRYDIR}/\${i}.chrY.vcf
done

# Extract chrY SVs from Manta paired somatic calls
${BCFTOOLS} view -r chrY ${MERGEDIR}/A549_paired.somatic.vcf > ${CHRYDIR}/A549_paired.somatic.chrY.vcf

# Find SVs private to CTRL2 (Y+) vs LOY (Y-)
# Compress and index for bcftools isec
for i in A549-CTRL2 A549-LOY; do
  ${BCFTOOLS} sort ${CHRYDIR}/\${i}.chrY.vcf -O z -o ${CHRYDIR}/\${i}.chrY.vcf.gz
  ${BCFTOOLS} index ${CHRYDIR}/\${i}.chrY.vcf.gz
done

# isec: find private and shared SVs
${BCFTOOLS} isec -p ${CHRYDIR}/isec_CTRL2_vs_LOY \
  ${CHRYDIR}/A549-CTRL2.chrY.vcf.gz \
  ${CHRYDIR}/A549-LOY.chrY.vcf.gz

# Output:
#   isec_CTRL2_vs_LOY/0000.vcf = private to A549-CTRL2 (Y+)
#   isec_CTRL2_vs_LOY/0001.vcf = private to A549-LOY (Y-)
#   isec_CTRL2_vs_LOY/0002.vcf = shared (in CTRL2)
#   isec_CTRL2_vs_LOY/0003.vcf = shared (in LOY)

echo '--- chrY SV counts ---'
echo 'CTRL2-private:' && grep -c -v '^#' ${CHRYDIR}/isec_CTRL2_vs_LOY/0000.vcf || true
echo 'LOY-private:' && grep -c -v '^#' ${CHRYDIR}/isec_CTRL2_vs_LOY/0001.vcf || true
echo 'Shared:' && grep -c -v '^#' ${CHRYDIR}/isec_CTRL2_vs_LOY/0002.vcf || true

" >> tmp.sh; sbatch tmp.sh
