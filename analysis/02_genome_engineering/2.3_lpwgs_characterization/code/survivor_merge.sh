BCFTOOLS=/global/home/hpc5173/miniconda3/envs/sv_analysis/bin/bcftools
SURVIVOR=/global/home/hpc5173/miniconda3/envs/sv_analysis/bin/SURVIVOR
DELLYDIR=/global/home/hpc5173/projects/LOY/WGS/sv/delly
MANTADIR=/global/home/hpc5173/projects/LOY/WGS/sv/manta
MERGEDIR=/global/home/hpc5173/projects/LOY/WGS/sv/merged

SAMPLES="A549-CTRL2 A549-LOY A549-BC14"

echo '#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=64GB
#SBATCH --cpus-per-task=16
#SBATCH --job-name=survivor_merge
#SBATCH --time=2-0
' > tmp.sh; echo "

for i in ${SAMPLES}; do

  # Convert Delly BCF to VCF
  ${BCFTOOLS} view ${DELLYDIR}/\${i}.geno.bcf > ${MERGEDIR}/\${i}.delly.vcf

  # Copy Manta VCF (diploidSV for per-sample calls)
  cp ${MANTADIR}/\${i}/results/variants/diploidSV.vcf.gz ${MERGEDIR}/\${i}.manta.vcf.gz
  gunzip -f ${MERGEDIR}/\${i}.manta.vcf.gz

  # Create file list for SURVIVOR
  echo ${MERGEDIR}/\${i}.delly.vcf > ${MERGEDIR}/\${i}.filelist.txt
  echo ${MERGEDIR}/\${i}.manta.vcf >> ${MERGEDIR}/\${i}.filelist.txt

  # Merge: 1kb max distance, require 2 callers, all SV types, min 50bp
  ${SURVIVOR} merge ${MERGEDIR}/\${i}.filelist.txt 1000 2 1 1 0 50 ${MERGEDIR}/\${i}.merged.vcf

done

# Also merge the Manta paired tumor/normal results
cp ${MANTADIR}/A549_paired/results/variants/somaticSV.vcf.gz ${MERGEDIR}/A549_paired.somatic.vcf.gz
gunzip -f ${MERGEDIR}/A549_paired.somatic.vcf.gz

" >> tmp.sh; sbatch tmp.sh
