#!/bin/bash
#SBATCH --job-name=mageck_count_ht1080
#SBATCH --output=/global/home/hpc5173/projects/CRISPR_screen/analysis/ht1080_screen_2/logs/mageck_count_%j.out
#SBATCH --error=/global/home/hpc5173/projects/CRISPR_screen/analysis/ht1080_screen_2/logs/mageck_count_%j.err
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

# HT1080 Screen 2 MAGeCK Count Script
# Purpose: Count sgRNA abundances in all samples using optimized parameters

# Set up environment
source /global/home/hpc5173/miniconda3/etc/profile.d/conda.sh
conda activate crispr-screen

# Define directories
DATA_DIR="/global/home/hpc5173/projects/CRISPR_screen/data/ht1080_screen_2_raw"
OUTPUT_DIR="/global/home/hpc5173/projects/CRISPR_screen/analysis/ht1080_screen_2/mageck/count"
LIBRARY="/global/home/hpc5173/projects/CRISPR_screen/data/brunello/broadgpp-brunello-library-corrected.txt"
CONTROL_SGRNA="/global/home/hpc5173/projects/CRISPR_screen/data/brunello/NTC_control_library.txt"

# Create output directory
mkdir -p ${OUTPUT_DIR}

# PARAMETERS - Based on barcode analysis
R1_TRIM=22   # sgRNA starts at position 22 in R1
# R2 not needed - R1 contains all sgRNA information

echo "==========================================="
echo "Starting MAGeCK Count for HT1080 Screen 2"
echo "Date: $(date)"
echo "==========================================="
echo "Parameters:"
echo "  R1 trim-5: ${R1_TRIM}"
echo "  R2 trim-5: ${R2_TRIM}"
echo "  Library: Brunello (77,441 sgRNAs)"
echo "  Normalization: Control sgRNA"
echo "==========================================="

# Process R1 reads
echo ""
echo "Processing R1 reads..."
echo "---------------------"

# Prepare file list for R1
R1_FILES="${DATA_DIR}/D0-CTRL_R1_001.fastq.gz"
R1_FILES="${R1_FILES},${DATA_DIR}/D0-LOY_R1_001.fastq.gz"
R1_FILES="${R1_FILES},${DATA_DIR}/CTRL-REP1_R1_001.fastq.gz"
R1_FILES="${R1_FILES},${DATA_DIR}/CTRL-REP2_R1_001.fastq.gz"
R1_FILES="${R1_FILES},${DATA_DIR}/LOY-REP1_R1_001.fastq.gz"
R1_FILES="${R1_FILES},${DATA_DIR}/LOY-REP2_R1_001.fastq.gz"

# Sample labels
LABELS="D0_CTRL,D0_LOY,POST_CTRL_REP1,POST_CTRL_REP2,POST_LOY_REP1,POST_LOY_REP2"

echo "Running MAGeCK count for R1 reads..."
mageck count \
    -l "${LIBRARY}" \
    --sample-label "${LABELS}" \
    --fastq "${R1_FILES}" \
    --day0-label "D0_CTRL,D0_LOY" \
    --pdf-report \
    --unmapped-to-file \
    -n ${OUTPUT_DIR}/ht1080_screen2_R1 \
    2>&1 | tee ${OUTPUT_DIR}/ht1080_screen2_R1_count.log

# Check if R1 processing was successful
if [ $? -eq 0 ]; then
    echo "R1 processing complete!"

    # Extract and display mapping statistics
    echo ""
    echo "R1 Mapping Statistics:"
    grep -A 10 "Summary of samples" ${OUTPUT_DIR}/ht1080_screen2_R1_count.log
else
    echo "ERROR: R1 processing failed!"
    exit 1
fi

# Note: R2 reads not processed - analysis showed sgRNA barcodes are only in R1
# R2 reads are reverse complement from opposite direction and don't contain usable sgRNA sequences

echo ""
echo "==========================================="
echo "Creating summary statistics..."
echo "==========================================="

# Create a summary file
SUMMARY_FILE="${OUTPUT_DIR}/count_summary.txt"
echo "MAGeCK Count Summary - HT1080 Screen 2" > ${SUMMARY_FILE}
echo "Date: $(date)" >> ${SUMMARY_FILE}
echo "" >> ${SUMMARY_FILE}
echo "Parameters used:" >> ${SUMMARY_FILE}
echo "  R1 trim-5: ${R1_TRIM}" >> ${SUMMARY_FILE}
echo "  R2 trim-5: ${R2_TRIM}" >> ${SUMMARY_FILE}
echo "" >> ${SUMMARY_FILE}

# Extract key statistics
echo "R1 Statistics:" >> ${SUMMARY_FILE}
echo "--------------" >> ${SUMMARY_FILE}
grep "Percentage:" ${OUTPUT_DIR}/ht1080_screen2_R1_count.log >> ${SUMMARY_FILE}
echo "" >> ${SUMMARY_FILE}

# R2 statistics not included as R2 reads were not processed

# Count detected sgRNAs
echo "sgRNA Detection:" >> ${SUMMARY_FILE}
echo "----------------" >> ${SUMMARY_FILE}
total_sgrnas=$(wc -l < ${LIBRARY})
r1_detected=$(awk 'NR>1 {sum=0; for(i=2;i<=NF;i++) if($i>0) {sum++; break}} END {print sum}' ${OUTPUT_DIR}/ht1080_screen2_R1.count.txt)
# r2_detected not calculated as R2 was not processed

echo "Total sgRNAs in library: ${total_sgrnas}" >> ${SUMMARY_FILE}
echo "sgRNAs detected (R1): ${r1_detected}" >> ${SUMMARY_FILE}
echo "Note: R2 reads not processed - sgRNA barcodes only in R1" >> ${SUMMARY_FILE}

echo ""
echo "Summary saved to: ${SUMMARY_FILE}"

echo ""
echo "==========================================="
echo "MAGeCK Count Complete!"
echo "End time: $(date)"
echo "==========================================="
echo ""
echo "Output files generated:"
echo "  Count matrix: ${OUTPUT_DIR}/ht1080_screen2_R1.count.txt"
echo "  Normalized counts: ${OUTPUT_DIR}/ht1080_screen2_R1.count_normalized.txt"
echo "  PDF report: ${OUTPUT_DIR}/ht1080_screen2_R1.countsummary.pdf"
echo "  Log file: ${OUTPUT_DIR}/ht1080_screen2_R1_count.log"
echo ""
echo "Next steps:"
echo "1. Review count QC metrics"
echo "2. Run MAGeCK test for differential analysis"
echo "3. Run MAGeCK-MLE for advanced analysis"