#!/bin/bash
#SBATCH --job-name=mageck_count_ht1080
#SBATCH --output=/global/home/hpc5173/projects/CRISPR_screen/analysis/ht1080_screen_2/logs/mageck_count_fixed_%j.out
#SBATCH --error=/global/home/hpc5173/projects/CRISPR_screen/analysis/ht1080_screen_2/logs/mageck_count_fixed_%j.err
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

# HT1080 Screen 2 MAGeCK Count Script - Fixed Version
# Purpose: Count sgRNA abundances in all samples using optimized parameters

# Set up environment
eval "$(conda shell.bash hook)"
conda activate crispr-screen

# Define directories
DATA_DIR="/global/home/hpc5173/projects/CRISPR_screen/data/ht1080_screen_2_raw"
OUTPUT_DIR="/global/home/hpc5173/projects/CRISPR_screen/analysis/ht1080_screen_2/mageck/count"
LIBRARY="/global/home/hpc5173/projects/CRISPR_screen/data/brunello/broadgpp-brunello-library-corrected.txt"
CONTROL_SGRNA="/global/home/hpc5173/projects/CRISPR_screen/data/brunello/NTC_control_library.txt"

# Create output directory
mkdir -p ${OUTPUT_DIR}

echo "==========================================="
echo "Starting MAGeCK Count for HT1080 Screen 2"
echo "Date: $(date)"
echo "==========================================="
echo "Processing all 6 R1 samples with trim-5=22"
echo "==========================================="

# Run MAGeCK count with space-separated files
mageck count \
    -l "${LIBRARY}" \
    --sample-label D0_CTRL,D0_LOY,POST_CTRL_REP1,POST_CTRL_REP2,POST_LOY_REP1,POST_LOY_REP2 \
    --fastq ${DATA_DIR}/D0-CTRL_R1_001.fastq.gz \
            ${DATA_DIR}/D0-LOY_R1_001.fastq.gz \
            ${DATA_DIR}/CTRL-REP1_R1_001.fastq.gz \
            ${DATA_DIR}/CTRL-REP2_R1_001.fastq.gz \
            ${DATA_DIR}/LOY-REP1_R1_001.fastq.gz \
            ${DATA_DIR}/LOY-REP2_R1_001.fastq.gz \
    --trim-5 22 \
    --norm-method control \
    --control-sgrna "${CONTROL_SGRNA}" \
    --pdf-report \
    -n ${OUTPUT_DIR}/ht1080_screen2_R1 \
    2>&1 | tee ${OUTPUT_DIR}/ht1080_screen2_R1_count.log

# Check if processing was successful
if [ $? -eq 0 ]; then
    echo "MAGeCK count completed successfully!"

    # Extract and display mapping statistics
    echo ""
    echo "Mapping Statistics:"
    echo "-------------------"
    grep -E "Total:|Mapped:|Percentage:" ${OUTPUT_DIR}/ht1080_screen2_R1_count.log

    echo ""
    echo "Files generated:"
    ls -lh ${OUTPUT_DIR}/ht1080_screen2_R1*
else
    echo "ERROR: MAGeCK count failed!"
    exit 1
fi

echo ""
echo "==========================================="
echo "MAGeCK Count Complete!"
echo "End time: $(date)"
echo "==========================================="
echo ""
echo "Next steps:"
echo "1. Review count QC metrics"
echo "2. Run MAGeCK test for differential analysis"
echo "3. Run MAGeCK-MLE for advanced analysis"