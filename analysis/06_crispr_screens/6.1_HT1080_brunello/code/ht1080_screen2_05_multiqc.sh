#!/bin/bash
#SBATCH --job-name=multiqc_ht1080
#SBATCH --output=/global/home/hpc5173/projects/CRISPR_screen/analysis/ht1080_screen_2/logs/multiqc_%j.out
#SBATCH --error=/global/home/hpc5173/projects/CRISPR_screen/analysis/ht1080_screen_2/logs/multiqc_%j.err
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G

# HT1080 Screen 2 MultiQC Script
# Purpose: Generate comprehensive quality control report

# Set up environment
eval "$(conda shell.bash hook)"
conda activate crispr-screen

# Define directories
ANALYSIS_DIR="/global/home/hpc5173/projects/CRISPR_screen/analysis/ht1080_screen_2"
FASTQC_DIR="${ANALYSIS_DIR}/fastqc"
MAGECK_DIR="${ANALYSIS_DIR}/mageck/count"
OUTPUT_DIR="${ANALYSIS_DIR}/multiqc"

# Create output directory
mkdir -p ${OUTPUT_DIR}

echo "==========================================="
echo "Starting MultiQC for HT1080 Screen 2"
echo "Date: $(date)"
echo "==========================================="

# Run MultiQC on all QC results
multiqc \
    ${FASTQC_DIR} \
    ${MAGECK_DIR} \
    -n ht1080_screen2_multiqc_report \
    -o ${OUTPUT_DIR} \
    --interactive \
    --export \
    2>&1 | tee ${OUTPUT_DIR}/multiqc.log

if [ $? -eq 0 ]; then
    echo "MultiQC completed successfully!"
    echo ""
    echo "Files generated:"
    ls -lh ${OUTPUT_DIR}/*

    echo ""
    echo "Key QC metrics summary:"
    echo "----------------------"
    # Extract key metrics from the report
    if [ -f "${OUTPUT_DIR}/multiqc_data/multiqc_general_stats.txt" ]; then
        echo "General statistics available in: ${OUTPUT_DIR}/multiqc_data/multiqc_general_stats.txt"
    fi
else
    echo "ERROR: MultiQC failed!"
    exit 1
fi

echo ""
echo "==========================================="
echo "MultiQC Complete!"
echo "End time: $(date)"
echo "==========================================="
echo ""
echo "Report location: ${OUTPUT_DIR}/ht1080_screen2_multiqc_report.html"
echo ""
echo "Next steps:"
echo "1. Review MultiQC report for quality metrics"
echo "2. Review MAGeCK count summary"
echo "3. Proceed with MAGeCK test and MLE if QC passes"