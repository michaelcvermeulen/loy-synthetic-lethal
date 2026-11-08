#!/bin/bash
#SBATCH --job-name=fastqc_ht1080_screen2
#SBATCH --output=/global/home/hpc5173/projects/CRISPR_screen/analysis/ht1080_screen_2/logs/fastqc_%j.out
#SBATCH --error=/global/home/hpc5173/projects/CRISPR_screen/analysis/ht1080_screen_2/logs/fastqc_%j.err
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G

# HT1080 Screen 2 FastQC Analysis Script
# Purpose: Run FastQC quality control on all 12 FASTQ files

# Set up environment
source /global/home/hpc5173/miniconda3/etc/profile.d/conda.sh
conda activate crispr-screen

# Define directories
DATA_DIR="/global/home/hpc5173/projects/CRISPR_screen/data/ht1080_screen_2_raw"
OUTPUT_DIR="/global/home/hpc5173/projects/CRISPR_screen/analysis/ht1080_screen_2/fastqc"

# Create output directory if it doesn't exist
mkdir -p ${OUTPUT_DIR}

echo "==========================================="
echo "Starting FastQC Analysis for HT1080 Screen 2"
echo "Date: $(date)"
echo "Data directory: ${DATA_DIR}"
echo "Output directory: ${OUTPUT_DIR}"
echo "==========================================="

# Run FastQC on all FASTQ files
echo "Running FastQC with ${SLURM_CPUS_PER_TASK} threads..."

fastqc \
    -t ${SLURM_CPUS_PER_TASK} \
    -o ${OUTPUT_DIR} \
    ${DATA_DIR}/*.fastq.gz

# Check if FastQC completed successfully
if [ $? -eq 0 ]; then
    echo "FastQC analysis completed successfully!"
    echo "Results saved to: ${OUTPUT_DIR}"

    # Count the number of reports generated
    num_reports=$(ls -1 ${OUTPUT_DIR}/*.html 2>/dev/null | wc -l)
    echo "Generated ${num_reports} FastQC reports"

    # List all generated reports
    echo ""
    echo "Generated reports:"
    ls -1 ${OUTPUT_DIR}/*.html
else
    echo "ERROR: FastQC analysis failed!"
    exit 1
fi

echo ""
echo "==========================================="
echo "FastQC Analysis Complete"
echo "End time: $(date)"
echo "==========================================="
echo ""
echo "Next steps:"
echo "1. Review individual FastQC reports in: ${OUTPUT_DIR}"
echo "2. Run MultiQC to aggregate results: sbatch code/ht1080_screen2_02_multiqc.sh"