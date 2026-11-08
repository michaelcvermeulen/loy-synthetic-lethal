#!/bin/bash
#SBATCH --job-name=multiqc_ht1080_screen2
#SBATCH --output=/global/home/hpc5173/projects/CRISPR_screen/analysis/ht1080_screen_2/logs/multiqc_%j.out
#SBATCH --error=/global/home/hpc5173/projects/CRISPR_screen/analysis/ht1080_screen_2/logs/multiqc_%j.err
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G

# HT1080 Screen 2 MultiQC Aggregation Script
# Purpose: Aggregate FastQC results into a single report

# Set up environment
source /global/home/hpc5173/miniconda3/etc/profile.d/conda.sh
conda activate crispr-screen

# Define directories
FASTQC_DIR="/global/home/hpc5173/projects/CRISPR_screen/analysis/ht1080_screen_2/fastqc"
OUTPUT_DIR="/global/home/hpc5173/projects/CRISPR_screen/analysis/ht1080_screen_2/multiqc"

# Create output directory if it doesn't exist
mkdir -p ${OUTPUT_DIR}

echo "==========================================="
echo "Starting MultiQC for HT1080 Screen 2"
echo "Date: $(date)"
echo "FastQC directory: ${FASTQC_DIR}"
echo "Output directory: ${OUTPUT_DIR}"
echo "==========================================="

# Check if FastQC results exist
num_fastqc=$(ls -1 ${FASTQC_DIR}/*_fastqc.zip 2>/dev/null | wc -l)
if [ ${num_fastqc} -eq 0 ]; then
    echo "ERROR: No FastQC results found in ${FASTQC_DIR}"
    echo "Please run FastQC first: sbatch code/ht1080_screen2_01_fastqc.sh"
    exit 1
fi

echo "Found ${num_fastqc} FastQC results to aggregate"

# Run MultiQC
multiqc \
    ${FASTQC_DIR} \
    -o ${OUTPUT_DIR} \
    -n ht1080_screen2_multiqc_report \
    --title "HT1080 Screen 2 Quality Control Report" \
    --comment "Brunello CRISPR Screen - HT1080 Cell Line (Screen 2)" \
    --force \
    --interactive

# Check if MultiQC completed successfully
if [ $? -eq 0 ]; then
    echo "MultiQC analysis completed successfully!"
    echo "Report saved to: ${OUTPUT_DIR}/ht1080_screen2_multiqc_report.html"

    # Create summary of key metrics
    echo ""
    echo "Creating summary of key QC metrics..."
    echo "Sample Summary:" > ${OUTPUT_DIR}/qc_summary.txt
    echo "===============" >> ${OUTPUT_DIR}/qc_summary.txt
    echo "Total samples: ${num_fastqc}" >> ${OUTPUT_DIR}/qc_summary.txt
    echo "" >> ${OUTPUT_DIR}/qc_summary.txt

    # Extract sample names from FastQC reports
    echo "Samples analyzed:" >> ${OUTPUT_DIR}/qc_summary.txt
    for report in ${FASTQC_DIR}/*_fastqc.zip; do
        basename ${report} _fastqc.zip >> ${OUTPUT_DIR}/qc_summary.txt
    done

    echo ""
    echo "Summary saved to: ${OUTPUT_DIR}/qc_summary.txt"
else
    echo "ERROR: MultiQC analysis failed!"
    exit 1
fi

echo ""
echo "==========================================="
echo "MultiQC Analysis Complete"
echo "End time: $(date)"
echo "==========================================="
echo ""
echo "Next steps:"
echo "1. Review aggregated report: ${OUTPUT_DIR}/ht1080_screen2_multiqc_report.html"
echo "2. Proceed with trim optimization: sbatch code/ht1080_screen2_03_trim_optimization.sh"