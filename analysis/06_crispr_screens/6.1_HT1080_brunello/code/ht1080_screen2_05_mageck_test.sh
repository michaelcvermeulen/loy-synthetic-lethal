#!/bin/bash
#SBATCH --job-name=mageck_test_ht1080
#SBATCH --output=/global/home/hpc5173/projects/CRISPR_screen/analysis/ht1080_screen_2/logs/mageck_test_%j.out
#SBATCH --error=/global/home/hpc5173/projects/CRISPR_screen/analysis/ht1080_screen_2/logs/mageck_test_%j.err
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

# HT1080 Screen 2 MAGeCK Test Script
# Purpose: Perform differential analysis between conditions

# Set up environment
source /global/home/hpc5173/miniconda3/etc/profile.d/conda.sh
conda activate crispr-screen

# Define directories
COUNT_DIR="/global/home/hpc5173/projects/CRISPR_screen/analysis/ht1080_screen_2/mageck/count"
OUTPUT_DIR="/global/home/hpc5173/projects/CRISPR_screen/analysis/ht1080_screen_2/mageck/test"

# Create output directory
mkdir -p ${OUTPUT_DIR}

echo "==========================================="
echo "Starting MAGeCK Test for HT1080 Screen 2"
echo "Date: $(date)"
echo "==========================================="

# Function to run MAGeCK test
run_mageck_test() {
    local count_file=$1
    local control_samples=$2
    local treatment_samples=$3
    local output_prefix=$4
    local description=$5

    echo ""
    echo "Running: ${description}"
    echo "Control: ${control_samples}"
    echo "Treatment: ${treatment_samples}"

    mageck test \
        -k ${count_file} \
        -t ${treatment_samples} \
        -c ${control_samples} \
        --norm-method median \
        --remove-zero both \
        --gene-test-fdr-threshold 0.05 \
        --adjust-method fdr \
        --variance-from-all-samples \
        --pdf-report \
        -n ${OUTPUT_DIR}/${output_prefix} \
        2>&1 | tee ${OUTPUT_DIR}/${output_prefix}.log

    if [ $? -eq 0 ]; then
        echo "✓ ${description} complete"

        # Extract top hits
        echo ""
        echo "Top 10 enriched genes:"
        head -n 11 ${OUTPUT_DIR}/${output_prefix}.gene_summary.txt | column -t

        echo ""
        echo "Top 10 depleted genes:"
        tail -n +2 ${OUTPUT_DIR}/${output_prefix}.gene_summary.txt | sort -k5 -g | head -n 10 | column -t
    else
        echo "✗ ${description} failed"
    fi
}

# Check if count files exist
echo "Checking for count files..."
if [ ! -f "${COUNT_DIR}/ht1080_screen2_R1.count.txt" ]; then
    echo "ERROR: Count file not found: ${COUNT_DIR}/ht1080_screen2_R1.count.txt"
    echo "Please run MAGeCK count first"
    exit 1
fi

echo "==========================================="
echo "Comparison 1: D0 LOY vs CTRL (Baseline)"
echo "==========================================="
run_mageck_test \
    "${COUNT_DIR}/ht1080_screen2_R1.count.txt" \
    "D0_CTRL" \
    "D0_LOY" \
    "D0_LOY_vs_CTRL_R1" \
    "D0 baseline comparison (R1)"

# Also run for R2 if available
if [ -f "${COUNT_DIR}/ht1080_screen2_R2.count.txt" ]; then
    run_mageck_test \
        "${COUNT_DIR}/ht1080_screen2_R2.count.txt" \
        "D0_CTRL" \
        "D0_LOY" \
        "D0_LOY_vs_CTRL_R2" \
        "D0 baseline comparison (R2)"
fi

echo ""
echo "==========================================="
echo "Comparison 2: POST LOY vs CTRL"
echo "==========================================="
run_mageck_test \
    "${COUNT_DIR}/ht1080_screen2_R1.count.txt" \
    "POST_CTRL_REP1,POST_CTRL_REP2" \
    "POST_LOY_REP1,POST_LOY_REP2" \
    "POST_LOY_vs_CTRL_R1" \
    "Post-selection comparison (R1)"

if [ -f "${COUNT_DIR}/ht1080_screen2_R2.count.txt" ]; then
    run_mageck_test \
        "${COUNT_DIR}/ht1080_screen2_R2.count.txt" \
        "POST_CTRL_REP1,POST_CTRL_REP2" \
        "POST_LOY_REP1,POST_LOY_REP2" \
        "POST_LOY_vs_CTRL_R2" \
        "Post-selection comparison (R2)"
fi

echo ""
echo "==========================================="
echo "Comparison 3: Temporal (POST vs D0) - CTRL"
echo "==========================================="
run_mageck_test \
    "${COUNT_DIR}/ht1080_screen2_R1.count.txt" \
    "D0_CTRL" \
    "POST_CTRL_REP1,POST_CTRL_REP2" \
    "POST_vs_D0_CTRL_R1" \
    "Temporal CTRL comparison (R1)"

echo ""
echo "==========================================="
echo "Comparison 4: Temporal (POST vs D0) - LOY"
echo "==========================================="
run_mageck_test \
    "${COUNT_DIR}/ht1080_screen2_R1.count.txt" \
    "D0_LOY" \
    "POST_LOY_REP1,POST_LOY_REP2" \
    "POST_vs_D0_LOY_R1" \
    "Temporal LOY comparison (R1)"

echo ""
echo "==========================================="
echo "Creating summary report..."
echo "==========================================="

# Create summary file
SUMMARY_FILE="${OUTPUT_DIR}/differential_analysis_summary.txt"
echo "MAGeCK Test Summary - HT1080 Screen 2" > ${SUMMARY_FILE}
echo "Date: $(date)" >> ${SUMMARY_FILE}
echo "" >> ${SUMMARY_FILE}

# Function to extract hit counts
extract_hits() {
    local file=$1
    local fdr_threshold=0.05

    if [ -f "${file}" ]; then
        enriched=$(awk -v thresh=${fdr_threshold} 'NR>1 && $9<=thresh && $3>0' ${file} | wc -l)
        depleted=$(awk -v thresh=${fdr_threshold} 'NR>1 && $10<=thresh && $3<0' ${file} | wc -l)
        echo "  Enriched genes (FDR < ${fdr_threshold}): ${enriched}"
        echo "  Depleted genes (FDR < ${fdr_threshold}): ${depleted}"
    else
        echo "  File not found: ${file}"
    fi
}

echo "Significant Hits Summary (FDR < 0.05):" >> ${SUMMARY_FILE}
echo "======================================" >> ${SUMMARY_FILE}
echo "" >> ${SUMMARY_FILE}

echo "D0 LOY vs CTRL:" >> ${SUMMARY_FILE}
extract_hits "${OUTPUT_DIR}/D0_LOY_vs_CTRL_R1.gene_summary.txt" >> ${SUMMARY_FILE}
echo "" >> ${SUMMARY_FILE}

echo "POST LOY vs CTRL:" >> ${SUMMARY_FILE}
extract_hits "${OUTPUT_DIR}/POST_LOY_vs_CTRL_R1.gene_summary.txt" >> ${SUMMARY_FILE}
echo "" >> ${SUMMARY_FILE}

echo "POST vs D0 (CTRL):" >> ${SUMMARY_FILE}
extract_hits "${OUTPUT_DIR}/POST_vs_D0_CTRL_R1.gene_summary.txt" >> ${SUMMARY_FILE}
echo "" >> ${SUMMARY_FILE}

echo "POST vs D0 (LOY):" >> ${SUMMARY_FILE}
extract_hits "${OUTPUT_DIR}/POST_vs_D0_LOY_R1.gene_summary.txt" >> ${SUMMARY_FILE}

echo ""
echo "Summary saved to: ${SUMMARY_FILE}"

echo ""
echo "==========================================="
echo "MAGeCK Test Complete!"
echo "End time: $(date)"
echo "==========================================="
echo ""
echo "Output files generated:"
echo "  Gene summaries: ${OUTPUT_DIR}/*gene_summary.txt"
echo "  sgRNA summaries: ${OUTPUT_DIR}/*sgrna_summary.txt"
echo "  PDF reports: ${OUTPUT_DIR}/*.pdf"
echo "  Log files: ${OUTPUT_DIR}/*.log"
echo ""
echo "Next steps:"
echo "1. Review gene summary files for significant hits"
echo "2. Examine PDF reports for quality metrics"
echo "3. Run MAGeCK-MLE for more sophisticated modeling"
echo "4. Generate visualizations (volcano plots, ranked plots)"