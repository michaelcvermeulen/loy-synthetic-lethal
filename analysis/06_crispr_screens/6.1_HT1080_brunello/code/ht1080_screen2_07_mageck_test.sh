#!/bin/bash
#SBATCH --job-name=mageck_test_ht1080
#SBATCH --output=/global/home/hpc5173/projects/CRISPR_screen/analysis/ht1080_screen_2/logs/mageck_test_%j.out
#SBATCH --error=/global/home/hpc5173/projects/CRISPR_screen/analysis/ht1080_screen_2/logs/mageck_test_%j.err
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

# HT1080 Screen 2 MAGeCK Test Script
# Purpose: Perform differential analysis between conditions

# Set up environment
eval "$(conda shell.bash hook)"
conda activate crispr-screen

# Define directories
COUNT_DIR="/global/home/hpc5173/projects/CRISPR_screen/analysis/ht1080_screen_2/mageck/count"
TEST_DIR="/global/home/hpc5173/projects/CRISPR_screen/analysis/ht1080_screen_2/mageck/test"

# Create output directory
mkdir -p ${TEST_DIR}

echo "==========================================="
echo "Starting MAGeCK Test for HT1080 Screen 2"
echo "Date: $(date)"
echo "==========================================="
echo ""

# Check if count matrix exists
if [ ! -f "${COUNT_DIR}/ht1080_screen2_R1.count.txt" ]; then
    echo "ERROR: Count matrix not found!"
    echo "Please run MAGeCK count first."
    exit 1
fi

echo "Running differential analysis comparisons..."
echo ""

# 1. LOY vs CTRL at POST timepoint (primary comparison)
echo "1. POST LOY vs CTRL comparison..."
mageck test \
    -k "${COUNT_DIR}/ht1080_screen2_R1.count.txt" \
    -t POST_LOY_REP1,POST_LOY_REP2 \
    -c POST_CTRL_REP1,POST_CTRL_REP2 \
    --norm-method control \
    --control-sgrna /global/home/hpc5173/projects/CRISPR_screen/data/brunello/NTC_control_library.txt \
    --pdf-report \
    --gene-lfc-method alphamedian \
    -n "${TEST_DIR}/POST_LOY_vs_CTRL" \
    2>&1 | tee "${TEST_DIR}/POST_LOY_vs_CTRL.log"

echo ""

# 2. POST_CTRL vs D0_CTRL (control time course)
echo "2. POST_CTRL vs D0_CTRL comparison..."
mageck test \
    -k "${COUNT_DIR}/ht1080_screen2_R1.count.txt" \
    -t POST_CTRL_REP1,POST_CTRL_REP2 \
    -c D0_CTRL \
    --norm-method control \
    --control-sgrna /global/home/hpc5173/projects/CRISPR_screen/data/brunello/NTC_control_library.txt \
    --pdf-report \
    --gene-lfc-method alphamedian \
    -n "${TEST_DIR}/POST_CTRL_vs_D0_CTRL" \
    2>&1 | tee "${TEST_DIR}/POST_CTRL_vs_D0_CTRL.log"

echo ""

# 3. POST_LOY vs D0_LOY (LOY time course)
echo "3. POST_LOY vs D0_LOY comparison..."
mageck test \
    -k "${COUNT_DIR}/ht1080_screen2_R1.count.txt" \
    -t POST_LOY_REP1,POST_LOY_REP2 \
    -c D0_LOY \
    --norm-method control \
    --control-sgrna /global/home/hpc5173/projects/CRISPR_screen/data/brunello/NTC_control_library.txt \
    --pdf-report \
    --gene-lfc-method alphamedian \
    -n "${TEST_DIR}/POST_LOY_vs_D0_LOY" \
    2>&1 | tee "${TEST_DIR}/POST_LOY_vs_D0_LOY.log"

echo ""

# 4. D0_LOY vs D0_CTRL (baseline comparison)
echo "4. D0_LOY vs D0_CTRL comparison..."
mageck test \
    -k "${COUNT_DIR}/ht1080_screen2_R1.count.txt" \
    -t D0_LOY \
    -c D0_CTRL \
    --norm-method control \
    --control-sgrna /global/home/hpc5173/projects/CRISPR_screen/data/brunello/NTC_control_library.txt \
    --pdf-report \
    --gene-lfc-method alphamedian \
    -n "${TEST_DIR}/D0_LOY_vs_D0_CTRL" \
    2>&1 | tee "${TEST_DIR}/D0_LOY_vs_D0_CTRL.log"

echo ""
echo "==========================================="
echo "MAGeCK Test Complete!"
echo "End time: $(date)"
echo "==========================================="
echo ""

# Summary of results
echo "Summary of significant hits:"
echo "----------------------------"

for comparison in POST_LOY_vs_CTRL POST_CTRL_vs_D0_CTRL POST_LOY_vs_D0_LOY D0_LOY_vs_D0_CTRL; do
    if [ -f "${TEST_DIR}/${comparison}.gene_summary.txt" ]; then
        echo ""
        echo "${comparison}:"
        # Count significant genes (FDR < 0.05)
        n_pos=$(awk '$3 < 0.05' "${TEST_DIR}/${comparison}.gene_summary.txt" | wc -l)
        n_neg=$(awk '$6 < 0.05' "${TEST_DIR}/${comparison}.gene_summary.txt" | wc -l)
        echo "  Positively selected genes (FDR < 0.05): ${n_pos}"
        echo "  Negatively selected genes (FDR < 0.05): ${n_neg}"

        # Show top 5 hits for each direction
        echo "  Top 5 positively selected:"
        awk 'NR>1 && $3 < 0.05' "${TEST_DIR}/${comparison}.gene_summary.txt" | \
            sort -k3,3g | head -5 | awk '{print "    " $1 " (FDR=" $3 ", LFC=" $15 ")"}'

        echo "  Top 5 negatively selected:"
        awk 'NR>1 && $6 < 0.05' "${TEST_DIR}/${comparison}.gene_summary.txt" | \
            sort -k6,6g | head -5 | awk '{print "    " $1 " (FDR=" $6 ", LFC=" $15 ")"}'
    fi
done

echo ""
echo "Output files generated in: ${TEST_DIR}"
echo ""
echo "Key output files:"
echo "  - *gene_summary.txt: Gene-level statistics"
echo "  - *sgrna_summary.txt: sgRNA-level statistics"
echo "  - *.pdf: QC and visualization plots"
echo ""
echo "Next steps:"
echo "1. Review gene_summary.txt files for significant hits"
echo "2. Examine PDF reports for QC metrics"
echo "3. Run MAGeCK-MLE for pathway analysis"
echo "4. Generate custom visualizations"