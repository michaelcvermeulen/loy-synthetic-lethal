#!/bin/bash
#SBATCH --job-name=mageck_mle_ht1080
#SBATCH --output=/global/home/hpc5173/projects/CRISPR_screen/analysis/ht1080_screen_2/logs/mageck_mle_%j.out
#SBATCH --error=/global/home/hpc5173/projects/CRISPR_screen/analysis/ht1080_screen_2/logs/mageck_mle_%j.err
#SBATCH --time=3:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

# HT1080 Screen 2 MAGeCK-MLE Script
# Purpose: Maximum Likelihood Estimation for more sophisticated modeling

# Set up environment
source /global/home/hpc5173/miniconda3/etc/profile.d/conda.sh
conda activate crispr-screen

# Define directories
COUNT_DIR="/global/home/hpc5173/projects/CRISPR_screen/analysis/ht1080_screen_2/mageck/count"
OUTPUT_DIR="/global/home/hpc5173/projects/CRISPR_screen/analysis/ht1080_screen_2/mageck_mle"
DESIGN_DIR="/global/home/hpc5173/projects/CRISPR_screen/analysis/ht1080_screen_2/design_matrices"

# Create directories
mkdir -p ${OUTPUT_DIR}
mkdir -p ${DESIGN_DIR}

echo "==========================================="
echo "Starting MAGeCK-MLE Analysis for HT1080 Screen 2"
echo "Date: $(date)"
echo "==========================================="

# Check if count file exists
if [ ! -f "${COUNT_DIR}/ht1080_screen2_R1.count.txt" ]; then
    echo "ERROR: Count file not found: ${COUNT_DIR}/ht1080_screen2_R1.count.txt"
    echo "Please run MAGeCK count first"
    exit 1
fi

# Create design matrix for LOY effect
echo "Creating design matrices..."

# Design 1: Simple LOY effect
cat > ${DESIGN_DIR}/design_loy_effect.txt << EOF
Samples	baseline	loy_effect
D0_CTRL	1	0
D0_LOY	1	1
POST_CTRL_REP1	1	0
POST_CTRL_REP2	1	0
POST_LOY_REP1	1	1
POST_LOY_REP2	1	1
EOF

# Design 2: Time-dependent LOY effect
cat > ${DESIGN_DIR}/design_time_loy.txt << EOF
Samples	baseline	post_selection	loy_effect	loy_post_interaction
D0_CTRL	1	0	0	0
D0_LOY	1	0	1	0
POST_CTRL_REP1	1	1	0	0
POST_CTRL_REP2	1	1	0	0
POST_LOY_REP1	1	1	1	1
POST_LOY_REP2	1	1	1	1
EOF

# Design 3: Batch effect control (if needed)
cat > ${DESIGN_DIR}/design_batch.txt << EOF
Samples	baseline	loy_effect	batch
D0_CTRL	1	0	1
D0_LOY	1	1	1
POST_CTRL_REP1	1	0	2
POST_CTRL_REP2	1	0	2
POST_LOY_REP1	1	1	2
POST_LOY_REP2	1	1	2
EOF

echo "Design matrices created"
echo ""

# Function to run MAGeCK MLE
run_mle() {
    local count_file=$1
    local design_file=$2
    local output_prefix=$3
    local description=$4

    echo "==========================================="
    echo "Running: ${description}"
    echo "Design: ${design_file}"
    echo "==========================================="

    mageck mle \
        -k ${count_file} \
        -d ${design_file} \
        -n ${OUTPUT_DIR}/${output_prefix} \
        --threads ${SLURM_CPUS_PER_TASK} \
        2>&1 | tee ${OUTPUT_DIR}/${output_prefix}_mle.log

    if [ $? -eq 0 ]; then
        echo "✓ ${description} complete"

        # Display summary statistics
        if [ -f "${OUTPUT_DIR}/${output_prefix}.gene_summary.txt" ]; then
            echo ""
            echo "Top 10 genes by LOY effect (beta score):"
            head -n 11 ${OUTPUT_DIR}/${output_prefix}.gene_summary.txt | column -t

            # Count significant genes
            sig_genes=$(awk 'NR>1 && $NF<0.05' ${OUTPUT_DIR}/${output_prefix}.gene_summary.txt | wc -l)
            echo ""
            echo "Significant genes (FDR < 0.05): ${sig_genes}"
        fi
    else
        echo "✗ ${description} failed"
    fi
}

# Run MLE with different designs
echo "Analysis 1: Simple LOY Effect"
run_mle \
    "${COUNT_DIR}/ht1080_screen2_R1.count.txt" \
    "${DESIGN_DIR}/design_loy_effect.txt" \
    "mle_loy_effect_R1" \
    "Simple LOY effect model (R1)"

echo ""
echo "Analysis 2: Time-Dependent LOY Effect"
run_mle \
    "${COUNT_DIR}/ht1080_screen2_R1.count.txt" \
    "${DESIGN_DIR}/design_time_loy.txt" \
    "mle_time_loy_R1" \
    "Time-dependent LOY effect with interaction (R1)"

# Also run for R2 if available
if [ -f "${COUNT_DIR}/ht1080_screen2_R2.count.txt" ]; then
    echo ""
    echo "Running MLE for R2 reads..."

    run_mle \
        "${COUNT_DIR}/ht1080_screen2_R2.count.txt" \
        "${DESIGN_DIR}/design_loy_effect.txt" \
        "mle_loy_effect_R2" \
        "Simple LOY effect model (R2)"
fi

echo ""
echo "==========================================="
echo "Creating MLE Summary Report..."
echo "==========================================="

# Create summary report
SUMMARY_FILE="${OUTPUT_DIR}/mle_analysis_summary.txt"
echo "MAGeCK-MLE Analysis Summary - HT1080 Screen 2" > ${SUMMARY_FILE}
echo "Date: $(date)" >> ${SUMMARY_FILE}
echo "" >> ${SUMMARY_FILE}

# Function to extract top genes
extract_top_genes() {
    local file=$1
    local effect=$2
    local num=10

    if [ -f "${file}" ]; then
        echo "${effect} - Top ${num} genes:" >> ${SUMMARY_FILE}
        echo "-----------------------------------" >> ${SUMMARY_FILE}

        # Get column number for the effect
        col_num=$(head -1 ${file} | tr '\t' '\n' | grep -n "${effect}" | cut -d: -f1)

        if [ ! -z "${col_num}" ]; then
            # Sort by beta score (descending for positive selection)
            tail -n +2 ${file} | sort -t$'\t' -k${col_num} -gr | head -${num} | \
                cut -f1,${col_num},$((col_num+1)) >> ${SUMMARY_FILE}
        fi
        echo "" >> ${SUMMARY_FILE}
    fi
}

# Extract results for each model
echo "Model 1: Simple LOY Effect" >> ${SUMMARY_FILE}
extract_top_genes "${OUTPUT_DIR}/mle_loy_effect_R1.gene_summary.txt" "loy_effect|beta"

echo "Model 2: Time-Dependent Effects" >> ${SUMMARY_FILE}
extract_top_genes "${OUTPUT_DIR}/mle_time_loy_R1.gene_summary.txt" "loy_effect|beta"
extract_top_genes "${OUTPUT_DIR}/mle_time_loy_R1.gene_summary.txt" "loy_post_interaction|beta"

# Compare R1 and R2 if both exist
if [ -f "${OUTPUT_DIR}/mle_loy_effect_R1.gene_summary.txt" ] && \
   [ -f "${OUTPUT_DIR}/mle_loy_effect_R2.gene_summary.txt" ]; then
    echo "" >> ${SUMMARY_FILE}
    echo "R1 vs R2 Concordance:" >> ${SUMMARY_FILE}
    echo "---------------------" >> ${SUMMARY_FILE}

    # Count overlapping significant genes
    awk 'NR>1 && $NF<0.05 {print $1}' ${OUTPUT_DIR}/mle_loy_effect_R1.gene_summary.txt | sort > /tmp/r1_sig.txt
    awk 'NR>1 && $NF<0.05 {print $1}' ${OUTPUT_DIR}/mle_loy_effect_R2.gene_summary.txt | sort > /tmp/r2_sig.txt

    r1_count=$(wc -l < /tmp/r1_sig.txt)
    r2_count=$(wc -l < /tmp/r2_sig.txt)
    overlap=$(comm -12 /tmp/r1_sig.txt /tmp/r2_sig.txt | wc -l)

    echo "R1 significant genes: ${r1_count}" >> ${SUMMARY_FILE}
    echo "R2 significant genes: ${r2_count}" >> ${SUMMARY_FILE}
    echo "Overlapping genes: ${overlap}" >> ${SUMMARY_FILE}

    rm /tmp/r1_sig.txt /tmp/r2_sig.txt
fi

echo ""
echo "Summary saved to: ${SUMMARY_FILE}"

echo ""
echo "==========================================="
echo "MAGeCK-MLE Analysis Complete!"
echo "End time: $(date)"
echo "==========================================="
echo ""
echo "Output files generated:"
echo "  Gene summaries: ${OUTPUT_DIR}/mle_*.gene_summary.txt"
echo "  sgRNA summaries: ${OUTPUT_DIR}/mle_*.sgrna_summary.txt"
echo "  Design matrices: ${DESIGN_DIR}/*.txt"
echo "  Summary report: ${SUMMARY_FILE}"
echo ""
echo "Next steps:"
echo "1. Review gene summary files for beta scores and significance"
echo "2. Compare results across different models"
echo "3. Generate visualizations of MLE results"
echo "4. Validate top hits with known essential genes"