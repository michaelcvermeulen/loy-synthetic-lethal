#!/bin/bash
#SBATCH --job-name=mageck_mle_ht1080
#SBATCH --output=/global/home/hpc5173/projects/CRISPR_screen/analysis/ht1080_screen_2/logs/mageck_mle_%j.out
#SBATCH --error=/global/home/hpc5173/projects/CRISPR_screen/analysis/ht1080_screen_2/logs/mageck_mle_%j.err
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

# HT1080 Screen 2 MAGeCK-MLE Script
# Purpose: Advanced analysis using Maximum Likelihood Estimation

# Set up environment
eval "$(conda shell.bash hook)"
conda activate crispr-screen

# Define directories
COUNT_DIR="/global/home/hpc5173/projects/CRISPR_screen/analysis/ht1080_screen_2/mageck/count"
MLE_DIR="/global/home/hpc5173/projects/CRISPR_screen/analysis/ht1080_screen_2/mageck/mle"
DESIGN_DIR="/global/home/hpc5173/projects/CRISPR_screen/data/mageck/design_matrices"

# Create output and design directories
mkdir -p ${MLE_DIR}
mkdir -p ${DESIGN_DIR}

echo "==========================================="
echo "Starting MAGeCK-MLE for HT1080 Screen 2"
echo "Date: $(date)"
echo "==========================================="
echo ""

# Check if count matrix exists
if [ ! -f "${COUNT_DIR}/ht1080_screen2_R1.count.txt" ]; then
    echo "ERROR: Count matrix not found!"
    echo "Please run MAGeCK count first."
    exit 1
fi

# Create design matrices
echo "Creating design matrices..."

# 1. Simple LOY effect design matrix
cat > "${DESIGN_DIR}/design_loy_effect.txt" << EOF
sample	baseline	loy_treatment
D0_CTRL	1	0
D0_LOY	1	0
POST_CTRL_REP1	1	0
POST_CTRL_REP2	1	0
POST_LOY_REP1	1	1
POST_LOY_REP2	1	1
EOF

# 2. Time and treatment interaction design matrix
cat > "${DESIGN_DIR}/design_time_treatment.txt" << EOF
sample	baseline	time	loy	time_loy_interaction
D0_CTRL	1	0	0	0
D0_LOY	1	0	1	0
POST_CTRL_REP1	1	1	0	0
POST_CTRL_REP2	1	1	0	0
POST_LOY_REP1	1	1	1	1
POST_LOY_REP2	1	1	1	1
EOF

# 3. Batch effect control (if needed)
cat > "${DESIGN_DIR}/design_with_batch.txt" << EOF
sample	baseline	loy_treatment	replicate
D0_CTRL	1	0	0
D0_LOY	1	0	0
POST_CTRL_REP1	1	0	1
POST_CTRL_REP2	1	0	2
POST_LOY_REP1	1	1	1
POST_LOY_REP2	1	1	2
EOF

echo "Design matrices created."
echo ""

# Run MAGeCK-MLE with different models
echo "Running MAGeCK-MLE analyses..."
echo ""

# 1. Simple LOY effect model
echo "1. Simple LOY effect model..."
mageck mle \
    -k "${COUNT_DIR}/ht1080_screen2_R1.count.txt" \
    -d "${DESIGN_DIR}/design_loy_effect.txt" \
    --norm-method control \
    --control-sgrna /global/home/hpc5173/projects/CRISPR_screen/data/brunello/NTC_control_library.txt \
    --threads 8 \
    -n "${MLE_DIR}/mle_loy_effect" \
    2>&1 | tee "${MLE_DIR}/mle_loy_effect.log"

echo ""

# 2. Time and treatment interaction model
echo "2. Time and treatment interaction model..."
mageck mle \
    -k "${COUNT_DIR}/ht1080_screen2_R1.count.txt" \
    -d "${DESIGN_DIR}/design_time_treatment.txt" \
    --norm-method control \
    --control-sgrna /global/home/hpc5173/projects/CRISPR_screen/data/brunello/NTC_control_library.txt \
    --threads 8 \
    -n "${MLE_DIR}/mle_time_treatment" \
    2>&1 | tee "${MLE_DIR}/mle_time_treatment.log"

echo ""

# 3. Model with batch effects (if relevant)
echo "3. Model with batch effects..."
mageck mle \
    -k "${COUNT_DIR}/ht1080_screen2_R1.count.txt" \
    -d "${DESIGN_DIR}/design_with_batch.txt" \
    --norm-method control \
    --control-sgrna /global/home/hpc5173/projects/CRISPR_screen/data/brunello/NTC_control_library.txt \
    --threads 8 \
    -n "${MLE_DIR}/mle_with_batch" \
    2>&1 | tee "${MLE_DIR}/mle_with_batch.log"

echo ""
echo "==========================================="
echo "MAGeCK-MLE Complete!"
echo "End time: $(date)"
echo "==========================================="
echo ""

# Summary of results
echo "Summary of MLE results:"
echo "-----------------------"

for model in loy_effect time_treatment with_batch; do
    result_file="${MLE_DIR}/mle_${model}.gene_summary.txt"
    if [ -f "$result_file" ]; then
        echo ""
        echo "${model} model:"

        # Count significant genes for each coefficient
        n_cols=$(head -1 "$result_file" | awk '{print NF}')

        # Parse beta coefficients and FDR columns
        head -1 "$result_file" | tr '\t' '\n' | nl | grep -E "beta|fdr" | while read num col; do
            if [[ "$col" == *"beta"* ]]; then
                beta_col=$num
                effect=$(echo "$col" | sed 's/|beta//g')

                # Find corresponding FDR column
                fdr_col=$(head -1 "$result_file" | tr '\t' '\n' | nl | grep "${effect}|fdr" | awk '{print $1}')

                if [ ! -z "$fdr_col" ]; then
                    # Count significant positive and negative effects
                    n_pos=$(awk -v bc=$beta_col -v fc=$fdr_col 'NR>1 && $fc < 0.1 && $bc > 0' "$result_file" | wc -l)
                    n_neg=$(awk -v bc=$beta_col -v fc=$fdr_col 'NR>1 && $fc < 0.1 && $bc < 0' "$result_file" | wc -l)

                    echo "  ${effect}:"
                    echo "    Significant positive effects (FDR < 0.1): ${n_pos}"
                    echo "    Significant negative effects (FDR < 0.1): ${n_neg}"

                    # Show top 3 hits
                    echo "    Top 3 effects:"
                    awk -v bc=$beta_col -v fc=$fdr_col 'NR>1 && $fc < 0.1' "$result_file" | \
                        sort -k${fdr_col},${fdr_col}g | head -3 | \
                        awk -v bc=$beta_col -v fc=$fdr_col '{print "      " $1 " (beta=" $bc ", FDR=" $fc ")"}'
                fi
            fi
        done
    fi
done

echo ""
echo "Output files generated in: ${MLE_DIR}"
echo ""
echo "Key output files:"
echo "  - *gene_summary.txt: Gene-level beta scores and significance"
echo "  - *sgrna_summary.txt: sgRNA-level statistics"
echo ""
echo "Design matrices saved in: ${DESIGN_DIR}"
echo ""
echo "Next steps:"
echo "1. Review gene_summary.txt files for each model"
echo "2. Compare results between different models"
echo "3. Perform pathway enrichment analysis on significant genes"
echo "4. Generate visualizations (volcano plots, heatmaps)"
echo "5. Export results for downstream analysis"