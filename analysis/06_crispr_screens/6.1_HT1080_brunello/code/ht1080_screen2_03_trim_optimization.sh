#!/bin/bash
#SBATCH --job-name=trim_opt_ht1080
#SBATCH --output=/global/home/hpc5173/projects/CRISPR_screen/analysis/ht1080_screen_2/logs/trim_opt_%j.out
#SBATCH --error=/global/home/hpc5173/projects/CRISPR_screen/analysis/ht1080_screen_2/logs/trim_opt_%j.err
#SBATCH --time=3:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G

# HT1080 Screen 2 Trim Optimization Script
# Purpose: Test multiple trim parameters to find optimal settings

# Set up environment
source /global/home/hpc5173/miniconda3/etc/profile.d/conda.sh
conda activate crispr-screen

# Define directories
DATA_DIR="/global/home/hpc5173/projects/CRISPR_screen/data/ht1080_screen_2_raw"
OUTPUT_DIR="/global/home/hpc5173/projects/CRISPR_screen/analysis/ht1080_screen_2/mageck/trim_optimization"
LIBRARY="/global/home/hpc5173/projects/CRISPR_screen/data/brunello/broadgpp-brunello-library-corrected.txt"
CONTROL_SGRNA="/global/home/hpc5173/projects/CRISPR_screen/data/brunello/NTC_control_library.txt"

# Create output directory
mkdir -p ${OUTPUT_DIR}

echo "==========================================="
echo "Starting Trim Parameter Optimization"
echo "Date: $(date)"
echo "==========================================="

# First, subsample the FASTQ files for faster testing
echo "Subsampling FASTQ files for optimization (100,000 reads per file)..."

# Create temporary directory for subsampled files
TEMP_DIR="${OUTPUT_DIR}/subsampled_fastq"
mkdir -p ${TEMP_DIR}

# Subsample one representative file from each condition
for sample in "D0-CTRL_R1_001.fastq.gz" "POST_CTRL_REP1_R1_001.fastq.gz" "POST_LOY_REP1_R1_001.fastq.gz"; do
    # Map original names to simplified names for testing
    if [[ $sample == *"D0-CTRL"* ]]; then
        output_name="D0_CTRL"
    elif [[ $sample == *"CTRL-REP1"* ]]; then
        output_name="POST_CTRL"
    elif [[ $sample == *"LOY-REP1"* ]]; then
        output_name="POST_LOY"
    fi

    echo "Subsampling ${sample}..."
    zcat ${DATA_DIR}/${sample} | head -n 400000 | gzip > ${TEMP_DIR}/${output_name}_R1_subsample.fastq.gz
done

# Also subsample R2 files
for sample in "D0-CTRL_R2_001.fastq.gz" "POST_CTRL_REP1_R2_001.fastq.gz" "POST_LOY_REP1_R2_001.fastq.gz"; do
    if [[ $sample == *"D0-CTRL"* ]]; then
        output_name="D0_CTRL"
    elif [[ $sample == *"CTRL-REP1"* ]]; then
        output_name="POST_CTRL"
    elif [[ $sample == *"LOY-REP1"* ]]; then
        output_name="POST_LOY"
    fi

    echo "Subsampling ${sample}..."
    zcat ${DATA_DIR}/${sample} | head -n 400000 | gzip > ${TEMP_DIR}/${output_name}_R2_subsample.fastq.gz
done

echo "Subsampling complete."
echo ""

# Test different trim parameters
echo "Testing trim parameters..."
echo "==========================================="

# Arrays for trim parameters to test
R1_TRIM_VALUES=(0 1 2 3)
R2_TRIM_VALUES=(20 21 22 23 24 25)

# Create results summary file
RESULTS_FILE="${OUTPUT_DIR}/trim_optimization_results.txt"
echo "Trim Optimization Results - $(date)" > ${RESULTS_FILE}
echo "Sample files: D0_CTRL, POST_CTRL, POST_LOY (100k reads each)" >> ${RESULTS_FILE}
echo "" >> ${RESULTS_FILE}
echo -e "Direction\tTrim_5\tTotal_Reads\tMapped_Reads\tMapping_Rate\tZero_Counts\tGini_Index" >> ${RESULTS_FILE}

# Test R1 trim parameters
echo "Testing R1 trim parameters..."
for trim_val in ${R1_TRIM_VALUES[@]}; do
    echo "  Testing R1 trim-5=${trim_val}..."

    # Create sample list for this test
    echo -e "D0_CTRL\t${TEMP_DIR}/D0_CTRL_R1_subsample.fastq.gz\tbaseline" > ${TEMP_DIR}/sample_list_R1.txt
    echo -e "POST_CTRL\t${TEMP_DIR}/POST_CTRL_R1_subsample.fastq.gz\tcontrol" >> ${TEMP_DIR}/sample_list_R1.txt
    echo -e "POST_LOY\t${TEMP_DIR}/POST_LOY_R1_subsample.fastq.gz\ttreatment" >> ${TEMP_DIR}/sample_list_R1.txt

    # Run MAGeCK count with this trim value
    mageck count \
        -l ${LIBRARY} \
        --sample-label D0_CTRL,POST_CTRL,POST_LOY \
        --fastq ${TEMP_DIR}/D0_CTRL_R1_subsample.fastq.gz ${TEMP_DIR}/POST_CTRL_R1_subsample.fastq.gz ${TEMP_DIR}/POST_LOY_R1_subsample.fastq.gz \
        --trim-5 ${trim_val} \
        --norm-method control \
        --control-sgrna ${CONTROL_SGRNA} \
        -n ${OUTPUT_DIR}/R1_trim${trim_val} \
        2>&1 | tee ${OUTPUT_DIR}/R1_trim${trim_val}.log

    # Extract mapping statistics from log
    total_reads=$(grep "Total:" ${OUTPUT_DIR}/R1_trim${trim_val}.log | tail -1 | awk '{print $2}')
    mapped_reads=$(grep "Mapped:" ${OUTPUT_DIR}/R1_trim${trim_val}.log | tail -1 | awk '{print $2}')
    mapping_rate=$(grep "Percentage:" ${OUTPUT_DIR}/R1_trim${trim_val}.log | tail -1 | awk '{print $2}')
    zero_counts=$(grep "with 0 counts:" ${OUTPUT_DIR}/R1_trim${trim_val}.count_normalized.txt | awk '{print $NF}')
    gini=$(grep "Gini index" ${OUTPUT_DIR}/R1_trim${trim_val}.count_normalized.txt | head -1 | awk '{print $NF}')

    echo -e "R1\t${trim_val}\t${total_reads}\t${mapped_reads}\t${mapping_rate}\t${zero_counts}\t${gini}" >> ${RESULTS_FILE}
done

# Test R2 trim parameters
echo "Testing R2 trim parameters..."
for trim_val in ${R2_TRIM_VALUES[@]}; do
    echo "  Testing R2 trim-5=${trim_val}..."

    # Create sample list for this test
    echo -e "D0_CTRL\t${TEMP_DIR}/D0_CTRL_R2_subsample.fastq.gz\tbaseline" > ${TEMP_DIR}/sample_list_R2.txt
    echo -e "POST_CTRL\t${TEMP_DIR}/POST_CTRL_R2_subsample.fastq.gz\tcontrol" >> ${TEMP_DIR}/sample_list_R2.txt
    echo -e "POST_LOY\t${TEMP_DIR}/POST_LOY_R2_subsample.fastq.gz\ttreatment" >> ${TEMP_DIR}/sample_list_R2.txt

    # Run MAGeCK count with this trim value
    mageck count \
        -l ${LIBRARY} \
        --sample-label D0_CTRL,POST_CTRL,POST_LOY \
        --fastq ${TEMP_DIR}/D0_CTRL_R2_subsample.fastq.gz ${TEMP_DIR}/POST_CTRL_R2_subsample.fastq.gz ${TEMP_DIR}/POST_LOY_R2_subsample.fastq.gz \
        --trim-5 ${trim_val} \
        --norm-method control \
        --control-sgrna ${CONTROL_SGRNA} \
        -n ${OUTPUT_DIR}/R2_trim${trim_val} \
        2>&1 | tee ${OUTPUT_DIR}/R2_trim${trim_val}.log

    # Extract mapping statistics from log
    total_reads=$(grep "Total:" ${OUTPUT_DIR}/R2_trim${trim_val}.log | tail -1 | awk '{print $2}')
    mapped_reads=$(grep "Mapped:" ${OUTPUT_DIR}/R2_trim${trim_val}.log | tail -1 | awk '{print $2}')
    mapping_rate=$(grep "Percentage:" ${OUTPUT_DIR}/R2_trim${trim_val}.log | tail -1 | awk '{print $2}')
    zero_counts=$(grep "with 0 counts:" ${OUTPUT_DIR}/R2_trim${trim_val}.count_normalized.txt | awk '{print $NF}')
    gini=$(grep "Gini index" ${OUTPUT_DIR}/R2_trim${trim_val}.count_normalized.txt | head -1 | awk '{print $NF}')

    echo -e "R2\t${trim_val}\t${total_reads}\t${mapped_reads}\t${mapping_rate}\t${zero_counts}\t${gini}" >> ${RESULTS_FILE}
done

echo ""
echo "==========================================="
echo "Trim Optimization Complete"
echo "==========================================="

# Display results
echo ""
echo "Results Summary:"
cat ${RESULTS_FILE}

# Find best parameters
echo ""
echo "Finding optimal parameters..."
best_r1=$(sort -t$'\t' -k5 -nr ${RESULTS_FILE} | grep "^R1" | head -1)
best_r2=$(sort -t$'\t' -k5 -nr ${RESULTS_FILE} | grep "^R2" | head -1)

echo "Best R1 trim: $(echo ${best_r1} | cut -f2) (mapping rate: $(echo ${best_r1} | cut -f5))"
echo "Best R2 trim: $(echo ${best_r2} | cut -f2) (mapping rate: $(echo ${best_r2} | cut -f5))"

# Save recommendations
echo "" >> ${RESULTS_FILE}
echo "RECOMMENDATIONS:" >> ${RESULTS_FILE}
echo "Best R1 trim-5: $(echo ${best_r1} | cut -f2)" >> ${RESULTS_FILE}
echo "Best R2 trim-5: $(echo ${best_r2} | cut -f2)" >> ${RESULTS_FILE}

echo ""
echo "Results saved to: ${RESULTS_FILE}"
echo ""
echo "Next steps:"
echo "1. Review the optimization results"
echo "2. Run full MAGeCK count with optimized parameters"
echo "   - Use the recommended trim values in the full analysis"