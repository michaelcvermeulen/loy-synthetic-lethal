#!/bin/bash

# Script to check which read (R1 or R2) contains the sgRNA barcodes
# This will help us identify which read to use for the actual analysis

echo "==========================================="
echo "Checking reads for sgRNA barcodes"
echo "Date: $(date)"
echo "==========================================="

# Define directories
DATA_DIR="/global/home/hpc5173/projects/CRISPR_screen/data/ht1080_screen_2_raw"
LIBRARY="/global/home/hpc5173/projects/CRISPR_screen/data/brunello/broadgpp-brunello-library-corrected.txt"
OUTPUT_DIR="/global/home/hpc5173/projects/CRISPR_screen/analysis/ht1080_screen_2/barcode_check"

# Create output directory
mkdir -p ${OUTPUT_DIR}

# Get a few example sgRNA sequences from the library
echo "Extracting example sgRNA sequences from library..."
head -20 ${LIBRARY} | tail -10 | cut -f2 > ${OUTPUT_DIR}/example_sgrnas.txt

echo "Example sgRNA sequences (20bp):"
cat ${OUTPUT_DIR}/example_sgrnas.txt

# Function to check reads for sgRNA presence
check_read_for_sgrna() {
    local fastq_file=$1
    local read_type=$2
    local sample_name=$3

    echo ""
    echo "----------------------------------------"
    echo "Checking: ${sample_name} - ${read_type}"
    echo "File: $(basename ${fastq_file})"

    # Extract first 10000 reads
    echo "Extracting first 10000 reads..."
    zcat ${fastq_file} | head -40000 | grep -A1 "^@" | grep -v "^@" | grep -v "^--" | head -1000 > ${OUTPUT_DIR}/${sample_name}_${read_type}_sequences.txt

    # Check different positions in the read for sgRNA matches
    echo "Checking for sgRNA sequences at different positions..."

    # Check position 0-20 (start of read)
    echo "  Position 0-20:"
    cut -c1-20 ${OUTPUT_DIR}/${sample_name}_${read_type}_sequences.txt | head -100 > ${OUTPUT_DIR}/${sample_name}_${read_type}_pos0-20.txt
    matches_0_20=$(grep -f ${OUTPUT_DIR}/example_sgrnas.txt ${OUTPUT_DIR}/${sample_name}_${read_type}_pos0-20.txt | wc -l)
    echo "    Matches found: ${matches_0_20}/100"

    # Check position 1-21 (offset by 1)
    echo "  Position 1-21:"
    cut -c2-21 ${OUTPUT_DIR}/${sample_name}_${read_type}_sequences.txt | head -100 > ${OUTPUT_DIR}/${sample_name}_${read_type}_pos1-21.txt
    matches_1_21=$(grep -f ${OUTPUT_DIR}/example_sgrnas.txt ${OUTPUT_DIR}/${sample_name}_${read_type}_pos1-21.txt | wc -l)
    echo "    Matches found: ${matches_1_21}/100"

    # Check position 23-43 (common for R2 with adapter)
    echo "  Position 23-43:"
    cut -c24-43 ${OUTPUT_DIR}/${sample_name}_${read_type}_sequences.txt | head -100 > ${OUTPUT_DIR}/${sample_name}_${read_type}_pos23-43.txt
    matches_23_43=$(grep -f ${OUTPUT_DIR}/example_sgrnas.txt ${OUTPUT_DIR}/${sample_name}_${read_type}_pos23-43.txt | wc -l)
    echo "    Matches found: ${matches_23_43}/100"

    # Show first 5 sequences from each position for visual inspection
    echo ""
    echo "Sample sequences from position 0-20:"
    head -5 ${OUTPUT_DIR}/${sample_name}_${read_type}_pos0-20.txt

    echo ""
    echo "Sample sequences from position 1-21:"
    head -5 ${OUTPUT_DIR}/${sample_name}_${read_type}_pos1-21.txt

    echo ""
    echo "Sample sequences from position 23-43:"
    head -5 ${OUTPUT_DIR}/${sample_name}_${read_type}_pos23-43.txt

    # Check for common CRISPR vector sequences
    echo ""
    echo "Checking for common vector elements..."

    # Check for U6 promoter remnants (often GGAAA or similar at start)
    echo "  Reads starting with G: $(grep "^G" ${OUTPUT_DIR}/${sample_name}_${read_type}_sequences.txt | head -100 | wc -l)/100"
    echo "  Reads starting with GG: $(grep "^GG" ${OUTPUT_DIR}/${sample_name}_${read_type}_sequences.txt | head -100 | wc -l)/100"

    # Check read length distribution
    echo ""
    echo "Read length (first 10 reads):"
    head -10 ${OUTPUT_DIR}/${sample_name}_${read_type}_sequences.txt | awk '{print length($0)}'
}

# Check one representative sample from each condition
echo ""
echo "=========================================="
echo "ANALYZING R1 READS"
echo "=========================================="

# Check D0-CTRL R1
check_read_for_sgrna "${DATA_DIR}/D0-CTRL_R1_001.fastq.gz" "R1" "D0-CTRL"

# Check POST sample R1
check_read_for_sgrna "${DATA_DIR}/CTRL-REP1_R1_001.fastq.gz" "R1" "CTRL-REP1"

echo ""
echo "=========================================="
echo "ANALYZING R2 READS"
echo "=========================================="

# Check D0-CTRL R2
check_read_for_sgrna "${DATA_DIR}/D0-CTRL_R2_001.fastq.gz" "R2" "D0-CTRL"

# Check POST sample R2
check_read_for_sgrna "${DATA_DIR}/CTRL-REP1_R2_001.fastq.gz" "R2" "CTRL-REP1"

# Now let's do a more comprehensive check using a larger sample
echo ""
echo "=========================================="
echo "COMPREHENSIVE sgRNA MAPPING TEST"
echo "=========================================="

# Test mapping with different trim parameters
echo "Testing mapping rates with 100,000 reads..."

# Create a subsample
echo "Creating subsampled files..."
zcat ${DATA_DIR}/D0-CTRL_R1_001.fastq.gz | head -400000 > ${OUTPUT_DIR}/test_R1.fastq
zcat ${DATA_DIR}/D0-CTRL_R2_001.fastq.gz | head -400000 > ${OUTPUT_DIR}/test_R2.fastq

# Activate conda environment for MAGeCK
source /global/home/hpc5173/miniconda3/etc/profile.d/conda.sh
conda activate crispr-screen

# Test R1 with different trim values
echo ""
echo "Testing R1 mapping rates:"
for trim in 0 1 2 3; do
    echo -n "  Trim-5=${trim}: "
    mageck count \
        -l ${LIBRARY} \
        --fastq ${OUTPUT_DIR}/test_R1.fastq \
        --trim-5 ${trim} \
        -n ${OUTPUT_DIR}/test_R1_trim${trim} \
        2>&1 | grep "Percentage:" | tail -1 | awk '{print $2}'
done

# Test R2 with different trim values
echo ""
echo "Testing R2 mapping rates:"
for trim in 0 1 20 21 22 23 24 25; do
    echo -n "  Trim-5=${trim}: "
    mageck count \
        -l ${LIBRARY} \
        --fastq ${OUTPUT_DIR}/test_R2.fastq \
        --trim-5 ${trim} \
        -n ${OUTPUT_DIR}/test_R2_trim${trim} \
        2>&1 | grep "Percentage:" | tail -1 | awk '{print $2}'
done

# Create summary
echo ""
echo "=========================================="
echo "SUMMARY"
echo "=========================================="

SUMMARY_FILE="${OUTPUT_DIR}/barcode_check_summary.txt"
echo "Barcode Check Summary - $(date)" > ${SUMMARY_FILE}
echo "" >> ${SUMMARY_FILE}
echo "Based on the analysis above, determine:" >> ${SUMMARY_FILE}
echo "1. Which read (R1 or R2) contains the sgRNA sequences" >> ${SUMMARY_FILE}
echo "2. What trim position gives the best mapping rate" >> ${SUMMARY_FILE}
echo "3. Whether we need to process both reads or just one" >> ${SUMMARY_FILE}
echo "" >> ${SUMMARY_FILE}
echo "Check the mapping rates and sequence patterns above to make this determination." >> ${SUMMARY_FILE}

cat ${SUMMARY_FILE}

echo ""
echo "Analysis complete. Check the output above to determine which read to use."
echo "Results saved to: ${OUTPUT_DIR}"