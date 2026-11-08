#!/bin/bash
# Simple test of different trim parameters

cd /global/home/hpc5173/projects/CRISPR_screen
source activate crispr-screen

# Create test directory
mkdir -p analysis/ht1080_screen_2/trim_test

# Extract 100k reads for testing
echo "Extracting test reads..."
zcat data/ht1080_screen_2_raw/D0-CTRL_R1_001.fastq.gz | head -400000 | gzip > analysis/ht1080_screen_2/trim_test/d0_test.fastq.gz
zcat data/ht1080_screen_2_raw/CTRL-REP1_R1_001.fastq.gz | head -400000 | gzip > analysis/ht1080_screen_2/trim_test/post_test.fastq.gz

echo ""
echo "Testing different trim-5 values on D0 sample..."
echo "==============================================="

for trim in 0 1 2 20 21 22 23 24 25; do
    echo ""
    echo "Testing trim-5=$trim on D0..."
    mageck count \
        -l data/brunello/broadgpp-brunello-library-corrected.txt \
        --sample-label D0_trim${trim} \
        --fastq analysis/ht1080_screen_2/trim_test/d0_test.fastq.gz \
        --trim-5 $trim \
        -n analysis/ht1080_screen_2/trim_test/d0_trim${trim} \
        2>&1 | grep -E "Total:|Mapped:" | tail -2
done

echo ""
echo ""
echo "Testing different trim-5 values on POST sample..."
echo "==============================================="

for trim in 0 1 2 20 21 22 23 24 25; do
    echo ""
    echo "Testing trim-5=$trim on POST..."
    mageck count \
        -l data/brunello/broadgpp-brunello-library-corrected.txt \
        --sample-label POST_trim${trim} \
        --fastq analysis/ht1080_screen_2/trim_test/post_test.fastq.gz \
        --trim-5 $trim \
        -n analysis/ht1080_screen_2/trim_test/post_trim${trim} \
        2>&1 | grep -E "Total:|Mapped:" | tail -2
done

echo ""
echo "Analysis complete. Check log files for details."