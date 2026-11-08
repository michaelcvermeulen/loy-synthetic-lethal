#!/bin/bash
#SBATCH --job-name=export_ht1080
#SBATCH --output=/global/home/hpc5173/projects/CRISPR_screen/analysis/ht1080_screen_2/logs/export_%j.out
#SBATCH --error=/global/home/hpc5173/projects/CRISPR_screen/analysis/ht1080_screen_2/logs/export_%j.err
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G

# HT1080 Screen 2 - Data Export Package Script
# Purpose: Create a comprehensive data package for sharing and local analysis

# Set up environment
source /global/home/hpc5173/miniconda3/etc/profile.d/conda.sh
conda activate crispr-screen

# Define directories
BASE_DIR="/global/home/hpc5173/projects/CRISPR_screen/analysis/ht1080_screen_2"
EXPORT_DIR="${BASE_DIR}/export_package"
DATE_STAMP=$(date +%Y%m%d)
PACKAGE_NAME="ht1080_screen2_results_${DATE_STAMP}"

echo "==========================================="
echo "Creating Export Package for HT1080 Screen 2"
echo "Date: $(date)"
echo "Package: ${PACKAGE_NAME}"
echo "==========================================="

# Create export directory structure
echo "Creating directory structure..."
mkdir -p ${EXPORT_DIR}/${PACKAGE_NAME}/{counts,qc,results,metadata,scripts,reports,visualization}

# Copy count matrices
echo "Copying count matrices..."
if [ -d "${BASE_DIR}/mageck/count" ]; then
    cp ${BASE_DIR}/mageck/count/*.count.txt ${EXPORT_DIR}/${PACKAGE_NAME}/counts/ 2>/dev/null
    cp ${BASE_DIR}/mageck/count/*.count_normalized.txt ${EXPORT_DIR}/${PACKAGE_NAME}/counts/ 2>/dev/null
    echo "  - Count matrices copied"
else
    echo "  - No count matrices found"
fi

# Copy QC reports
echo "Copying QC reports..."
if [ -d "${BASE_DIR}/fastqc" ]; then
    # Copy only the HTML reports (not the zip files to save space)
    cp ${BASE_DIR}/fastqc/*.html ${EXPORT_DIR}/${PACKAGE_NAME}/qc/ 2>/dev/null
    echo "  - FastQC reports copied"
fi

if [ -d "${BASE_DIR}/multiqc" ]; then
    cp ${BASE_DIR}/multiqc/*.html ${EXPORT_DIR}/${PACKAGE_NAME}/qc/ 2>/dev/null
    echo "  - MultiQC report copied"
fi

if [ -d "${BASE_DIR}/qc" ]; then
    cp ${BASE_DIR}/qc/*.txt ${EXPORT_DIR}/${PACKAGE_NAME}/qc/ 2>/dev/null
    cp ${BASE_DIR}/qc/*.pdf ${EXPORT_DIR}/${PACKAGE_NAME}/qc/ 2>/dev/null
    echo "  - Count QC files copied"
fi

# Copy MAGeCK results
echo "Copying MAGeCK results..."
if [ -d "${BASE_DIR}/mageck/test" ]; then
    cp ${BASE_DIR}/mageck/test/*.gene_summary.txt ${EXPORT_DIR}/${PACKAGE_NAME}/results/ 2>/dev/null
    cp ${BASE_DIR}/mageck/test/*.sgrna_summary.txt ${EXPORT_DIR}/${PACKAGE_NAME}/results/ 2>/dev/null
    cp ${BASE_DIR}/mageck/test/*.pdf ${EXPORT_DIR}/${PACKAGE_NAME}/results/ 2>/dev/null
    echo "  - MAGeCK test results copied"
fi

if [ -d "${BASE_DIR}/mageck_mle" ]; then
    cp ${BASE_DIR}/mageck_mle/*.gene_summary.txt ${EXPORT_DIR}/${PACKAGE_NAME}/results/ 2>/dev/null
    cp ${BASE_DIR}/mageck_mle/*.sgrna_summary.txt ${EXPORT_DIR}/${PACKAGE_NAME}/results/ 2>/dev/null
    echo "  - MAGeCK-MLE results copied"
fi

# Copy visualization plots
echo "Copying visualization plots..."
if [ -d "${BASE_DIR}/visualization" ]; then
    cp ${BASE_DIR}/visualization/*.pdf ${EXPORT_DIR}/${PACKAGE_NAME}/visualization/ 2>/dev/null
    echo "  - Visualization plots copied"
fi

# Copy metadata
echo "Copying metadata files..."
cp ${BASE_DIR}/sample_metadata.txt ${EXPORT_DIR}/${PACKAGE_NAME}/metadata/ 2>/dev/null
cp ${BASE_DIR}/mageck_labels.txt ${EXPORT_DIR}/${PACKAGE_NAME}/metadata/ 2>/dev/null
cp ${BASE_DIR}/PIPELINE_STATUS.md ${EXPORT_DIR}/${PACKAGE_NAME}/metadata/ 2>/dev/null

# Copy analysis scripts
echo "Copying analysis scripts..."
cp /global/home/hpc5173/projects/CRISPR_screen/code/ht1080_screen2_*.sh ${EXPORT_DIR}/${PACKAGE_NAME}/scripts/
cp /global/home/hpc5173/projects/CRISPR_screen/code/ht1080_screen2_*.R ${EXPORT_DIR}/${PACKAGE_NAME}/scripts/

# Create summary report
echo "Creating summary report..."
REPORT_FILE="${EXPORT_DIR}/${PACKAGE_NAME}/reports/ANALYSIS_SUMMARY.md"

cat > ${REPORT_FILE} << EOF
# HT1080 Screen 2 - Analysis Summary

## Overview
- **Analysis Date**: $(date)
- **Dataset**: HT1080 CRISPR Screen (Screen 2)
- **Library**: Brunello (77,441 sgRNAs, 19,114 genes)
- **Samples**: 6 (D0_CTRL, D0_LOY, POST_CTRL_REP1/2, POST_LOY_REP1/2)

## Pipeline Steps Completed

### 1. Quality Control
- FastQC analysis of all 12 FASTQ files
- MultiQC aggregated report
- Count matrix QC metrics

### 2. sgRNA Counting
- MAGeCK count with optimized trim parameters
- R1 trim-5: ${R1_TRIM:-TBD}
- R2 trim-5: ${R2_TRIM:-TBD}

### 3. Differential Analysis
- D0 LOY vs CTRL comparison
- POST LOY vs CTRL comparison
- Temporal analysis (POST vs D0)

### 4. Advanced Analysis
- MAGeCK-MLE with multiple models
- Visualization of results

## Key Findings
EOF

# Add statistics if available
if [ -f "${BASE_DIR}/mageck/test/differential_analysis_summary.txt" ]; then
    echo "" >> ${REPORT_FILE}
    echo "### Differential Analysis Results" >> ${REPORT_FILE}
    echo '```' >> ${REPORT_FILE}
    cat ${BASE_DIR}/mageck/test/differential_analysis_summary.txt >> ${REPORT_FILE}
    echo '```' >> ${REPORT_FILE}
fi

# Add file inventory
echo "" >> ${REPORT_FILE}
echo "## File Inventory" >> ${REPORT_FILE}
echo "" >> ${REPORT_FILE}

for dir in counts qc results metadata scripts reports visualization; do
    count=$(ls -1 ${EXPORT_DIR}/${PACKAGE_NAME}/${dir}/ 2>/dev/null | wc -l)
    echo "- **${dir}**: ${count} files" >> ${REPORT_FILE}
done

echo "" >> ${REPORT_FILE}
echo "## Usage Instructions" >> ${REPORT_FILE}
echo "" >> ${REPORT_FILE}
echo "1. Extract the archive to your local machine" >> ${REPORT_FILE}
echo "2. Count matrices are in \`counts/\` directory" >> ${REPORT_FILE}
echo "3. Differential analysis results are in \`results/\` directory" >> ${REPORT_FILE}
echo "4. QC reports are in \`qc/\` directory" >> ${REPORT_FILE}
echo "5. Visualization plots are in \`visualization/\` directory" >> ${REPORT_FILE}
echo "6. Analysis scripts are in \`scripts/\` for reproducibility" >> ${REPORT_FILE}

echo "Summary report created"

# Create README
echo "Creating README..."
cat > ${EXPORT_DIR}/${PACKAGE_NAME}/README.txt << EOF
HT1080 CRISPR Screen 2 - Analysis Results
==========================================

This package contains the complete analysis results for the HT1080 CRISPR screen (Screen 2).

Directory Structure:
--------------------
- counts/         : Raw and normalized count matrices from MAGeCK
- qc/            : Quality control reports (FastQC, MultiQC, count QC)
- results/       : Differential analysis results (MAGeCK test and MLE)
- metadata/      : Sample information and experimental design
- scripts/       : Analysis scripts used in the pipeline
- reports/       : Summary reports and documentation
- visualization/ : Plots and figures

Key Files:
----------
- counts/ht1080_screen2_R1.count.txt : Raw count matrix (R1 reads)
- results/*gene_summary.txt : Gene-level differential analysis results
- qc/ht1080_screen2_multiqc_report.html : Aggregated QC report
- reports/ANALYSIS_SUMMARY.md : Complete analysis summary

Contact:
--------
Analysis performed by: Claude (AI Assistant)
Date: $(date)
Project: CRISPR_screen

For questions about the analysis, refer to the scripts/ directory
for complete reproducibility information.
EOF

echo "README created"

# Create compressed archive
echo ""
echo "Creating compressed archive..."
cd ${EXPORT_DIR}
tar -czf ${PACKAGE_NAME}.tar.gz ${PACKAGE_NAME}/

# Get file sizes
PACKAGE_SIZE=$(du -h ${PACKAGE_NAME}.tar.gz | cut -f1)
UNCOMPRESSED_SIZE=$(du -sh ${PACKAGE_NAME}/ | cut -f1)

echo ""
echo "==========================================="
echo "Export Package Creation Complete!"
echo "==========================================="
echo ""
echo "Package details:"
echo "  Name: ${PACKAGE_NAME}.tar.gz"
echo "  Location: ${EXPORT_DIR}/"
echo "  Compressed size: ${PACKAGE_SIZE}"
echo "  Uncompressed size: ${UNCOMPRESSED_SIZE}"
echo ""
echo "File counts by category:"
for dir in counts qc results metadata scripts reports visualization; do
    count=$(ls -1 ${EXPORT_DIR}/${PACKAGE_NAME}/${dir}/ 2>/dev/null | wc -l)
    printf "  %-15s : %3d files\n" ${dir} ${count}
done
echo ""
echo "To download the package, use:"
echo "  scp ${USER}@$(hostname):${EXPORT_DIR}/${PACKAGE_NAME}.tar.gz ."
echo ""
echo "Export complete at: $(date)"