#!/bin/bash
# Merge fastq files from multiple single cells at the same developmental stage
# Usage: bash merge_fastq.sh <stage_name> <SRR1_R1.fq.gz,SRR2_R1.fq.gz,...> <SRR1_R2.fq.gz,SRR2_R2.fq.gz,...>
#
# Example:
#   bash merge_fastq.sh WT_2-cell "SRR2089603,SRR2089604,SRR2089605,SRR2089606,SRR2089607" "R1" "R2"
#
# Note: This script is provided as a reference. In practice, you may need to adapt
# the file paths and naming conventions according to your data organization.

STAGE=$1
R1_FILES=$2   # Comma-separated list of R1 fastq files
R2_FILES=$3   # Comma-separated list of R2 fastq files

# Convert comma-separated strings to space-separated
R1_LIST=$(echo $R1_FILES | tr ',' ' ')
R2_LIST=$(echo $R2_FILES | tr ',' ' ')

echo "Merging R1 files for ${STAGE}..."
cat ${R1_LIST} > ${STAGE}_1.fastq.gz

echo "Merging R2 files for ${STAGE}..."
cat ${R2_LIST} > ${STAGE}_2.fastq.gz

echo "Merged files: ${STAGE}_1.fastq.gz, ${STAGE}_2.fastq.gz"
