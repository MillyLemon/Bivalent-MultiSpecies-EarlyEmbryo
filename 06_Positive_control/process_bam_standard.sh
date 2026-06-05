#!/bin/bash
# Standard BAM processing: blacklist removal + PCR duplicate removal
# Usage: bash process_bam_standard.sh <sample_name> <genome_blacklist>
# Example: bash process_bam_standard.sh GSM691591_ES.H3K27 /path/to/mm10.blacklist.bed

SAMPLE=$1
BLACKLIST=$2

echo "=== Processing ${SAMPLE} ==="

# 1. Blacklist removal
echo "Removing blacklist regions..."
intersectBed -a ${SAMPLE}.fixed.bam -b ${BLACKLIST} -v > ${SAMPLE}.final.bam

# 2. PCR duplicate removal
echo "Removing PCR duplicates..."
picard MarkDuplicates \
    INPUT=${SAMPLE}.final.bam \
    OUTPUT=${SAMPLE}.pmd.bam \
    REMOVE_DUPLICATES=true \
    ASSUME_SORTED=true \
    METRICS_FILE=${SAMPLE}.pmd.out \
    2>${SAMPLE}.picard.log

echo "Done! Output: ${SAMPLE}.pmd.bam"