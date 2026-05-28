#!/bin/bash
# ChIP-seq pipeline for single-end reads
# Usage: bash chip_seq_pipeline_SE.sh <reads.fastq.gz> <sample_name> <genome_index> <blacklist_bed>

READS=$1
SAMPLE=$2
GENOME_INDEX=$3
BLACKLIST=$4

# 1. Quality control
fastp -i ${READS} \
      -o ${SAMPLE}_trim.fq.gz \
      -l 25 \
      -h ${SAMPLE}_fastp.html -j ${SAMPLE}_fastp.json

# 2. Genome mapping
hisat2 -X 500 -p 4 \
       --no-spliced-alignment --no-temp-splicesite \
       -x ${GENOME_INDEX} \
       -U ${SAMPLE}_trim.fq.gz \
       -3 1 --summary-file ${SAMPLE}_aln_sum.txt | \
       samtools view -ShuF 4 -q 30 - | \
       samtools sort -o ${SAMPLE}.sorted.bam

# 3. Blacklist removal
intersectBed -a ${SAMPLE}.sorted.bam -b ${BLACKLIST} -v > ${SAMPLE}.final.bam

# 4. PCR duplicate removal
picard MarkDuplicates \
       INPUT=${SAMPLE}.final.bam \
       OUTPUT=${SAMPLE}.pmd.bam \
       REMOVE_DUPLICATES=true ASSUME_SORTED=true \
       METRICS_FILE=${SAMPLE}_pmd.out \
       2>${SAMPLE}_picard.log

echo "Pipeline finished for ${SAMPLE}. Output: ${SAMPLE}.pmd.bam"