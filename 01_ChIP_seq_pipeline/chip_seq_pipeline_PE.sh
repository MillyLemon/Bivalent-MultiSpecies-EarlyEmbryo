#!/bin/bash
# ChIP-seq pipeline for paired-end reads
# Usage: bash chip_seq_pipeline_PE.sh <R1.fastq.gz> <R2.fastq.gz> <sample_name> <genome_index> <blacklist_bed>

R1=$1
R2=$2
SAMPLE=$3
GENOME_INDEX=$4
BLACKLIST=$5

# 1. Quality control
fastp -i ${R1} -I ${R2} \
      -o ${SAMPLE}_1_trim.fq.gz -O ${SAMPLE}_2_trim.fq.gz \
      --detect_adapter_for_pe -l 25 \
      -h ${SAMPLE}_fastp.html -j ${SAMPLE}_fastp.json

# 2. Genome mapping
hisat2 -X 500 -p 4 \
       --no-spliced-alignment --no-temp-splicesite \
       -x ${GENOME_INDEX} \
       -1 ${SAMPLE}_1_trim.fq.gz -2 ${SAMPLE}_2_trim.fq.gz \
       -3 1 --summary-file ${SAMPLE}_aln_sum.txt | \
       samtools view -ShuF 4 -f 2 -q 30 - | \
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