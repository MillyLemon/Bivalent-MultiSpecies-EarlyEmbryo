#!/bin/bash
# Merge biological replicates for ChIP-seq
# Usage: bash merge_replicates.sh <output_name> <rep1.bam> <rep2.bam> [rep3.bam ...]

OUTPUT=$1
shift  # 移除第一个参数，剩余的都是输入BAM文件

samtools merge -@ 16 ${OUTPUT}.bam $@
samtools index ${OUTPUT}.bam

echo "Merged replicates into ${OUTPUT}.bam"