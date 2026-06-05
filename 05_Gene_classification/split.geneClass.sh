#!/bin/bash

# 参考BED文件
#REF_BED="mm10_ensGene_WithStrand.bed"
REF_BED="gencode.vM23.annotation.genes.bed"

# 遍历所有样本目录
for sample_dir in */; do
    # 去除目录名的斜杠
    sample_dir=${sample_dir%/}
    result_file="${sample_dir}/Gene_classification.txt"
    
    if [[ -f "$result_file" ]]; then
        echo "处理样本: $sample_dir"
        
        # 为当前样本创建三个BED文件
        broad_bed="${sample_dir}/broad_genes.bed"
        promoter_bed="${sample_dir}/promoter_genes.bed" 
        tss_bed="${sample_dir}/tss_genes.bed"
        
        # 清空输出文件
        > "$broad_bed"
        > "$promoter_bed"
        > "$tss_bed"
        
        # 提取Broad基因
        awk '$4 == "Broad" {print $1}' "$result_file" | while read gene_id; do
            grep "$gene_id" "$REF_BED" >> "$broad_bed"
        done
        
        # 提取Promoter基因  
        awk '$4 == "Promoter" {print $1}' "$result_file" | while read gene_id; do
            grep "$gene_id" "$REF_BED" >> "$promoter_bed"
        done
        
        # 提取TSS基因
        awk '$4 == "TSS" {print $1}' "$result_file" | while read gene_id; do
            grep "$gene_id" "$REF_BED" >> "$tss_bed"
        done
        
        echo "完成: $sample_dir - Broad:$(wc -l < "$broad_bed") Promoter:$(wc -l < "$promoter_bed") TSS:$(wc -l < "$tss_bed")"
    fi
done
