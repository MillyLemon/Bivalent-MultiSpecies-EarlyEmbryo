#!/bin/bash

dir="/data/lixm/My_project/Analysis/mouse/paper1/ChIP_seq/picard_bam/merged_bams"

# 遍历所有样本
for sample_dir in  `cat sample.name`; do
    sample=${sample_dir%/}
    echo "处理样本: $sample"

 # 检查文件是否存在
    echo "检查BigWig文件:"
    ls "${dir}/merged_${sample}_uniform_25bp_CPM.bw"
    echo "检查BED文件:"
    ls "${sample}/broad_genes.trans.bed" "${sample}/promoter_genes.trans.bed" "${sample}/tss_genes.trans.bed"

    # TSS 分析
    computeMatrix reference-point \
        --referencePoint TSS \
        -S "${dir}/merged_${sample}_uniform_25bp_CPM.bw" \
        -R "${sample}/broad_genes.trans.bed" "${sample}/promoter_genes.trans.bed" "${sample}/tss_genes.trans.bed" \
        -b 30000 -a 30000 \
        --skipZeros \
        -o ${sample}/${sample}_matrix_trans_TSS.gz
    
    # TES 分析  
    computeMatrix scale-regions \
        -S "${dir}/merged_${sample}_uniform_25bp_CPM.bw" \
        -R "${sample}/broad_genes.trans.bed" "${sample}/promoter_genes.trans.bed" "${sample}/tss_genes.trans.bed" \
        -b 30000 -a 30000 \
        --binSize  25  \
        --regionBodyLength 30000 \
        --skipZeros \
        -o ${sample}/${sample}_matrix_trans_TES.gz
    
    # 绘图
    plotProfile -m ${sample}/${sample}_matrix_trans_TSS.gz \
        -out ${sample}/${sample}_trans_TSS_profile.pdf \
	--regionsLabel "Broad" "Promoter" "TSS" \
	--legendLocation best \
	--samplesLabel "" \
        --plotTitle "TSS Signal - ${sample}"
    
    plotProfile -m ${sample}/${sample}_matrix_trans_TES.gz \
        -out ${sample}/${sample}_trans_TES_profile.pdf \
        --regionsLabel "Broad" "Promoter" "TSS" \
        --legendLocation best \
	--samplesLabel "" \
        --plotTitle "TES Signal - ${sample}"

    plotHeatmap -m ${sample}/${sample}_matrix_trans_TSS.gz \
	-out "${sample}/heatmap_${sample}.trans.TSS.pdf" \
	--colorMap RdBu \
	--whatToShow 'heatmap and colorbar' \
	--regionsLabel "Broad" "Promoter" "TSS" \
	--legendLocation none

    plotHeatmap -m ${sample}/${sample}_matrix_trans_TES.gz \
	-out "${sample}/heatmap_${sample}.trans.TES.pdf" \
	--colorMap RdBu \
	--whatToShow 'heatmap and colorbar' \
	--regionsLabel "Broad" "Promoter" "TSS" \
	--legendLocation none

done
