#!/bin/bash
# 保存为 run_mergePeaks_sicer.sh

echo "Converting SICER island files to broadPeak format..."

# 处理每个发育阶段
stages=("2_cell" "4_cell" "8_cell" "morula" "ICM" "TE")

for stage in "${stages[@]}"; do
    k27_input="merged_${stage}_H3K27me3-W200-G600-FDR0.01-island.bed"
    k4_input="merged_${stage}_H3K4me3-W200-G600-FDR0.01-island.bed"
    
    k27_output="${stage}_H3K27me3_peaks.broadPeak"
    k4_output="${stage}_H3K4me3_peaks.broadPeak"
    
    if [[ -f "$k27_input" && -f "$k4_input" ]]; then
        echo "========================================"
        echo "Processing $stage..."
        echo "========================================"
        
        # 转换H3K27me3文件 (SICER island -> broadPeak)
        echo "Converting $k27_input to $k27_output..."
        # SICER格式: chr start end read_count
        # broadPeak格式: chr start end name score strand signalValue pValue qValue
        awk -v stage="$stage" -v mod="H3K27me3" '
        BEGIN {OFS="\t"; count=1}
        {
            # 创建唯一的peak名称
            peak_name = stage "_" mod "_peak_" count
            # 使用read_count作为score（缩放一下，使其在0-1000范围内）
            score = ($4 > 1000 ? 1000 : $4)
            # 输出broadPeak格式
            print $1, $2, $3, peak_name, score, ".", $4, -1, -1
            count++
        }' "$k27_input" > "$k27_output"
        
        echo "  K27: $(wc -l < "$k27_output") peaks created"
        
        # 转换H3K4me3文件
        echo "Converting $k4_input to $k4_output..."
        awk -v stage="$stage" -v mod="H3K4me3" '
        BEGIN {OFS="\t"; count=1}
        {
            peak_name = stage "_" mod "_peak_" count
            score = ($4 > 1000 ? 1000 : $4)
            print $1, $2, $3, peak_name, score, ".", $4, -1, -1
            count++
        }' "$k4_input" > "$k4_output"
        
        echo "  K4: $(wc -l < "$k4_output") peaks created"
        
        # 运行mergePeaks
        echo "Running mergePeaks..."
        mergePeaks -d given \
           -venn ${stage}_given.venn \
           -prefix ${stage}_given \
           -matrix ${stage}_given \
           "$k27_output" \
           "$k4_output"
        
        # 检查输出
        echo "Checking output files..."
        for ext in "_given.txt" "_given.venn.txt"; do
            if [[ -f "${stage}${ext}" ]]; then
                echo "  Created: ${stage}${ext}"
            fi
        done
        
        echo "Completed $stage"
        echo ""
    else
        echo "Warning: Missing files for $stage"
        echo ""
    fi
done

echo "All done!"
