#!/bin/bash
# 保存为 extract_all_genes_ensembl.sh

echo "开始提取基因..."
echo "文件\t类型\t基因数" > gene_counts.tsv

for file in *genic.annotate; do
    [ ! -f "$file" ] && continue
    
    # 判断类型和阶段
    if [[ $file == *given* ]]; then
        stage=$(echo $file | cut -d'_' -f1)
        type="bivalent"
    elif [[ $file == *K27only* ]]; then
        stage=$(echo $file | sed 's/_K27only.*//')
        type="K27only"
    elif [[ $file == *K4only* ]]; then
        stage=$(echo $file | sed 's/_K4only.*//')
        type="K4only"
    else
        continue
    fi
    
    # 提取基因
    output="${stage}_${type}_ensembl.txt"
    awk -F'\t' 'NR>1 {
        ensembl = $15
        gsub(/\.[0-9]+$/, "", ensembl)
        if (ensembl != "") print ensembl
    }' "$file" | sort -u > "$output"
    
    count=$(wc -l < "$output")
    echo -e "$file\t$type\t$count" >> gene_counts.tsv
    echo "生成: $output ($count 个基因)"
done

echo "完成！基因数量统计见 gene_counts.tsv"
