# 循环处理当前目录下所有以 .pmd.bam 结尾的文件
for bam_file in *.pmd.bam; do
    echo "正在处理: $bam_file"
    # 定义输出文件名（在原文件名基础上添加 .standard.bam 后缀）
    output_file="${bam_file%.pmd.bam}.standard.bam"
    
    # 这里插入您已验证成功的方案B脚本核心部分，并将输入输出换成变量
    samtools view -h "$bam_file" | \
    awk 'BEGIN {OFS="\t"}
        /^@SQ/ {
            if ($2 ~ /SN:MT\./ || $2 ~ /SN:X\./ || $2 ~ /SN:Y\./ || $2 ~ /SN:[0-9]+\./) {
                split($2, parts, /[:\.]/);
                chrom_id = parts[2];
                if (chrom_id == "MT") { new_name = "chrM"; }
                else if (chrom_id == "X") { new_name = "chrX"; }
                else if (chrom_id == "Y") { new_name = "chrY"; }
                else { new_name = "chr" chrom_id; }
                sub(/SN:[^ \t]+/, "SN:" new_name, $2);
            }
            print $0;
            next;
        }
        !/^@/ {
            if ($3 != "*") {
                split($3, a, ".");
                chrom_id = a[1];
                if (chrom_id == "M" || chrom_id == "MT") { $3 = "chrM"; }
                else if (chrom_id == "X") { $3 = "chrX"; }
                else if (chrom_id == "Y") { $3 = "chrY"; }
                else { $3 = "chr" chrom_id; }
            }
            print $0;
            next;
        }
        { print }' | \
    samtools view -bS -o "$output_file" -
    
    echo "已生成: $output_file"
    echo "-------------------------"
done
