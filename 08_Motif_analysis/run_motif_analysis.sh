#!/bin/bash
# 文件名：run_all_motif_analysis.sh
# 用途：为每个样本的两种基因列表，提取对应Peak并独立进行Motif分析
# 核心：使用.annotate文件坐标，保证分析区域与注释区域完全一致

echo "=== 开始批量基因-Peak关联Motif分析 ==="
echo "基于讨论结论，将使用.annotate文件坐标进行分析。"

# ========== 1. 核心路径与参数设置 ==========
# 【请根据你的实际路径修改】
ANNOTATE_DIR="../../ChIP_seq/sicer2"
SAMPLE_LIST_FILE="sample.name" # 当前目录下的样本列表
GENOME="mm10"
HOMER_TOOL_PATH="/pathto/software/bin/Miniconda/envs/sicer2-env/bin/findMotifsGenome.pl"

# 输出主目录（脚本所在目录下生成）
MAIN_OUTPUT_DIR="./MotifAnalysis_Results_$(date +%Y%m%d_%H%M%S)"
echo -e "\n所有结果将输出至: $(pwd)/${MAIN_OUTPUT_DIR}/"

# ========== 2. 检查环境和输入文件 ==========
echo -e "\n--- 步骤1: 检查输入与工具 ---"
if [ ! -f "$SAMPLE_LIST_FILE" ]; then
    echo "错误: 未找到样本列表文件 $SAMPLE_LIST_FILE"
    exit 1
fi

if [ ! -f "$HOMER_TOOL_PATH" ]; then
    echo "错误: 未找到HOMER工具 $HOMER_TOOL_PATH"
    exit 1
fi

SAMPLE_COUNT=$(wc -l < "$SAMPLE_LIST_FILE")
echo "找到 $SAMPLE_COUNT 个待处理样本。"

# ========== 3. 主分析循环 ==========
echo -e "\n--- 步骤2: 开始循环处理每个样本 ---"

while IFS= read -r SAMPLE_NAME; do
    if [ -z "$SAMPLE_NAME" ]; then continue; fi
    echo -e "\n>>>>>> 正在处理样本：$SAMPLE_NAME <<<<<<"

    # 3.1 定义关键文件路径
    ANNOTATE_FILE="${ANNOTATE_DIR}/${SAMPLE_NAME}_peaks.annotate"
    SAMPLE_GENE_DIR="./${SAMPLE_NAME}" # 假设基因列表在样本同名子目录
    BROAD_GENE_LIST="${SAMPLE_GENE_DIR}/broad_genes.bed"
    PROMOTER_TSS_GENE_LIST="${SAMPLE_GENE_DIR}/promoter_tss_genes.bed"
    
    # 检查必需文件
    MISSING_FILE=0
    for FILE in "$ANNOTATE_FILE" "$BROAD_GENE_LIST" "$PROMOTER_TSS_GENE_LIST"; do
        if [ ! -f "$FILE" ]; then
            echo "  [警告] 文件不存在: $FILE"
            MISSING_FILE=1
        fi
    done
    if [ $MISSING_FILE -eq 1 ]; then
        echo "  [跳过] 样本 $SAMPLE_NAME 因文件缺失被跳过。"
        continue
    fi

    # 3.2 为当前样本创建专属输出目录
    SAMPLE_OUTPUT_DIR="${MAIN_OUTPUT_DIR}/${SAMPLE_NAME}"
    mkdir -p "${SAMPLE_OUTPUT_DIR}/logs" "${SAMPLE_OUTPUT_DIR}/tmp"

    # 3.3 处理两种基因列表
    for GENE_LIST_TYPE in "broad" "promoter_tss"; do
        echo -e "\n  ** 处理基因列表类型: ${GENE_LIST_TYPE} **"
        GENE_LIST_FILE="${SAMPLE_GENE_DIR}/${GENE_LIST_TYPE}_genes.bed"
        OUTPUT_BED="${SAMPLE_OUTPUT_DIR}/${SAMPLE_NAME}_${GENE_LIST_TYPE}_peaks.bed"
        MOTIF_OUTPUT_DIR="${SAMPLE_OUTPUT_DIR}/${GENE_LIST_TYPE}_motif"

        # A. 从基因列表中提取Ensembl ID (第4列) - 去掉版本号
        cut -f4 "$GENE_LIST_FILE" | awk -F'.' '{print $1}' | sort -u > "${SAMPLE_OUTPUT_DIR}/tmp/${GENE_LIST_TYPE}_gene_ids.txt"
        GENE_NUM=$(wc -l < "${SAMPLE_OUTPUT_DIR}/tmp/${GENE_LIST_TYPE}_gene_ids.txt")
        echo "    从列表中找到 $GENE_NUM 个唯一基因ID (已去除版本号)。"

        # B. 从.annotate文件中匹配基因并获取Peak ID (排除Intergenic)
        # 匹配第15列(Nearest Ensembl)，并过滤掉第8列为'Intergenic'的行
        awk -F'\t' -v gene_id_file="${SAMPLE_OUTPUT_DIR}/tmp/${GENE_LIST_TYPE}_gene_ids.txt" '
            BEGIN {
                while((getline < gene_id_file) > 0) {genes[$1]=1}
                close(gene_id_file)
            }
            NR>1 && ($15 in genes) && ($8 != "Intergenic") {print $1}
        ' "$ANNOTATE_FILE" | sort -u > "${SAMPLE_OUTPUT_DIR}/tmp/${GENE_LIST_TYPE}_peak_ids.txt"

        PEAK_ID_NUM=$(wc -l < "${SAMPLE_OUTPUT_DIR}/tmp/${GENE_LIST_TYPE}_peak_ids.txt")
        echo "    在.annotate中关联到 $PEAK_ID_NUM 个非Intergenic的Peak。"

        # C. 核心：直接从.annotate文件提取Peak坐标 (BED6格式)
        if [ "$PEAK_ID_NUM" -lt 10 ]; then
            echo "    [跳过] Peak数量($PEAK_ID_NUM)少于10，跳过Motif分析。"
            continue
        fi

        echo "    正在提取Peak坐标并生成BED文件..."
        awk -F'\t' -v peak_id_file="${SAMPLE_OUTPUT_DIR}/tmp/${GENE_LIST_TYPE}_peak_ids.txt" '
            BEGIN {
                while((getline < peak_id_file) > 0) {pid[$1]=1}
                close(peak_id_file)
            }
            $1 in pid {
                # 直接从.annotate文件输出BED6格式: chr, start, end, name, score, strand
                # 使用第2,3,4,1,5列，第6列固定为"."
                print $2 "\t" $3 "\t" $4 "\t" $1 "\t" $5 "\t."
            }
        ' "$ANNOTATE_FILE" > "$OUTPUT_BED"

        # D. 运行HOMER Motif分析
        echo "    正在运行findMotifsGenome.pl (此步骤较慢)..."
        mkdir -p "$MOTIF_OUTPUT_DIR"
        
        # 记录开始时间
        START_TIME=$(date +%s)
        
        # 运行motif分析，后台记录日志
        "$HOMER_TOOL_PATH" "$OUTPUT_BED" "$GENOME" "$MOTIF_OUTPUT_DIR" \
            -size 200 \
            -mask \
            -p 4 > "${SAMPLE_OUTPUT_DIR}/logs/${GENE_LIST_TYPE}_motif.log" 2>&1
        
        # 检查运行结果
        if [ -f "${MOTIF_OUTPUT_DIR}/knownResults.txt" ]; then
            END_TIME=$(date +%s)
            DURATION=$((END_TIME - START_TIME))
            TOP_MOTIF=$(head -n 2 "${MOTIF_OUTPUT_DIR}/knownResults.txt" | tail -n 1 | cut -f1)
            echo -e "    [成功] Motif分析完成，耗时 ${DURATION} 秒。"
            echo "    结果目录: $MOTIF_OUTPUT_DIR"
            echo "    最显著Motif: $TOP_MOTIF"
        else
            echo -e "    [失败] Motif分析未成功完成，请检查日志: ${SAMPLE_OUTPUT_DIR}/logs/${GENE_LIST_TYPE}_motif.log"
        fi
    done

    # 3.4 清理当前样本的临时文件
    rm -rf "${SAMPLE_OUTPUT_DIR}/tmp"

done < "$SAMPLE_LIST_FILE"

# ========== 4. 生成最终汇总报告 ==========
echo -e "\n--- 步骤3: 生成分析摘要报告 ---"
SUMMARY_FILE="${MAIN_OUTPUT_DIR}/Analysis_Summary.txt"
{
    echo "基因-Peak关联Motif分析汇总报告"
    echo "==============================="
    echo "分析时间: $(date)"
    echo "样本列表: ${SAMPLE_LIST_FILE}"
    echo "基因组: ${GENOME}"
    echo "总样本数: ${SAMPLE_COUNT}"
    echo "输出主目录: ${MAIN_OUTPUT_DIR}"
    echo ""
    echo "各样本分析状态:"
    echo "-----------------------------------------------------------------------"
} > "$SUMMARY_FILE"

# 遍历所有样本输出目录，统计结果
for SAMPLE_DIR in "${MAIN_OUTPUT_DIR}"/*/; do
    if [ -d "$SAMPLE_DIR" ] && [ "$(basename "$SAMPLE_DIR")" != "logs" ]; then
        SAMPLE_NAME=$(basename "$SAMPLE_DIR")
        echo -n "[${SAMPLE_NAME}] " >> "$SUMMARY_FILE"
        
        for TYPE in "broad" "promoter_tss"; do
            BED_FILE="${SAMPLE_DIR}/${SAMPLE_NAME}_${TYPE}_peaks.bed"
            MOTIF_DIR="${SAMPLE_DIR}/${TYPE}_motif"
            if [ -f "$BED_FILE" ]; then
                PEAK_NUM=$(wc -l < "$BED_FILE" 2>/dev/null || echo "0")
                if [ -f "${MOTIF_DIR}/knownResults.txt" ]; then
                    echo -n "${TYPE}:${PEAK_NUM}peaks(成功) " >> "$SUMMARY_FILE"
                else
                    echo -n "${TYPE}:${PEAK_NUM}peaks(失败) " >> "$SUMMARY_FILE"
                fi
            else
                echo -n "${TYPE}:无 " >> "$SUMMARY_FILE"
            fi
        done
        echo "" >> "$SUMMARY_FILE" # 换行
    fi
done

{
    echo ""
    echo "后续操作指南:"
    echo "1. 所有Motif结果可在各样本的 'broad_motif/' 和 'promoter_tss_motif/' 目录下查看。"
    echo "2. 打开 homerResults.html 文件在浏览器中可交互式查看富集结果。"
    echo "3. 查看 knownResults.txt 文件获取已知转录因子Motif的详细统计信息。"
    echo "4. 具体运行日志在各样本目录下的 'logs/' 子目录中。"
} >> "$SUMMARY_FILE"

echo -e "\n=== 批量分析全部完成！ ==="
echo "结果总目录: $MAIN_OUTPUT_DIR"
echo "详细摘要已保存至: $SUMMARY_FILE"
echo -e "\n提示：可以使用以下命令快速查看所有成功分析的样本："
echo "grep '成功' ${MAIN_OUTPUT_DIR}/Analysis_Summary.txt"
