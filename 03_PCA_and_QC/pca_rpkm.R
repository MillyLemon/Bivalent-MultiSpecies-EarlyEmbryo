library(ggplot2)
library(ggrepel)
library(FactoMineR)
library(factoextra)
library(pheatmap)

# 设置过滤阈值
rpkm_threshold <- 0.5  # RPKM过滤阈值

# 定义处理单个mark数据的函数（直接接收文件路径）
process_single_mark_pca <- function(coverage_file, mark_name, color_palette = "Set1") {
  cat("处理", mark_name, "数据...\n")
  
  # 1. 读取覆盖度矩阵
  cat("读取文件:", coverage_file, "\n")
  if (!file.exists(coverage_file)) {
    stop("文件不存在: ", coverage_file)
  }
  coverage_data <- read.table(coverage_file, header = TRUE, sep = "\t",
                              check.names = FALSE, stringsAsFactors = FALSE)
  
  # 检查数据
  cat("矩阵维度:", dim(coverage_data), "\n")
  
  # 2. 数据预处理
  # 提取样本数据（排除peak_id列）
  count_matrix <- coverage_data[, -1]  # 去掉第一列peak_id
  rownames(count_matrix) <- coverage_data$peak_id
  
  # 从peak_id计算每个peak的长度（单位：bp）
  cat("计算peak长度...\n")
  calculate_peak_length <- function(peak_ids) {
    lengths <- numeric(length(peak_ids))
    for(i in seq_along(peak_ids)) {
      parts <- strsplit(peak_ids[i], "_")[[1]]
      if(length(parts) >= 3) {
        start <- as.numeric(parts[2])
        end <- as.numeric(parts[3])
        lengths[i] <- end - start
      } else {
        lengths[i] <- NA
      }
    }
    return(lengths)
  }
  
  peak_lengths <- calculate_peak_length(coverage_data$peak_id)
  cat("peak长度统计(bp):\n")
  print(summary(peak_lengths))
  
  # 检查是否有无效的长度
  if(any(is.na(peak_lengths) | peak_lengths <= 0)) {
    cat("警告：发现", sum(is.na(peak_lengths) | peak_lengths <= 0), "个无效长度，将被过滤\n")
    valid_indices <- !is.na(peak_lengths) & peak_lengths > 0
    count_matrix <- count_matrix[valid_indices, ]
    peak_lengths <- peak_lengths[valid_indices]
    cat("过滤无效长度后保留", nrow(count_matrix), "个peaks\n")
  }
  
  # 3. 数据标准化（使用RPKM）
  cat("数据标准化...\n")
  # 计算每列的总reads数（单位：million）
  total_reads_million <- colSums(count_matrix) / 1e6
  
  # 计算每个peak的长度（单位：kb）
  peak_lengths_kb <- peak_lengths / 1000
  
  # 计算RPKM矩阵
  rpkm_matrix <- matrix(0, nrow = nrow(count_matrix), ncol = ncol(count_matrix))
  rownames(rpkm_matrix) <- rownames(count_matrix)
  colnames(rpkm_matrix) <- colnames(count_matrix)
  
  for(i in 1:ncol(count_matrix)) {
    rpkm_matrix[, i] <- count_matrix[, i] / total_reads_million[i] / peak_lengths_kb
  }
  
  # 4. 基于RPKM值过滤低表达peaks
  cat("基于RPKM值过滤低表达peaks...\n")
  # 过滤条件：保留至少在2个样本中RPKM >= 0.5的peaks
  keep_peaks <- apply(rpkm_matrix, 1, function(x) sum(x >= rpkm_threshold) >= 2)
  
  # 应用过滤
  rpkm_matrix_filtered <- rpkm_matrix[keep_peaks, ]
  count_matrix_filtered <- count_matrix[keep_peaks, ]
  peak_lengths_filtered <- peak_lengths[keep_peaks]
  
  cat("基于RPKM过滤后保留", nrow(rpkm_matrix_filtered), "个peaks (",
      round(sum(keep_peaks) / length(keep_peaks) * 100, 1), "%)\n")
  cat("过滤阈值: RPKM >=", rpkm_threshold, "\n")
  
  # 输出过滤前后RPKM值分布对比
  cat("\nRPKM值分布对比:\n")
  cat("过滤前RPKM统计:\n")
  print(summary(as.vector(rpkm_matrix)))
  cat("过滤后RPKM统计:\n")
  print(summary(as.vector(rpkm_matrix_filtered)))
  
  # 5. 对数转换
  cat("进行对数转换...\n")
  logrpkm_matrix <- log2(rpkm_matrix_filtered + 1)
  
  # 6. PCA分析
  cat("进行PCA分析...\n")
  pca_result <- prcomp(t(logrpkm_matrix), scale. = TRUE, center = TRUE)
  
  # 7. 提取样本信息
  sample_names <- colnames(count_matrix_filtered)
  pca_data <- as.data.frame(pca_result$x)
  pca_data$sample <- sample_names
  
  # 解析样本信息（假设样本命名格式为：stage_replicate，例如 "2cell_1"）
  # 可根据实际命名调整分隔符和提取位置
  pca_data$stage <- sapply(strsplit(sample_names, "_"), function(x) x[1])
  pca_data$replicate <- sapply(strsplit(sample_names, "_"), function(x) x[2])
  
  correct_order <- c("2cell", "4cell", "8cell", "morula", "ICM", "TE")
  pca_data$stage <- factor(pca_data$stage, levels = correct_order)
  
  cat("样本分组信息:\n")
  print(table(pca_data$stage))
  
  # 8. 绘制PCA图
  cat("绘制PCA图...\n")
  
  # 计算方差解释比例
  variance_explained <- round(summary(pca_result)$importance[2, ] * 100, 1)
  
  # PCA图：按阶段着色
  p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = stage)) +
    geom_point(size = 4, alpha = 0.8) +
    geom_text_repel(aes(label = sample), size = 3, max.overlaps = 20) +
    theme_minimal() +
    labs(title = paste0("PCA - ", mark_name, " Samples (RPKM >= ", rpkm_threshold, ")"),
         x = paste0("PC1 (", variance_explained[1], "%)"),
         y = paste0("PC2 (", variance_explained[2], "%)")) +
    scale_color_brewer(palette = color_palette) +
    stat_ellipse(aes(group = stage), level = 0.8, linetype = 2, alpha = 0.5)
  
  print(p)
  
  # 保存图片
  output_file <- paste0("PCA_", mark_name, "_RPKM", rpkm_threshold, ".pdf")
  ggsave(output_file, width = 10, height = 8)
  cat("保存PCA图:", output_file, "\n")
  
  # 9. 样本相关性热图
  cat("绘制样本相关性热图...\n")
  cor_matrix <- cor(logrpkm_matrix)
  
  # 创建样本注释信息
  annotation_df <- data.frame(
    Stage = pca_data$stage,
    row.names = pca_data$sample
  )
  
  heatmap_file <- paste0("sample_correlation_heatmap_", mark_name, "_RPKM", rpkm_threshold, ".pdf")
  pheatmap(cor_matrix,
           annotation_col = annotation_df,
           annotation_row = annotation_df,
           main = paste0(mark_name, " Sample Correlation - RPKM >= ", rpkm_threshold),
           fontsize = 8,
           filename = heatmap_file)
  cat("保存相关性热图:", heatmap_file, "\n")
  
  # 10. 保存分析结果
  cat("保存分析结果...\n")
  
  # 保存PCA坐标
  pca_coord_file <- paste0("pca_coordinates_", mark_name, "_RPKM", rpkm_threshold, ".csv")
  write.csv(pca_data, pca_coord_file, row.names = FALSE)
  cat("保存PCA坐标:", pca_coord_file, "\n")
  
  # 保存RPKM矩阵
  rpkm_file <- paste0("logrpkm_matrix_", mark_name, "_RPKM", rpkm_threshold, ".csv")
  write.csv(logrpkm_matrix, rpkm_file)
  cat("保存RPKM矩阵:", rpkm_file, "\n")
  
  # 返回结果
  return(list(
    pca_result = pca_result,
    pca_data = pca_data,
    logrpkm_matrix = logrpkm_matrix,
    mark_name = mark_name
  ))
}

# ========== 主程序 ==========
# 假设 k4.coverage_matrix.txt 和 k27.coverage_matrix.txt 位于当前工作目录
# 如果需要指定其他路径，请修改下面的文件路径
k4_file <- "k4.coverage_matrix.txt"
k27_file <- "k27.coverage_matrix.txt"

cat("开始分别处理K4和K27数据...\n")

# 处理K4数据
k4_result <- process_single_mark_pca(k4_file, "H3K4me3", "Set1")

# 处理K27数据
k27_result <- process_single_mark_pca(k27_file, "H3K27me3", "Set1")

# 保存完整的分析结果
save(k4_result, k27_result,
     file = paste0("separated_pca_analysis_results_RPKM", rpkm_threshold, ".RData"))

cat("\n分离PCA分析完成！\n")
cat("生成的文件:\n")
cat("- PCA_H3K4me3_RPKM", rpkm_threshold, ".pdf: H3K4me3样本的独立PCA图\n", sep="")
cat("- PCA_H3K27me3_RPKM", rpkm_threshold, ".pdf: H3K27me3样本的独立PCA图\n", sep="")
cat("- sample_correlation_heatmap_H3K4me3_RPKM", rpkm_threshold, ".pdf: H3K4me3样本相关性热图\n", sep="")
cat("- sample_correlation_heatmap_H3K27me3_RPKM", rpkm_threshold, ".pdf: H3K27me3样本相关性热图\n", sep="")
cat("- pca_coordinates_H3K4me3_RPKM", rpkm_threshold, ".csv: H3K4me3样本PCA坐标\n", sep="")
cat("- pca_coordinates_H3K27me3_RPKM", rpkm_threshold, ".csv: H3K27me3样本PCA坐标\n", sep="")
cat("- logrpkm_matrix_H3K4me3_RPKM", rpkm_threshold, ".csv: H3K4me3样本对数RPKM矩阵\n", sep="")
cat("- logrpkm_matrix_H3K27me3_RPKM", rpkm_threshold, ".csv: H3K27me3样本对数RPKM矩阵\n", sep="")
cat("- separated_pca_analysis_results_RPKM", rpkm_threshold, ".RData: 完整的分离PCA分析结果\n", sep="")