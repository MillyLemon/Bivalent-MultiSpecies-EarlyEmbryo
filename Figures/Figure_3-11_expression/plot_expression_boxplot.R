library(tidyverse)
library(ggplot2)


species_name <- "Mouse"
output_prefix <- "mouse1"  # 文件前缀

stage_names <- c("2_cell", "4_cell", "8_cell", "morula", "ICM", "TE")
rna_samples <- c("WT_2-cell", "WT_4-cell", "WT_8-cell", "WT_Morula", "WT_ICM", "WT_TE")

# 基因列表文件命名规则
# K27only: 2_cell_K27only_ensembl.txt, 4_cell_K27only_ensembl.txt, ...
# K4only: 2_cell_K4only_ensembl.txt, 4_cell_K4only_ensembl.txt, ...
# bivalent: 2_bivalent_ensembl.txt, 4_bivalent_ensembl.txt, ...（注意去掉了_cell）

# ==================== 读取TPM矩阵 ====================

tpm <- read.csv("tpm_matrix.csv", row.names = 1, check.names = FALSE)

# ==================== 定义函数读取基因列表 ====================

get_expr_for_group <- function(stage, group) {
  # group 可以是 "bivalent", "K27only", "K4only"
  if (group == "bivalent") {
    file <- paste0(gsub("_cell", "", stage), "_bivalent_ensembl.txt")
  } else {
    file <- paste0(stage, "_", group, "_ensembl.txt")
  }
  if (!file.exists(file)) {
    warning(paste("File not found:", file))
    return(NULL)
  }
  genes <- readLines(file)
  genes <- genes[genes != ""]
  if (length(genes) == 0) return(NULL)
  
  # 找到对应的RNA样本
  rna_sample <- rna_samples[match(stage, stage_names)]
  common <- intersect(genes, rownames(tpm))
  
  if (length(common) == 0) return(NULL)
  
  data.frame(
    Gene = common,
    Stage = stage,
    Group = group,
    TPM = tpm[common, rna_sample],
    stringsAsFactors = FALSE
  )
}

# ==================== 收集所有数据 ====================

cat("正在处理", species_name, "...\n")

all_data <- list()
for (st in stage_names) {
  for (grp in c("bivalent", "K27only", "K4only")) {
    df <- get_expr_for_group(st, grp)
    if (!is.null(df)) {
      all_data[[paste(st, grp)]] <- df
      cat("  ", st, grp, ": ", nrow(df), " genes\n", sep = "")
    }
  }
}

# 合并数据
combined_df <- do.call(rbind, all_data)
combined_df$logTPM <- log2(combined_df$TPM + 1)

# 确保阶段顺序正确
combined_df$Stage <- factor(combined_df$Stage, levels = stage_names)
combined_df$Group <- factor(combined_df$Group, levels = c("K27only", "K4only", "bivalent"))

# ==================== 保存统计文件（用于后续拼图）====================

stats_df <- combined_df %>%
  group_by(Stage, Group) %>%
  summarise(
    gene_count = n(),
    mean_TPM = mean(TPM),
    median_TPM = median(TPM),
    mean_logTPM = mean(logTPM),
    median_logTPM = median(logTPM),
    q1_logTPM = quantile(logTPM, 0.25),
    q3_logTPM = quantile(logTPM, 0.75),
    .groups = "drop"
  )

# 添加物种信息
stats_df$species <- species_name
stats_df$species_prefix <- output_prefix

# 保存统计文件
stats_output <- paste0(output_prefix, "_expression_stats.csv")
write.csv(stats_df, stats_output, row.names = FALSE)
cat("\n统计文件已保存至:", stats_output, "\n")

# ==================== 计算用于绘图的基因数量标签 ====================

# 自定义颜色：深红、深蓝、紫
group_colors <- c("K27only"   = "#E41A1C",   # 深红
                  "K4only"    = "#377EB8",   # 深蓝
                  "bivalent"  = "#984EA3")   # 紫

# 计算每组基因数量，并定位到上须顶端
group_counts <- combined_df %>%
  group_by(Stage, Group) %>%
  summarise(
    count = n(),
    q3 = quantile(logTPM, 0.75),
    iqr = IQR(logTPM),
    max_val = max(logTPM),
    # 上须端点 = 不超过 Q3 + 1.5*IQR 的最大值
    upper_whisker = min(max_val, q3 + 1.5 * iqr),
    # 标签放在上须上方 0.15 处（可根据需要微调）
    y_pos = upper_whisker + 0.15,
    .groups = "drop"
  )
group_counts$label <- as.character(group_counts$count)

# ==================== 绘制箱线图 ====================

p <- ggplot(combined_df, aes(x = Stage, y = logTPM, fill = Group)) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.7,
    position = position_dodge(0.8)
  ) +
  geom_text(
    data = group_counts,
    aes(x = Stage, y = y_pos, label = label, group = Group),
    position = position_dodge(0.8),
    size = 5,                    # 数字标签大小
    vjust = 0,
    color = "black"
  ) +
  scale_fill_manual(values = group_colors) +
  theme_minimal() +
  labs(
    title = species_name,
    y = "log2(TPM+1)",
    x = "Developmental Stage"
  ) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 14),  # x轴平放，字号14
    axis.text.y = element_text(size = 14),                         # y轴刻度字号14
    axis.title = element_text(size = 16),                          # 轴标题字号16
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5), # 图标题字号18
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 14),                         # 图例文字字号14
    panel.grid.minor = element_blank()
  )

# ==================== 保存图片（增大画布尺寸） ====================

pdf_output <- paste0(output_prefix, "_boxplot.pdf")
png_output <- paste0(output_prefix, "_boxplot.png")

ggsave(pdf_output, p, width = 12, height = 7.5)  # 原8x5 → 12x7.5
ggsave(png_output, p, width = 12, height = 7.5, dpi = 300)

cat("PDF已保存至:", pdf_output, "\n")
cat("PNG已保存至:", png_output, "\n")

# ==================== 显示图片 ====================

print(p)

# ==================== 打印汇总信息 ====================

cat("\n=== 数据汇总 ===\n")
cat("总基因数:", nrow(combined_df), "\n")
cat("样本数:", length(unique(combined_df$Stage)), "个阶段 × 3个组\n")
cat("\n各阶段基因数量分布:\n")
print(table(combined_df$Stage, combined_df$Group))

cat("\n=== 处理完成 ===\n")
