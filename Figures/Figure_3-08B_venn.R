# ============================================================
# Figure 3-8B: Venn diagram of H3K4me3 and H3K27me3 peak overlaps
#
# 数据来源：
#   本脚本使用的数值来自手动统计文件 mouse1.filter.stat。
#   该文件中的每一行数据（如 2_cell_K27Me3: 12651）来源于对
#   过滤 intergenic 后的 .annotate 文件执行 wc -l 命令：
#     wc -l 2_cell_K27only_genic.annotate
#     wc -l 2_cell_K4only_genic.annotate
#     wc -l 2_cell_biv_genic.annotate
#   这些 *_genic.annotate 文件由 04_Bivalent_analysis/filter_all_intergenic.sh
#   从 HOMER annotatePeaks.pl 的输出中过滤生成。
#
# 输入文件: mouse1.filter.stat (手动整理)
# 输出文件: Figure_3-8B_venn.pdf
# ============================================================

# ==================== 加载必要的包 ====================
if (!require("VennDiagram")) install.packages("VennDiagram")
if (!require("gridExtra")) install.packages("gridExtra")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("dplyr")) install.packages("dplyr")
if (!require("tidyr")) install.packages("tidyr")

library(VennDiagram)
library(gridExtra)
library(ggplot2)
library(dplyr)
library(tidyr)

# ==================== 数据设置 ====================
stages <- c("2_cell", "4_cell", "8_cell", "morula", "ICM", "TE")

# 使用提供的实际数据
##all_counts <- list(
#  "2_cell" = list(k27_only = 50459, k4_only = 76243, bivalent = 8299),
#  "4_cell" = list(k27_only = 59824, k4_only = 70733, bivalent = 3070),
#  "8_cell" = list(k27_only = 63891, k4_only = 58077, bivalent = 3480),
#  "morula" = list(k27_only = 46591, k4_only = 66876, bivalent = 4538),
#  "ICM" = list(k27_only = 76480, k4_only = 52071, bivalent = 4727),
#  "TE" = list(k27_only = 72831, k4_only = 53109, bivalent = 4352)
#)
#filter
all_counts <- list(
  "2_cell" = list(k27_only = 12651, k4_only = 37866, bivalent = 3627),
  "4_cell" = list(k27_only = 16075, k4_only = 41164, bivalent = 1339),
  "8_cell" = list(k27_only = 16660, k4_only = 34317, bivalent = 1655),
  "morula" = list(k27_only = 11674, k4_only = 38249, bivalent = 1837),
  "ICM" = list(k27_only = 21828, k4_only = 30347, bivalent = 2519),
  "TE" = list(k27_only = 20853, k4_only = 31080, bivalent = 2522)
)
# ==================== 创建数据汇总表 ====================
# 创建数据框
plot_data <- list()
for (stage in stages) {
  counts <- all_counts[[stage]]
  plot_data[[stage]] <- data.frame(
    Stage = stage,
    H3K27me3_only = counts$k27_only,
    H3K4me3_only = counts$k4_only,
    Bivalent = counts$bivalent,
    Total = counts$k27_only + counts$k4_only + counts$bivalent
  )
}

# 合并数据
summary_df <- do.call(rbind, plot_data)

# 计算百分比
summary_df <- summary_df %>%
  mutate(
    H3K27me3_only_pct = round(H3K27me3_only / Total * 100, 2),
    H3K4me3_only_pct = round(H3K4me3_only / Total * 100, 2),
    Bivalent_pct = round(Bivalent / Total * 100, 2),
    H3K27me3_total = H3K27me3_only + Bivalent,
    H3K4me3_total = H3K4me3_only + Bivalent
  )

# ==================== 1. 保存CSV文件 ====================
write.csv(summary_df, "peak_summary.csv", row.names = FALSE)
cat("✅ 已保存数据文件: peak_summary.csv\n")

# ==================== 2. Venn 图网格（PDF） ====================
pdf("venn_diagrams_grid.pdf", width = 12, height = 8)
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 3)))

for (i in 1:length(stages)) {
  row <- (i-1) %/% 3 + 1
  col <- (i-1) %% 3 + 1
  
  pushViewport(viewport(layout.pos.row = row, layout.pos.col = col))
  
  counts <- all_counts[[stages[i]]]
  k27_total <- counts$k27_only + counts$bivalent
  k4_total <- counts$k4_only + counts$bivalent
  
  draw.pairwise.venn(
    area1 = k27_total,
    area2 = k4_total,
    cross.area = counts$bivalent,
    category = c("H3K27me3", "H3K4me3"),
    fill = c("#E41A1C", "#377EB8"),
    alpha = 0.6,
    cex = 1.5,
    col =NA,
    cat.cex = 1.2,
    cat.pos = c(0, 0),
    cat.dist = 0.05,
    margin = 0.1
  )
  
  # 添加标题
  grid.text(stages[i], x = 0.5, y = 0.95, 
            gp = gpar(fontsize = 12, fontface = "bold"))
  
  # 添加总计
  grid.text(paste("Total:", sum(unlist(counts))), 
            x = 0.5, y = 0.05,
            gp = gpar(fontsize = 10))
  
  popViewport()
}
dev.off()
cat("✅ 已生成图表文件: venn_diagrams_grid.pdf\n")

# ==================== 3. 百分比堆叠柱状图（PDF） ====================
# 准备长格式数据
summary_df_long <- summary_df %>%
  select(Stage, H3K27me3_only, H3K4me3_only, Bivalent, Total) %>%
  mutate(
    H3K27me3_only_pct = H3K27me3_only / Total * 100,
    H3K4me3_only_pct = H3K4me3_only / Total * 100,
    Bivalent_pct = Bivalent / Total * 100
  ) %>%
  select(Stage, ends_with("_pct")) %>%
  pivot_longer(
    cols = ends_with("_pct"),
    names_to = "Peak_Type",
    values_to = "Percentage"
  ) %>%
  mutate(
    Peak_Type = gsub("_pct", "", Peak_Type),
    Peak_Type = factor(Peak_Type,
                       levels = c("H3K27me3_only", "H3K4me3_only", "Bivalent"),
                       labels = c("H3K27me3-only", "H3K4me3-only", "Bivalent"))
  )

# 创建百分比堆叠柱状图
summary_plot_stacked <- ggplot(summary_df_long, 
                                aes(x = Stage, y = Percentage, fill = Peak_Type)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  scale_fill_manual(values = c("H3K27me3-only" = "#ff9999", 
                               "H3K4me3-only" = "#66b3ff", 
                               "Bivalent" = "#99ff99")) +
  labs(title = "Percentage of H3K27me3, H3K4me3 and Bivalent Peaks\nAcross Developmental Stages",
       y = "Percentage (%)",
       fill = "Peak Type",
       x = "Developmental Stage") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 12),
    legend.position = "right",
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  ) +
  # 添加百分比标签（数值较大时显示）
  geom_text(aes(label = ifelse(Percentage > 5, sprintf("%.1f%%", Percentage), "")),
            position = position_stack(vjust = 0.5),
            size = 3.5, color = "black", fontface = "bold")

# 保存为PDF
ggsave("peak_distribution_stacked_percent.pdf", 
       summary_plot_stacked, 
       width = 10, 
       height = 7,
       device = "pdf")
cat("已生成图表文件: peak_distribution_stacked_percent.pdf\n")

# ==================== 输出总结 ====================
cat("\n")
cat(rep("=", 50), sep = "")
cat("\n")
cat("🎉 分析完成！共生成3个文件：\n")
cat(rep("-", 50), sep = "")
cat("\n")
cat("1. peak_summary.csv              - 数据汇总表格\n")
cat("2. venn_diagrams_grid.pdf        - 韦恩图网格（6个阶段）\n")
cat("3. peak_distribution_stacked_percent.pdf - 百分比堆叠柱状图\n")
cat(rep("=", 50), sep = "")
cat("\n\n")

# 显示数据预览
cat("数据预览（前6行）：\n")
print(head(summary_df))

cat("\n总计峰值数：", sum(summary_df$Total), "\n")
cat("各阶段峰值分布：\n")
for (i in 1:nrow(summary_df)) {
  cat(sprintf("  %-10s: %d peaks\n", summary_df$Stage[i], summary_df$Total[i]))
}

cat("\n✅ 所有文件已成功生成！\n")
