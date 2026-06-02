library(ggplot2)
library(cowplot)

# ==================== 数据读取和清理 ====================
read_clean_stat <- function(file) {
  df <- read.table(file, header = TRUE, stringsAsFactors = FALSE)
  df <- df[!is.na(df$sample) & df$sample != "", ]
  df <- df[!is.na(df$GeneCounts), ]
  df$GeneCounts <- as.numeric(df$GeneCounts)
  cat(file, ": 读取", nrow(df), "行数据\n")
  return(df)
}

cat("=== 数据读取情况 ===\n")
peak_m1 <- read_clean_stat("mouse1.gene.stat")
peak_m2 <- read_clean_stat("mouse2.gene.stat")
peak_hs <- read_clean_stat("human.gene.stat")
peak_rat <- read_clean_stat("rat.gene.stat")
peak_bov <- read_clean_stat("bovine.gene.stat")
peak_pig <- read_clean_stat("porcine.gene.stat")  # 修正文件名

# ==================== 自定义三色配色 ====================
my_colors <- c("k27_gene" = "#E41A1C",   # 深红
               "k4_gene"  = "#377EB8",   # 深蓝
               "biv_gene" = "#984EA3")   # 紫
# ======================================================

# ==================== 设置因子顺序 ====================
peak_m1$Group <- factor(peak_m1$Group, levels = c("k27_gene", "k4_gene", "biv_gene"))
peak_m2$Group <- factor(peak_m2$Group, levels = c("k27_gene", "k4_gene", "biv_gene"))
peak_hs$Group <- factor(peak_hs$Group, levels = c("k27_gene", "k4_gene", "biv_gene"))
peak_rat$Group <- factor(peak_rat$Group, levels = c("k27_gene", "k4_gene", "biv_gene"))
peak_bov$Group <- factor(peak_bov$Group, levels = c("k27_gene", "k4_gene", "biv_gene"))
peak_pig$Group <- factor(peak_pig$Group, levels = c("k27_gene", "k4_gene", "biv_gene"))

# 设置sample因子顺序（保持不变）
peak_m1$sample <- factor(peak_m1$sample, levels = c(
  '2_cell_K27Me3', '2_cell_K4Me3', '2_biv',
  '4_cell_K27Me3', '4_cell_K4Me3', '4_biv',
  '8_cell_K27Me3', '8_cell_K4Me3', '8_biv',
  'morula_K27Me3', 'morula_K4Me3', 'morula_biv',
  'ICM_K27Me3', 'ICM_K4Me3', 'ICM_biv',
  'TE_K27Me3', 'TE_K4Me3', 'TE_biv'
))
# ... 其余sample因子设置与原脚本相同（此处省略以节省篇幅，实际运行时请保留全部）

# ==================== 统一大字号主题 ====================
big_theme <- theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    plot.title = element_text(size = 18, face = "bold"),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(),
    legend.position = "none"
  )
# ======================================================

# ==================== 绘图 ====================
pdf("Total.gene.stats.color.pdf", width = 18, height = 9)

# mouse1
p1 <- ggplot(data = peak_m1, aes(x = sample, y = GeneCounts, fill = Group)) + 
  geom_bar(stat = 'identity', position = 'dodge') +
  scale_fill_manual(values = my_colors) +                     # 指定颜色
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +  # 上方扩展10%
  ggtitle("mouse1") + xlab("") + ylab("GeneCounts") +
  scale_x_discrete(labels = c(
    '2C_K27', '2C_K4', '2C_biv',
    '4C_K27', '4C_K4', '4C_biv',
    '8C_K27', '8C_K4', '8C_biv',
    'Mo_K27', 'Mo_K4', 'Mo_biv',
    'ICM_K27', 'ICM_K4', 'ICM_biv',
    'TE_K27', 'TE_K4', 'TE_biv'
  )) +
  big_theme

# mouse2
p2 <- ggplot(data = peak_m2, aes(x = sample, y = GeneCounts, fill = Group)) + 
  geom_bar(stat = 'identity', position = 'dodge') +
  scale_fill_manual(values = my_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  ggtitle("mouse2") + xlab("") + ylab("GeneCounts") +
  scale_x_discrete(labels = c(
    'PN5_K27', 'PN5_K4', 'PN5_biv',
    '2C_early_K27', '2C_early_K4', '2C_early_biv',
    '2C_late_K27', '2C_late_K4', '2C_late_biv',
    '8C_K27', '8C_K4', '8C_biv',
    'ICM_K27', 'ICM_K4', 'ICM_biv'
  )) +
  big_theme

# human
p3 <- ggplot(data = peak_hs, aes(x = sample, y = GeneCounts, fill = Group)) + 
  geom_bar(stat = 'identity', position = 'dodge') +
  scale_fill_manual(values = my_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  ggtitle("human") + xlab("") + ylab("GeneCounts") +
  scale_x_discrete(labels = c(
    '4C_3PN_K27', '4C_3PN_K4', '4C_3PN_biv',
    '8C_2PN_K27', '8C_2PN_K4', '8C_2PN_biv',
    '8C_3PN_K27', '8C_3PN_K4', '8C_3PN_biv',
    'ICM_2PN_K27', 'ICM_2PN_K4', 'ICM_2PN_biv'
  )) +
  big_theme

# rat
p4 <- ggplot(data = peak_rat, aes(x = sample, y = GeneCounts, fill = Group)) + 
  geom_bar(stat = 'identity', position = 'dodge') +
  scale_fill_manual(values = my_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  ggtitle("rat") + xlab("") + ylab("GeneCounts") +
  scale_x_discrete(labels = c(
    'FGO_K27', 'FGO_K4', 'FGO_biv',
    'MII_K27', 'MII_K4', 'MII_biv',
    '1C_K27', '1C_K4', '1C_biv',
    '2C_K27', '2C_K4', '2C_biv',
    '4C_K27', '4C_K4', '4C_biv',
    '8C_K27', '8C_K4', '8C_biv',
    'Bl_K27', 'Bl_K4', 'Bl_biv'
  )) +
  big_theme

# bovine
p5 <- ggplot(data = peak_bov, aes(x = sample, y = GeneCounts, fill = Group)) + 
  geom_bar(stat = 'identity', position = 'dodge') +
  scale_fill_manual(values = my_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  ggtitle("bovine") + xlab("") + ylab("GeneCounts") +
  scale_x_discrete(labels = c(
    'FGO_K27', 'FGO_K4', 'FGO_biv',
    'MII_K27', 'MII_K4', 'MII_biv',
    '4C_K27', '4C_K4', '4C_biv',
    '8C_K27', '8C_K4', '8C_biv',
    '16C_K27', '16C_K4', '16C_biv',
    'Bl_K27', 'Bl_K4', 'Bl_biv'
  )) +
  big_theme

# porcine
p6 <- ggplot(data = peak_pig, aes(x = sample, y = GeneCounts, fill = Group)) + 
  geom_bar(stat = 'identity', position = 'dodge') +
  scale_fill_manual(values = my_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  ggtitle("porcine") + xlab("") + ylab("GeneCounts") +
  scale_x_discrete(labels = c(
    'FGO_K27', 'FGO_K4', 'FGO_biv',
    'MII_K27', 'MII_K4', 'MII_biv',
    '2C_K27', '2C_K4', '2C_biv',
    '4C_K27', '4C_K4', '4C_biv',
    '8C_K27', '8C_K4', '8C_biv',
    'Bl_K27', 'Bl_K4', 'Bl_biv',
    'Mo_K27', 'Mo_K4', 'Mo_biv'
  )) +
  big_theme

# ==================== 图例 ====================
legend_plot <- ggplot(data = peak_pig, aes(x = sample, y = GeneCounts, fill = Group)) + 
  geom_bar(stat = 'identity', position = 'dodge') +
  scale_fill_manual(values = my_colors,
                    name = "Group",
                    labels = c("H3K27me3 monovalent gene", "H3K4me3 monovalent gene", "Bivalent gene")) +
  theme_bw() + 
  theme(legend.position = "bottom",
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.box.margin = margin(5, 5, 5, 5))

legend <- get_legend(legend_plot)

# ==================== 组合 ====================
plots_grid <- plot_grid(p1, p2, p3, p4, p5, p6, 
                        nrow = 2, ncol = 3, 
                        labels = LETTERS[1:6], 
                        label_size = 21)

final_plot <- plot_grid(plots_grid, legend, 
                        ncol = 1, 
                        rel_heights = c(1, 0.1))

print(final_plot)
dev.off()

cat("\n=== 完成 ===\n")
cat("PDF已生成: Total.gene.stats.pdf (宽度18，高度9，各自Y轴，上方扩展10%)\n")
