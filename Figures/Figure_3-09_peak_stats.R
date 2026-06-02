library(ggplot2)
library(cowplot)

# 读取数据
peak_m1 <- read.table("mouse1.filter.stat", head = TRUE)
peak_m2 <- read.table("mouse2.filter.stat", head = TRUE)
peak_hs <- read.table("human.filter.stat", head = TRUE)
peak_rat <- read.table("rat.filter.stat", head = TRUE)
peak_bov <- read.table("bovine.filter.stat", head = TRUE)
peak_pig <- read.table("porcine.filter.stat", head = TRUE)

# ========== 定义简化标签的函数 ==========
simplify_label <- function(x) {
  x <- gsub("K27Me3", "K27", x, fixed = TRUE)
  x <- gsub("K4Me3", "K4", x, fixed = TRUE)
  return(x)
}
# ========================================

# ========== 三色配色（深红、深蓝、紫）==========
my_colors <- c("k27_monovalent" = "#E41A1C",   # 深红
               "k4_monovalent"  = "#377EB8",   # 深蓝
               "bivalent"       = "#984EA3")   # 紫
# ==============================================

# ========== 统一大字号主题（按新要求调整）=========
big_theme <- theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),  # x轴标签 14
    axis.text.y = element_text(size = 14),                         # y轴标签 14（与x轴一致）
    axis.title.y = element_text(size = 16),                        # y轴标题 16
    plot.title = element_text(size = 18, face = "bold"),           # 子图标题 18
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(),
    legend.position = "none"
  )
# =================================================

# 设置 Group 因子水平
peak_m1$Group <- factor(peak_m1$Group, levels = c("k27_monovalent", "k4_monovalent", "bivalent"))
peak_m2$Group <- factor(peak_m2$Group, levels = c("k27_monovalent", "k4_monovalent", "bivalent"))
peak_hs$Group <- factor(peak_hs$Group, levels = c("k27_monovalent", "k4_monovalent", "bivalent"))
peak_rat$Group <- factor(peak_rat$Group, levels = c("k27_monovalent", "k4_monovalent", "bivalent"))
peak_bov$Group <- factor(peak_bov$Group, levels = c("k27_monovalent", "k4_monovalent", "bivalent"))
peak_pig$Group <- factor(peak_pig$Group, levels = c("k27_monovalent", "k4_monovalent", "bivalent"))

# ========== 处理各数据框的 sample 列（简化并保留顺序）==========

# mouse1
orig_levels_m1 <- c('2_cell_K27Me3','2_cell_K4Me3','2_cell_biv','4_cell_K27Me3','4_cell_K4Me3','4_cell_biv',
                    '8_cell_K27Me3','8_cell_K4Me3','8_cell_biv','morula_K27Me3','morula_K4Me3','morula_biv',
                    'ICM_K27Me3','ICM_K4Me3','ICM_biv','TE_K27Me3','TE_K4Me3','TE_biv')
simp_levels_m1 <- simplify_label(orig_levels_m1)
peak_m1$sample <- simplify_label(peak_m1$sample)
peak_m1$sample <- factor(peak_m1$sample, levels = simp_levels_m1)

# mouse2
orig_levels_m2 <- c('PN5zygote_K27Me3','PN5zygote_K4Me3','PN5zygote_biv','2-cellearly_K27Me3','2-cellearly_K4Me3','2-cellearly_biv',
                    '2-celllate_K27Me3','2-celllate_K4Me3','2-celllate_biv','8-cell_K27Me3','8-cell_K4Me3','8-cell_biv',
                    'ICM_K27Me3','ICM_K4Me3','ICM_biv')
simp_levels_m2 <- simplify_label(orig_levels_m2)
peak_m2$sample <- simplify_label(peak_m2$sample)
peak_m2$sample <- factor(peak_m2$sample, levels = simp_levels_m2)

# human
orig_levels_hs <- c('4cell_3PN_K27Me3','4cell_3PN_K4Me3','4cell_3PN_biv','8cell_2PN_K27Me3','8cell_2PN_K4Me3','8cell_2PN_biv',
                    '8cell_3PN_K27Me3','8cell_3PN_K4Me3','8cell_3PN_biv','ICM_2PN_K27Me3','ICM_2PN_K4Me3','ICM_2PN_biv')
simp_levels_hs <- simplify_label(orig_levels_hs)
peak_hs$sample <- simplify_label(peak_hs$sample)
peak_hs$sample <- factor(peak_hs$sample, levels = simp_levels_hs)

# rat
orig_levels_rat <- c('FGO_K27Me3','FGO_K4Me3','FGO_biv','MII_K27Me3','MII_K4Me3','MII_biv','1C_K27Me3','1C_K4Me3','1C_biv',
                     '2C_K27Me3','2C_K4Me3','2C_biv','4C_K27Me3','4C_K4Me3','4C_biv','8C_K27Me3','8C_K4Me3','8C_biv',
                     'Bl_K27Me3','Bl_K4Me3','Bl_biv')
simp_levels_rat <- simplify_label(orig_levels_rat)
peak_rat$sample <- simplify_label(peak_rat$sample)
peak_rat$sample <- factor(peak_rat$sample, levels = simp_levels_rat)

# bovine
orig_levels_bov <- c('FGO_K27Me3','FGO_K4Me3','FGO_biv','MII_K27Me3','MII_K4Me3','MII_biv','4C_K27Me3','4C_K4Me3','4C_biv',
                     '8C_K27Me3','8C_K4Me3','8C_biv','16C_K27Me3','16C_K4Me3','16C_biv','Bl_K27Me3','Bl_K4Me3','Bl_biv')
simp_levels_bov <- simplify_label(orig_levels_bov)
peak_bov$sample <- simplify_label(peak_bov$sample)
peak_bov$sample <- factor(peak_bov$sample, levels = simp_levels_bov)

# pig
orig_levels_pig <- c('FGO_K27Me3','FGO_K4Me3','FGO_biv','MII_K27Me3','MII_K4Me3','MII_biv','2C_K27Me3','2C_K4Me3','2C_biv',
                     '4C_K27Me3','4C_K4Me3','4C_biv','8C_K27Me3','8C_K4Me3','8C_biv',
                     'Mo_K27Me3','Mo_K4Me3','Mo_biv','Bl_K27Me3','Bl_K4Me3','Bl_biv')
simp_levels_pig <- simplify_label(orig_levels_pig)
peak_pig$sample <- simplify_label(peak_pig$sample)
peak_pig$sample <- factor(peak_pig$sample, levels = simp_levels_pig)

# ========== 创建子图（使用新字号设置）==========

p1 <- ggplot(data = peak_m1, mapping = aes(x = sample, y = PeakCounts, fill = Group)) +
  geom_bar(stat = 'identity', position = 'dodge', size = 5) +
  scale_fill_manual(values = my_colors) +
  ggtitle("mouse1") + xlab("") +
  big_theme

p2 <- ggplot(data = peak_m2, mapping = aes(x = sample, y = PeakCounts, fill = Group)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  scale_fill_manual(values = my_colors) +
  ggtitle("mouse2") + xlab("") +
  big_theme +
  theme(axis.text.x = element_text(angle = 52, hjust = 1, size = 14))  # 保持角度，字号14

p3 <- ggplot(data = peak_hs, mapping = aes(x = sample, y = PeakCounts, fill = Group)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  scale_fill_manual(values = my_colors) +
  ggtitle("human") + xlab("") +
  big_theme

p4 <- ggplot(data = peak_rat, mapping = aes(x = sample, y = PeakCounts, fill = Group)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  scale_fill_manual(values = my_colors) +
  ggtitle("rat") + xlab("") +
  big_theme

p5 <- ggplot(data = peak_bov, mapping = aes(x = sample, y = PeakCounts, fill = Group)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  scale_fill_manual(values = my_colors) +
  ggtitle("bovine") + xlab("") +
  big_theme

p6 <- ggplot(data = peak_pig, mapping = aes(x = sample, y = PeakCounts, fill = Group)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  scale_fill_manual(values = my_colors) +
  ggtitle("porcine") + xlab("") +
  big_theme

# ========== 创建用于提取图例的图（图例标签已修改）==========
legend_plot <- ggplot(data = peak_m1, mapping = aes(x = sample, y = PeakCounts, fill = Group)) + 
  geom_bar(stat = 'identity', position = 'dodge') +
  scale_fill_manual(values = my_colors,
		    name = "Group",
                    labels = c("H3K27me3 monovalent",
                               "H3K4me3 monovalent",
                               "Bivalent")) +
  theme_bw() +
  theme(legend.position = "bottom",
        #legend.title = element_blank(),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.box.margin = margin(5, 5, 5, 5))

legend <- cowplot::get_legend(legend_plot)

# ========== 组合子图网格（2行3列）==========
plots_grid <- cowplot::plot_grid(p1, p2, p3, p4, p5, p6, 
                                  nrow = 2, ncol = 3, 
                                  labels = LETTERS[1:6], 
                                  label_size = 21)   # 图序号 21

# ========== 将网格和图例纵向组合（图例在底部）==========
final_plot <- cowplot::plot_grid(plots_grid, legend, 
                                  ncol = 1, 
                                  rel_heights = c(1, 0.1))

# ========== 输出PDF（增大画布尺寸）==========
pdf("k4_k27_bivalent_peakCounts.color.pdf", width = 18, height = 9)
print(final_plot)
dev.off()
