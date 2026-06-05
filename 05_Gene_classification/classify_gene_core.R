# ============================================================
# H3K27me3 基因分类核心函数
# 基于 Young et al. (2011) Nucleic Acids Research
# 适配 R 4.x 环境
# 依赖：
#   - GenomicRanges (GRanges, coverage, strand, seqnames 等)
#   - caTools (runmean 滑动窗口平滑)
#   - getRegion() 函数，定义在 functions.R 中
#   - annot 为 GRanges 对象，包含 gene_id 列
#   - cover 为 coverage() 返回的 RleList 对象
#
# 使用前请先：
#   library(GenomicRanges)
#   library(caTools)
#   source("functions.R")
# ============================================================

classify_gene <- function(gene, annot, cover, 
                          peak_window = 200, 
                          min_gene_length = 5000,
                          Promoter_region = c(-3000, -100),
                          TSS_region = c(-99, 1000),
                          Broad_region = c(1001, NA),
                          Promoter_ratio = 1.25,
                          TSS_ratio = 1.25,
                          Broad_ratio = 0.35) {
  
  # 获取基因注释信息
  gene_info <- annot[match(gene, values(annot)$gene_id)]
  
  # 过滤过短基因
  if (width(gene_info) < min_gene_length) {
    return(list(classification = "Unclassified"))
  }
  
  # 动态确定 Broad 区域终止位置（至基因末端）
  if (is.na(Broad_region[2])) {
    Broad_region[2] <- width(gene_info)
  }
  
  # 确定数据提取区间（包含峰值检测窗口）
  data_region <- c(min(Promoter_region, TSS_region, Broad_region) - peak_window,
                   max(Promoter_region, TSS_region, Broad_region) + peak_window)
  
  # 根据链向确定基因组坐标提取区间
  forward <- as.character(strand(gene_info)) == "+"
  if (forward) {
    fetch_region <- start(gene_info) + data_region
  } else {
    fetch_region <- rev(end(gene_info) - data_region)
  }
  
  # 提取覆盖度数据
  data <- getRegion(cover, chr = as.character(seqnames(gene_info)),
                    start = fetch_region[1], end = fetch_region[2])
  if (!forward) data <- rev(data)
  
  # 滑动窗口平滑用于峰检测
  peak_data <- runmean(data, k = peak_window, endrule = "constant")
  
  # 提取指定区域的信号数据
  get_region_signal <- function(region_bounds) {
    w <- (region_bounds[1]:region_bounds[2]) - data_region[1] + 1
    list(data = data[w], peak = peak_data[w])
  }
  
  Promoter <- get_region_signal(Promoter_region)
  TSS      <- get_region_signal(TSS_region)
  Broad    <- get_region_signal(Broad_region)
  
  # 计算各区域平均覆盖度和最大峰值
  mean_cover <- c(mean(Promoter$data), mean(TSS$data), mean(Broad$data))
  max_peak   <- c(max(Promoter$peak), max(TSS$peak), max(Broad$peak))
  
  # 判定分类
  max_region <- which.max(mean_cover)
  classification <- "Unclassified"
  
  if (length(max_region) == 1) {
    if (max_region == 1 && max_peak[1] / max(max_peak[-1]) > Promoter_ratio) {
      classification <- "Promoter"
    } else if (max_region == 2 && max_peak[2] / max(max_peak[-2]) > TSS_ratio) {
      classification <- "TSS"
    } else if (max_region == 3) {
      above_mean <- sum(Broad$peak > mean(Broad$peak)) / length(Broad$peak)
      if (above_mean > Broad_ratio) {
        classification <- "Broad"
      }
    }
  }
  
  return(list(classification = classification))
}