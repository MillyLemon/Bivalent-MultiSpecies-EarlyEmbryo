# GO分析脚本：按基因分类(Promoter/TSS vs Broad)分别分析
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)
library(enrichplot)

# 设置输出目录
output_dir <- "GO_analysis_by_classification_2"
dir.create(output_dir, showWarnings = FALSE)
dir.create(file.path(output_dir, "results"), showWarnings = FALSE)
dir.create(file.path(output_dir, "plots"), showWarnings = FALSE)

# 读取样本列表
samples <- readLines("sample.name")
cat("发现样本数:", length(samples), "\n")

# 定义GO分析函数
perform_classified_go_analysis <- function(sample_name) {
  cat("\n>>> 处理样本:", sample_name, "\n")
  
  # 1. 读取数据文件
  input_file <- file.path(sample_name, "Gene_classification.entrez.txt")
  if (!file.exists(input_file)) {
    cat("  文件不存在，跳过\n")
    return(NULL)
  }
  
  data <- read.table(input_file, header = TRUE, sep = "\t", 
                     stringsAsFactors = FALSE, na.strings = c("", "NA"))
  
  # 2. 提取Promoter+TSS基因 (Promoter和可能的TSS类别)
  promoter_tss_genes <- data[data$Classification %in% c("Promoter", "TSS"), ]
  promoter_tss_entrez <- na.omit(promoter_tss_genes$Entrez_ID)
  promoter_tss_entrez <- unique(promoter_tss_entrez)
  
  # 3. 提取Broad基因
  broad_genes <- data[data$Classification == "Broad", ]
  broad_entrez <- na.omit(broad_genes$Entrez_ID)
  broad_entrez <- unique(broad_entrez)
  
  # 4. 统计信息
  cat("  总基因数:", nrow(data), "\n")
  cat("  Promoter+TSS基因数:", length(promoter_tss_entrez), 
      " (分类为Promoter/TSS:", nrow(promoter_tss_genes), ")\n")
  cat("  Broad基因数:", length(broad_entrez), 
      " (分类为Broad:", nrow(broad_genes), ")\n")
  
  results <- list()
  
  # 5. 对Promoter+TSS基因进行GO分析
  if (length(promoter_tss_entrez) >= 10) {
    cat("  分析Promoter+TSS基因...")
    
    go_promoter <- enrichGO(gene = promoter_tss_entrez,
                           OrgDb = org.Mm.eg.db,
                           keyType = "ENTREZID",
                           ont = "BP",           # 生物过程
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.2,
                           readable = TRUE)
    
    if (!is.null(go_promoter) && nrow(go_promoter) > 0) {
      cat(" 发现", nrow(go_promoter), "个显著条目\n")
      
      # 保存结果
      result_file <- file.path(output_dir, "results", 
                              paste0(sample_name, "_promoter_TSS_GO.csv"))
      write.csv(go_promoter@result, result_file, row.names = FALSE)
      
      # 生成可视化图片
      plot_file <- file.path(output_dir, "plots", 
                            paste0(sample_name, "_promoter_TSS_GO.pdf"))
      p1 <- tryCatch({
        dotplot(go_promoter, showCategory=10, font.size=10, 
                title=paste(sample_name, "Promoter/TSS Genes GO Enrichment"))
      }, error = function(e) {
        barplot(go_promoter, showCategory=10, font.size=10, 
                title=paste(sample_name, "Promoter/TSS Genes GO Enrichment"))
      })
      
      ggsave(plot_file, p1, width = 12, height = 8, dpi = 300)
      
      results$promoter <- list(
        go_object = go_promoter,
        gene_count = length(promoter_tss_entrez),
        significant_terms = nrow(go_promoter),
        result_file = result_file,
        plot_file = plot_file
      )
    } else {
      cat(" 未发现显著富集条目\n")
    }
  } else {
    cat("  Promoter+TSS基因数不足(<10)，跳过分析\n")
  }
  
  # 6. 对Broad基因进行GO分析
  if (length(broad_entrez) >= 10) {
    cat("  分析Broad基因...")
    
    go_broad <- enrichGO(gene = broad_entrez,
                        OrgDb = org.Mm.eg.db,
                        keyType = "ENTREZID",
                        ont = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.2,
                        readable = TRUE)
    
    if (!is.null(go_broad) && nrow(go_broad) > 0) {
      cat(" 发现", nrow(go_broad), "个显著条目\n")
      
      # 保存结果
      result_file <- file.path(output_dir, "results", 
                              paste0(sample_name, "_broad_GO.csv"))
      write.csv(go_broad@result, result_file, row.names = FALSE)
      
      # 生成可视化图片
      plot_file <- file.path(output_dir, "plots", 
                            paste0(sample_name, "_broad_GO.pdf"))
      p2 <- tryCatch({
        dotplot(go_broad, showCategory=10, font.size=10, 
                title=paste(sample_name, "Broad Genes GO Enrichment"))
      }, error = function(e) {
        barplot(go_broad, showCategory=10, font.size=10, 
                title=paste(sample_name, "Broad Genes GO Enrichment"))
      })
      
      ggsave(plot_file, p2, width = 12, height = 8, dpi = 300)
      
      results$broad <- list(
        go_object = go_broad,
        gene_count = length(broad_entrez),
        significant_terms = nrow(go_broad),
        result_file = result_file,
        plot_file = plot_file
      )
    } else {
      cat(" 未发现显著富集条目\n")
    }
  } else {
    cat("  Broad基因数不足(<10)，跳过分析\n")
  }
  
  # 7. 保存基因列表文件（用于后续分析）
  if (length(promoter_tss_entrez) > 0) {
    list_file <- file.path(output_dir, "results", 
                          paste0(sample_name, "_promoter_TSS_genes.txt"))
    write.table(data.frame(Entrez_ID = promoter_tss_entrez),
                list_file, row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
  
  if (length(broad_entrez) > 0) {
    list_file <- file.path(output_dir, "results", 
                          paste0(sample_name, "_broad_genes.txt"))
    write.table(data.frame(Entrez_ID = broad_entrez),
                list_file, row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
  
  return(results)
}

# 7. 批量处理所有样本
all_results <- list()
summary_data <- data.frame()

cat("\n>>> 开始批量GO分析...\n")

for (sample in samples) {
  results <- perform_classified_go_analysis(sample)
  
  if (!is.null(results)) {
    all_results[[sample]] <- results
    
    # 收集统计信息
    promoter_terms <- if (!is.null(results$promoter)) results$promoter$significant_terms else 0
    broad_terms <- if (!is.null(results$broad)) results$broad$significant_terms else 0
    
    summary_data <- rbind(summary_data, data.frame(
      Sample = sample,
      Promoter_TSS_Genes = if (!is.null(results$promoter)) results$promoter$gene_count else 0,
      Broad_Genes = if (!is.null(results$broad)) results$broad$gene_count else 0,
      Promoter_TSS_Sig_Terms = promoter_terms,
      Broad_Sig_Terms = broad_terms,
      stringsAsFactors = FALSE
    ))
  }
}

# 8. 生成总体统计报告
if (nrow(summary_data) > 0) {
  cat("\n>>> 生成总体统计报告...\n")
  
  # 保存详细统计
  summary_file <- file.path(output_dir, "GO_analysis_summary.csv")
  write.csv(summary_data, summary_file, row.names = FALSE)
  cat("  统计报告已保存至:", summary_file, "\n")
  
  # 打印简要统计
  cat("\n>>> 样本分析概况:\n")
  print(summary_data, row.names = FALSE)
  
  # 统计成功分析的样本
  successful_promoter <- sum(summary_data$Promoter_TSS_Sig_Terms > 0)
  successful_broad <- sum(summary_data$Broad_Sig_Terms > 0)
  
  cat("\n>>> 分析成功率:\n")
  cat("  Promoter/TSS基因分析成功样本数:", successful_promoter, "/", nrow(summary_data), "\n")
  cat("  Broad基因分析成功样本数:", successful_broad, "/", nrow(summary_data), "\n")
}

# 9. 生成合并比较图（可选）
if (length(all_results) >= 2) {
  cat("\n>>> 生成样本间比较图...\n")
  
  # 收集所有显著的GO条目
  all_go_terms <- list()
  
  for (sample in names(all_results)) {
    if (!is.null(all_results[[sample]]$promoter) && 
        all_results[[sample]]$promoter$significant_terms > 0) {
      promoter_go <- all_results[[sample]]$promoter$go_object@result
      promoter_go$Sample <- sample
      promoter_go$Gene_Type <- "Promoter_TSS"
      all_go_terms[[paste(sample, "promoter")]] <- head(promoter_go, 10)
    }
    
    if (!is.null(all_results[[sample]]$broad) && 
        all_results[[sample]]$broad$significant_terms > 0) {
      broad_go <- all_results[[sample]]$broad$go_object@result
      broad_go$Sample <- sample
      broad_go$Gene_Type <- "Broad"
      all_go_terms[[paste(sample, "broad")]] <- head(broad_go, 10)
    }
  }
  
  if (length(all_go_terms) > 0) {
    combined_terms <- do.call(rbind, all_go_terms)
    
    # 保存合并数据
    combined_file <- file.path(output_dir, "all_significant_GO_terms.csv")
    write.csv(combined_terms, combined_file, row.names = FALSE)
    cat("  所有显著GO条目已保存至:", combined_file, "\n")
  }
}

cat("\n>>> 分析完成！\n")
cat("输出目录结构:\n")
cat(output_dir, "/\n")
cat("├── results/          # CSV结果文件和基因列表\n")
cat("│   ├── {样本}_promoter_TSS_GO.csv\n")
cat("│   ├── {样本}_broad_GO.csv\n")
cat("│   ├── {样本}_promoter_TSS_genes.txt\n")
cat("│   └── {样本}_broad_genes.txt\n")
cat("├── plots/            # 可视化图片\n")
cat("│   ├── {样本}_promoter_TSS_GO.pdf\n")
cat("│   └── {样本}_broad_GO.pdf\n")
cat("└── GO_analysis_summary.csv  # 总体统计报告\n")
