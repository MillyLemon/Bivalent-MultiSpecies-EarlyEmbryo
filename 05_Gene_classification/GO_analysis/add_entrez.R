# 脚本：使用新的映射表为样本文件添加Entrez ID列 (修改版：未匹配基因输出为NA)
cat(">>> 开始批量处理样本文件，添加Entrez ID列...\n")

# 1. 读取新的参考映射文件 ------------------------------------------------------
cat("1. 加载参考映射表...\n")
ref_map <- read.table("gencode_vM23_ensembl_to_entrez_simple.txt", 
                      header = FALSE, 
                      sep = "\t",
                      stringsAsFactors = FALSE,
                      col.names = c("Ensembl_Base", "Entrez_ID"))

# 转换为快速查找的命名向量
ref_dict <- setNames(ref_map$Entrez_ID, ref_map$Ensembl_Base)

cat("   参考映射表大小:", length(ref_dict), "个映射关系\n")
cat("   前3个映射示例:\n")
print(head(ref_map, 3))

# 2. 读取样本列表 --------------------------------------------------------------
samples <- readLines("sample.name")
cat("\n2. 发现样本数:", length(samples), "\n")
cat("   样本列表:", paste(samples, collapse=", "), "\n")

# 3. 定义处理函数 ------------------------------------------------------------
process_sample <- function(sample_name, ref_dict) {
  input_file <- file.path(sample_name, "Gene_classification.txt")
  output_file <- file.path(sample_name, "Gene_classification.entrez.txt")
  
  # 检查文件是否存在
  if (!file.exists(input_file)) {
    cat("   [跳过] 样本", sample_name, "：输入文件不存在\n")
    return(NULL)
  }
  
  cat("   处理样本:", sample_name, "...")
  
  # 读取样本文件
  sample_data <- tryCatch({
    read.table(input_file, 
               header = TRUE, 
               sep = "\t",
               stringsAsFactors = FALSE,
               quote = "", 
               comment.char = "",
               check.names = FALSE)
  }, error = function(e) {
    cat("   读取文件失败:", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(sample_data)) return(NULL)
  
  # 检查列结构
  if (ncol(sample_data) < 4) {
    cat("   文件列数不足，跳过\n")
    return(NULL)
  }
  
  # 获取原始基因ID（第一列）
  original_ids <- sample_data[[1]]
  
  # 提取基础Ensembl ID（去除版本号）
  base_ids <- sub("\\.\\d+$", "", original_ids)
  
  # 在参考字典中查找匹配的Entrez ID
  # 关键修改：这里的结果已经包含NA
  entrez_ids <- ref_dict[base_ids]
  names(entrez_ids) <- NULL  # 移除名称属性
  
  # 统计匹配结果
  matched_count <- sum(!is.na(entrez_ids))
  total_count <- length(entrez_ids)
  match_rate <- round(matched_count / total_count * 100, 2)
  
  cat(" 匹配", matched_count, "/", total_count, 
      "(", match_rate, "%)\n")
  
  # 构建新的数据框
  # 获取原始列名
  original_colnames <- colnames(sample_data)
  
  # 创建新数据框：第一列原始ID，第二列Entrez ID（包含NA），然后是其他列
  # 关键修改：确保未匹配的值在数据框中就是NA
  new_data <- data.frame(
    gene_id = original_ids,           # 保持原始ID（带版本号）
    Entrez_ID = entrez_ids,           # 新增的Entrez ID列（包含NA）
    stringsAsFactors = FALSE
  )
  
  # 添加原始的其他列（从第二列开始）
  for (i in 2:ncol(sample_data)) {
    new_data[[original_colnames[i]]] <- sample_data[[i]]
  }
  
  # 保存结果文件
  # 关键修改：移除 `na = ""` 参数，让NA值在文件中自然显示为NA
  write.table(new_data, 
              file = output_file,
              sep = "\t",
              row.names = FALSE,
              quote = FALSE)
              # 注意：这里不再设置 `na = ""`，默认NA就会输出为"NA"字符串
  
  return(list(
    sample = sample_name,
    total = total_count,
    matched = matched_count,
    rate = match_rate,
    unmatched_count = total_count - matched_count
  ))
}

# 4. 批量处理所有样本 ---------------------------------------------------------
cat("\n3. 开始批量处理所有样本...\n\n")

results <- list()
failed_samples <- c()

for (sample in samples) {
  result <- process_sample(sample, ref_dict)
  if (!is.null(result)) {
    results[[sample]] <- result
  } else {
    failed_samples <- c(failed_samples, sample)
  }
}

# 5. 生成详细的统计报告 -------------------------------------------------------
cat("\n4. 生成处理统计报告...\n")

if (length(results) > 0) {
  # 创建汇总数据框
  summary_data <- data.frame(
    Sample = sapply(results, function(x) x$sample),
    Total_Genes = sapply(results, function(x) x$total),
    Matched = sapply(results, function(x) x$matched),
    Match_Rate = sapply(results, function(x) x$rate),
    Unmatched = sapply(results, function(x) x$unmatched_count),
    stringsAsFactors = FALSE
  )
  
  # 按匹配率排序
  summary_data <- summary_data[order(summary_data$Match_Rate, decreasing = TRUE), ]
  
  # 打印汇总表
  cat("\n>>> 各样本匹配统计:\n")
  print(summary_data, row.names = FALSE)
  
  # 计算总体统计
  overall_total <- sum(summary_data$Total_Genes)
  overall_matched <- sum(summary_data$Matched)
  overall_rate <- round(overall_matched / overall_total * 100, 2)
  
  cat("\n>>> 总体统计:\n")
  cat("   处理样本数:", nrow(summary_data), "\n")
  cat("   总基因数:", overall_total, "\n")
  cat("   成功匹配基因数:", overall_matched, "\n")
  cat("   总体匹配率:", overall_rate, "%\n")
  
  # 分析匹配率分布
  cat("\n>>> 匹配率分布:\n")
  cat("   最低匹配率:", min(summary_data$Match_Rate), "%\n")
  cat("   最高匹配率:", max(summary_data$Match_Rate), "%\n")
  cat("   平均匹配率:", round(mean(summary_data$Match_Rate), 2), "%\n")
  cat("   中位数匹配率:", median(summary_data$Match_Rate), "%\n")
  
  # 保存详细统计结果
  write.csv(summary_data, 
            "sample_entrez_matching_summary.csv", 
            row.names = FALSE)
  cat("\n   详细统计已保存至: sample_entrez_matching_summary.csv\n")
  
  # 对于匹配率较低的样本，提供诊断建议
  low_rate_samples <- summary_data[summary_data$Match_Rate < 80, ]
  if (nrow(low_rate_samples) > 0) {
    cat("\n>>> 注意：以下样本匹配率较低(<80%):\n")
    print(low_rate_samples, row.names = FALSE)
    cat("\n   建议检查这些样本的基因ID是否与参考文件格式一致。\n")
  }
  
} else {
  cat("   错误：没有成功处理任何样本\n")
}

# 报告失败的样本
if (length(failed_samples) > 0) {
  cat("\n>>> 处理失败的样本:\n")
  cat("   ", paste(failed_samples, collapse=", "), "\n")
  cat("   请检查这些样本目录是否存在或文件格式是否正确。\n")
}

# 6. 生成未匹配基因的汇总报告（可选）------------------------------------------
cat("\n5. 生成未匹配基因分析...\n")

if (length(results) > 0 && any(sapply(results, function(x) x$unmatched_count > 0))) {
  # 收集所有样本的未匹配基因
  all_unmatched <- list()
  
  for (sample in names(results)) {
    input_file <- file.path(sample, "Gene_classification.txt")
    sample_data <- read.table(input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    original_ids <- sample_data[[1]]
    base_ids <- sub("\\.\\d+$", "", original_ids)
    
    # 找出未匹配的基因
    unmatched_mask <- is.na(ref_dict[base_ids])
    if (any(unmatched_mask)) {
      all_unmatched[[sample]] <- data.frame(
        Sample = sample,
        Original_ID = original_ids[unmatched_mask],
        Base_ID = base_ids[unmatched_mask],
        stringsAsFactors = FALSE
      )
    }
  }
  
  if (length(all_unmatched) > 0) {
    unmatched_df <- do.call(rbind, all_unmatched)
    
    # 按Base_ID统计哪些基因在多个样本中都未匹配
    common_unmatched <- table(unmatched_df$Base_ID)
    common_unmatched <- common_unmatched[order(common_unmatched, decreasing = TRUE)]
    
    cat("   发现", length(common_unmatched), "个唯一未匹配基因ID\n")
    
    if (length(common_unmatched) > 0) {
      # 保存未匹配基因列表
      write.table(data.frame(Ensembl_Base = names(common_unmatched),
                            Frequency = as.numeric(common_unmatched)),
                  "common_unmatched_genes.txt",
                  sep = "\t", row.names = FALSE, quote = FALSE)
      cat("   常见未匹配基因列表已保存至: common_unmatched_genes.txt\n")
      
      # 显示最常见的10个未匹配基因
      cat("\n   最常见的10个未匹配基因:\n")
      print(head(data.frame(Ensembl_ID = names(common_unmatched),
                           Frequency = as.numeric(common_unmatched)), 10))
    }
  }
}

cat("\n>>> 处理完成！\n")
cat(">>> 每个样本目录下已生成: Gene_classification.entrez.txt\n")
cat(">>> 文件格式：第一列原始基因ID，第二列Entrez_ID（未匹配基因为NA），然后是原始其他列\n")
cat(">>> 后续GO分析请注意使用 na.omit() 过滤NA值。\n")
