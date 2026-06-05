library(tidyverse)

# 读取所有quant.genes.sf文件
tpm_files <- list.files("../../../../RNA_seq/salmon_output", 
                        pattern = "quant.genes.sf", 
                        recursive = TRUE, 
                        full.names = TRUE)

tpm_list <- list()

for (file in tpm_files) {
  # 提取样本名
  sample <- basename(dirname(file))
  
  # 读取TPM文件
  df <- read_tsv(file, show_col_types = FALSE)
  
  # 处理基因名（去除版本号）
  df <- df %>%
    mutate(Gene = str_remove(Name, "\\.[0-9]+$")) %>%
    select(Gene, TPM)
  
  # 添加到列表
  tpm_list[[sample]] <- df
}

# 合并为矩阵
tpm_matrix <- tpm_list %>%
  reduce(full_join, by = "Gene") %>%
  column_to_rownames("Gene")

# 设置列名
colnames(tpm_matrix) <- names(tpm_list)

# 保存
write.csv(tpm_matrix, "tpm_matrix.csv")
cat("TPM矩阵已保存: tpm_matrix.csv\n")
