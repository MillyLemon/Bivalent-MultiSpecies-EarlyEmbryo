# ============================================================
# 完整整合分析代码
# 输出：GO富集表，网络边文件，TF靶基因GO气泡图（跳过错误）
# ============================================================

library(clusterProfiler)
library(org.Mm.eg.db)
library(igraph)       # 用于快速绘制网络
library(ggplot2)
library(dplyr)


# 2. 读入基因列表 --------------------------------------------------
nr4a2  <- readLines("NR4A2_gene_list.txt")
prdm1  <- readLines("PRDM1_gene_list.txt")
zbtb18 <- readLines("ZBTB18_gene_list.txt")
znf135 <- readLines("ZNF135_gene_list.txt")

gene_lists <- list(NR4A2 = nr4a2, PRDM1 = prdm1, ZBTB18 = zbtb18, ZNF135 = znf135)
cat("Gene list sizes:\n")
print(sapply(gene_lists, length))

# 3. 基因ID转换 ----------------------------------------------------
cat("\nConverting gene symbols to Entrez IDs...\n")
entrez_lists <- lapply(gene_lists, function(genes) {
  bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)$ENTREZID
})
cat("ID conversion done.\n")

# 4. GO BP富集 -----------------------------------------------------
go_result_list <- list()
for (tf in names(entrez_lists)) {
  ego <- enrichGO(
    gene          = entrez_lists[[tf]],
    OrgDb         = org.Mm.eg.db,
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2
  )
  if (!is.null(ego) && nrow(ego@result) > 0) {
    ego <- setReadable(ego, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
    # 清理无效基因符号
    ego@result <- ego@result[!is.na(ego@result$geneID) & ego@result$geneID != "", ]
    if (nrow(ego@result) > 0) {
      go_result_list[[tf]] <- ego
      cat(tf, "GO terms:", nrow(ego@result), "\n")
    }
  }
}
saveRDS(go_result_list, "TF_target_GO_results.rds")

# 5. 安全绘制GO气泡图（跳过可能的错误） ----------------------------
for (tf in names(go_result_list)) {
  tryCatch({
    ego <- go_result_list[[tf]]
    # 如果条目过多，仅取前200个最显著的，防止绘图崩溃
    if (nrow(ego@result) > 200) {
      ego@result <- ego@result[order(ego@result$p.adjust), ][1:200, ]
    }
    p <- dotplot(ego, showCategory = 15, title = paste(tf, "target GO BP"))
    ggsave(paste0(tf, "_GO_dotplot.pdf"), p, width = 10, height = 8)
    cat("Saved dotplot for", tf, "\n")
  }, error = function(e) {
    cat("Dotplot failed for", tf, ":", e$message, "\nSaving table instead.\n")
    write.csv(go_result_list[[tf]]@result, paste0(tf, "_GO_table.csv"), row.names = FALSE)
  })
}

# 6. 整合所有TF的GO表 ----------------------------------------------
go_all <- bind_rows(lapply(names(go_result_list), function(tf) {
  df <- go_result_list[[tf]]@result
  df$TF <- tf
  df
}))
write.csv(go_all, "TF_GO_enrichment_all.csv", row.names = FALSE)

# 7. 构建网络 edge 文件（TF -> 基因 or GO） -------------------------
edges <- data.frame(Source = character(), Target = character(), Type = character(), stringsAsFactors = FALSE)

# TF -> 基因
for (tf in names(gene_lists)) {
  edges <- rbind(edges, data.frame(Source = tf, Target = gene_lists[[tf]], Type = "targets"))
}
# TF -> 显著GO (前5个)
for (tf in names(go_result_list)) {
  top_go <- head(go_result_list[[tf]]@result$ID, 5)
  if (length(top_go) > 0) {
    edges <- rbind(edges, data.frame(Source = tf, Target = top_go, Type = "GO"))
  }
}
write.csv(edges, "network_edges.csv", row.names = FALSE)
cat("Network edges saved:", nrow(edges), "rows.\n")

# 8. 快速绘制网络图（可选） ----------------------------------------
tryCatch({
  g <- graph_from_data_frame(edges, directed = TRUE)
  pdf("network_preview.pdf", width = 12, height = 10)
  plot(g, vertex.size = 4, vertex.label.cex = 0.4, layout = layout_with_fr,
       vertex.color = ifelse(grepl("GO:", V(g)$name), "lightblue", "orange"),
       main = "Motif-Target Gene-GO Network")
  dev.off()
  cat("Network plot saved as network_preview.pdf\n")
}, error = function(e) { cat("Network plot failed:", e$message, "\n") })

# 9. 如果仍需要交叉分析矩阵 -----------------------------------------
icm_go <- read.csv("ICM_H3K27me3_promoter_TSS_GO.csv")
eight_go <- read.csv("8_cell_H3K27me3_promoter_TSS_GO.csv")
icm_sig <- icm_go %>% filter(p.adjust < 0.05)
eight_sig <- eight_go %>% filter(p.adjust < 0.05)

calc_overlap <- function(tf_genes, go_genes) sum(tf_genes %in% go_genes)
extract_genes <- function(go_df) {
  gl <- strsplit(go_df$geneID, "/")
  names(gl) <- go_df$ID
  gl
}
icm_genes <- extract_genes(icm_sig)
eight_genes <- extract_genes(eight_sig)

icm_tfs <- gene_lists[c("NR4A2", "PRDM1", "ZBTB18")]
icm_overlap <- sapply(icm_tfs, function(tf_genes) sapply(icm_genes, calc_overlap, tf_genes = tf_genes))
rownames(icm_overlap) <- names(icm_genes)
icm_overlap <- icm_overlap[rowSums(icm_overlap) > 0, , drop = FALSE]
write.csv(as.data.frame(icm_overlap), "ICM_overlap_matrix.csv")

eight_overlap <- sapply(gene_lists["ZNF135"], function(tf_genes) sapply(eight_genes, calc_overlap, tf_genes = tf_genes))
rownames(eight_overlap) <- names(eight_genes)
eight_overlap <- eight_overlap[rowSums(eight_overlap) > 0, , drop = FALSE]
write.csv(as.data.frame(eight_overlap), "8cell_overlap_matrix.csv")

cat("Analysis complete.\n")
