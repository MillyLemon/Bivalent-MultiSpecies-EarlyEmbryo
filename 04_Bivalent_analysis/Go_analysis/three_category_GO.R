# ============================================================
# 三类修饰基因 GO-BP 富集分析（通用版）
# 用于论文图 3-12 和附录图 A16-A19
#
# 输入: *_K4only_ensembl.txt, *_K27only_ensembl.txt, *_bivalent_ensembl.txt
#       每行一个 Ensembl gene ID（无版本号）
# 输出: GO_<stage>_<species>.pdf（三类基因的 GO 气泡图纵向排列）
# ============================================================

library(clusterProfiler)
library(ggplot2)

# ---------- 参数设置（根据物种和阶段修改）----------
SPECIES       <- "Human"          # 物种名（用于标题）
STAGE         <- "8-cell 2PN"     # 发育阶段（用于标题和文件名）
ORG_DB        <- org.Hs.eg.db     # 物种注释包：人=org.Hs.eg.db, 小鼠=org.Mm.eg.db
KEY_TYPE      <- "ENSEMBL"        # 基因 ID 类型

FILE_K4  <- "8cell2PN_K4only_ensembl.txt"
FILE_K27 <- "8cell2PN_K27only_ensembl.txt"
FILE_BIV <- "8cell2PN_bivalent_ensembl.txt"

# ---------- GO 富集函数 ----------
run_go_single <- function(gene_file, category_name, color_label) {
  genes <- readLines(gene_file)
  go <- enrichGO(
    gene          = genes,
    OrgDb         = ORG_DB,
    keyType       = KEY_TYPE,
    ont           = "BP",
    pAdjustMethod = "BH",
    qvalueCutoff  = 0.05
  )

  if (!is.null(go) && nrow(go@result) > 0) {
    p <- dotplot(go, showCategory = 10) +
      ggtitle(paste(SPECIES, STAGE, "-", category_name))
    return(list(go = go, plot = p, n_terms = nrow(go@result)))
  } else {
    message("  No significant terms for ", category_name)
    return(NULL)
  }
}

# ---------- 主程序 ----------
cat("Processing", SPECIES, STAGE, "...\n")

res_k4  <- run_go_single(FILE_K4,  "H3K4me3 monovalent")
res_k27 <- run_go_single(FILE_K27, "H3K27me3 monovalent")
res_biv <- run_go_single(FILE_BIV, "Bivalent")

# 纵向排列三类基因的 GO 图
plots <- list()
if (!is.null(res_k4))  plots$k4  <- res_k4$plot
if (!is.null(res_k27)) plots$k27 <- res_k27$plot
if (!is.null(res_biv)) plots$biv <- res_biv$plot

if (length(plots) > 0) {
  combined <- wrap_plots(plots, ncol = 1)
  outfile <- paste0("GO_", gsub(" ", "_", STAGE), "_", SPECIES, ".pdf")
  ggsave(outfile, combined, width = 10, height = 4 * length(plots))
  cat("Saved:", outfile, "\n")
}

cat("Done.\n")