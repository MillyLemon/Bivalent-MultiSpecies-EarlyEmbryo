# ============================================================
# Dispersion 参数敏感性测试
#
# 测试方法：
#   在 ClassifyGenes.R 中，修改 exactTest 的 dispersion 参数：
#     pvals <- exactTest(d, colnames(region_counts)[2:1], dispersion = 0.05)
#     pvals <- exactTest(d, colnames(region_counts)[2:1], dispersion = 0.10)
#     pvals <- exactTest(d, colnames(region_counts)[2:1], dispersion = 0.20)
#   每次修改后重新运行 ClassifyGenes.R，统计 Gene_classification.txt 中各类基因数量。
#
# 测试结果：
#   详细结果表见 Figures/data/dispersion_test_results.csv
#   三种 dispersion 值下，Broad/Promoter/TSS 基因数量差异极小，分类结果稳定。
#
# 结论：
#   分类结果对 dispersion 参数不敏感，最终选定 dispersion = 0.05。
# ============================================================