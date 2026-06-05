# Figure 3-14: H3K27me3 distribution patterns and expression (4‑cell stage)

This sub‑directory contains the scripts for generating the three component panels of Thesis Figure 3-14.

## Workflow

1. **Gene count bar charts** (`plot_gene_counts.py`)  
   统计每个样本中 Broad、Promoter、TSS 三类基因的数量  
   **输出**: `combined_gene_counts_plot.fix.pdf`

2. **TSS/TES signal profiles** (`plot_tss_tes_profiles.sh`)  
   使用 deepTools 绘制三类基因的 TSS 和 TES 信号分布图  
   **输出**: `*_trans_TES_profile.pdf`（正文选用 TES profile）

3. **Expression boxplots** (`plot_expression_boxplot.py`)  
   绘制三类基因的 TPM 表达箱线图，含 Mann‑Whitney U 检验显著性标注  
   **输出**: `combined_sample_boxplots_log2_plus1_with_zero.fix.pdf`

4. **Assemble final figure** (Adobe Illustrator, manual)  
   将 A、B、C 三个子图拼接为论文图 3-14

## Scripts

| Script | Language | Description |
|--------|----------|-------------|
| `plot_gene_counts.py` | Python | Generate bar charts of gene counts for Broad, Promoter, and TSS categories per sample |
| `plot_tss_tes_profiles.sh` | Shell | Generate TSS/TES signal profiles and heatmaps using deepTools |
| `plot_expression_boxplot.py` | Python | Generate boxplots with Mann‑Whitney U significance tests between categories |

## Notes

- **Final assembly**: Individual panels were combined into Figure 3-14 using **Adobe Illustrator**.
- **Stage selection**: The 4‑cell stage was selected for the main figure; other stages are shown in Appendix Figures A4–A15.
- **Panel B selection**: The script generates both TSS and TES profiles. **TES profiles** were used in the main figure.