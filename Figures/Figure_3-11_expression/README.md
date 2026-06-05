# Figure 3-11: Expression boxplots of three gene categories

This sub‑directory contains the scripts and example data for generating the expression boxplots shown in Thesis Figure 3-11.

## Workflow
1. **Extract gene IDs** (`extract_ensembl_genes.sh`)  
   从 `*_genic.annotate` 文件中提取三类基因的 Ensembl ID 列表  
   **输出**: `*_K27only_ensembl.txt`, `*_K4only_ensembl.txt`, `*_bivalent_ensembl.txt`

2. **Build TPM matrix** (`build_tpm_matrix.R`)  
   从 salmon 输出目录读取所有 `quant.genes.sf` 文件，合并为 TPM 表达矩阵  
   **输出**: `tpm_matrix.csv`

3. **Plot single‑species boxplot** (`plot_expression_boxplot.R`)  
   读取基因列表和 TPM 矩阵，生成单个物种的箱线图  
   **输出**: `*_boxplot.pdf`

4. **Assemble final figure** (Adobe Illustrator, manual)  
   将多个物种的单图拼接为论文图 3-11

## Scripts

| Script | Description |
|--------|-------------|
| `extract_ensembl_genes.sh` | Extract unique Ensembl gene IDs (version numbers removed) from filtered `.annotate` files for K27only, K4only, and bivalent categories |
| `build_tpm_matrix.R` | Read all `quant.genes.sf` files from salmon output and merge into a single TPM matrix |
| `plot_expression_boxplot.R` | Read gene lists and TPM matrix, generate a boxplot for a single species with gene count labels |

## Input

| Data | Source | Description |
|------|--------|-------------|
| `*_genic.annotate` | `04_Bivalent_analysis/` | Filtered peak annotation files (intergenic removed) |
| `quant.genes.sf` | `02_RNA_seq_pipeline/` | Salmon gene‑level quantification output |

## Output

| File | Description |
|------|-------------|
| `*_K27only_ensembl.txt` | Ensembl gene IDs for H3K27me3‑only genes |
| `*_K4only_ensembl.txt` | Ensembl gene IDs for H3K4me3‑only genes |
| `*_bivalent_ensembl.txt` | Ensembl gene IDs for bivalent genes |
| `tpm_matrix.csv` | TPM expression matrix for all samples |
| `*_boxplot.pdf` | Single‑species boxplot (ready for AI assembly) |

## Notes

- **Final figure assembly**: Individual species boxplots were combined into the final Figure 3-11 using **Adobe Illustrator**. No R script is provided for this step.
- **Cross‑species application**: The same scripts were applied to all five species. The `plot_expression_boxplot.R` script should be edited to set the correct species name, stage names, and RNA sample names for each species.
- **Example statistics file**: `mouse1_expression_stats.csv` is provided as an example of the output from `plot_expression_boxplot.R`.
