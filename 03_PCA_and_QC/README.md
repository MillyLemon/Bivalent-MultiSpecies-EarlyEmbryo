# 03_PCA_and_QC
PCA‑based quality control for ChIP‑seq samples across five species. This module evaluates replicate reproducibility, identifies outlier samples, and determines filtering thresholds for downstream analysis.
## Workflow overview
| 步骤 | 脚本/工具 | 输出 |
|------|-----------|------|
| 1 | SICER2 | scoreisland files |
| 2 | prepare_k4_matrix.sh / prepare_k27_matrix.sh | sorted BED files |
| 3 | bedtools merge | merged_peak_union.bed |
| 4 | bedtools coverage | coverage_matrix.txt |
| 5 | pca_rpkm.R | RPKM calculation → PCA plots |
| 6 | histogram_rpkm.R | RPKM histograms |


## Directory contents

| File | Description |
|------|-------------|
| `prepare_k4_matrix.sh` | Preprocess H3K4me3 data: scoreisland → union BED → coverage matrix |
| `prepare_k27_matrix.sh` | Preprocess H3K27me3 data: scoreisland → union BED → coverage matrix |
| `run_all_preprocessing.sh` | Master script to run both K4 and K27 preprocessing |
| `pca_rpkm.R` | Main PCA analysis: RPKM calculation, PCA decomposition, and plotting |
| `histogram_rpkm.R` | RPKM distribution histograms for threshold determination |

## Usage

### Step 1: Preprocessing

The preprocessing step converts SICER2 scoreisland files into coverage matrices.

**Option A – Run all at once:**
```bash
bash run_all_preprocessing.sh
```
**Option B – Run separately for each mark:**
```bash
bash prepare_k4_matrix.sh    # For H3K4me3
bash prepare_k27_matrix.sh   # For H3K27me3
```
Input: SICER2 `*_f2q30_pmd-W200-G600.scoreisland` files
Output: `k4.coverage_matrix.txt` and `k27.coverage_matrix.txt`

The coverage matrix has one column per sample and one row per peak region in the union set. Each cell contains the number of reads from that sample overlapping that peak.

### Step 2: PCA analysis
`r
source("pca_rpkm.R")
`
### Step 3: Histogram analysis
`r
source("histogram_rpkm.R")
`

### RPKM calculation
To normalise for both sequencing depth and peak length:
```
RPKM = reads_within_region / ((region_length / 1000) × (total_mapped_reads / 10⁶))
```
This ensures broad domains (e.g. H3K27me3) and sharp peaks (e.g. H3K4me3) are comparable on the same scale.

### Filtering strategy
|Criterion|Threshold|Rationale|
|------|------|------|
|RPKM ≥ 0.5	|Background noise filter|Removes low‑abundance unreliable signals|
|Present in ≥ 2 samples|	Replicate consistency	|Excludes peaks detected in only one sample|
These thresholds were chosen based on RPKM distribution histograms (see `histogram_rpkm.R`), combined with published ChIP‑seq filtering standards.
### Key QC findings
Outlier removed: 4cell_H3K4me3_rep2 in mouse dataset GSE73952 was excluded based on PCA position and IGV visual inspection showing consistently low genome‑wide signal

Single‑replicate stages: Retained after confirming their PCA positions clustered with similar‑stage samples in independent datasets

Cross‑species comparison: The same filtering criteria were applied uniformly to all five species

### Software requirements
Shell: bash, bedtools v2.30.0

R: ≥ 4.0.2 with packages: ggplot2, ggrepel, FactoMineR, factoextra, pheatmap, reshape2, dplyr, tidyr

### Related figure
Figure 3-1: PCA clustering of H3K4me3 and H3K27me3 samples

Appendix Figures A2-A3: RPKM distribution histograms
