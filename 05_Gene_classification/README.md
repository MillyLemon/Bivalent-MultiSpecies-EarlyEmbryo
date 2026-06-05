# 05_Gene_classification

Classification of H3K27me3‑modified genes into **Promoter**, **TSS**, and **Broad** enrichment patterns based on ChIP‑seq coverage profiles.

This method is adapted from **Young et al. (2011) *Nucleic Acids Research***. It has been systematically repaired and validated for compatibility with R ≥ 4.0.2 and SICER2 peak calling output.

## Overview

Each gene is divided into three regions:

| Region   | Definition                                             | Typical signal pattern            |
|----------|--------------------------------------------------------|-----------------------------------|
| Promoter | TSS upstream 3000 bp to TSS upstream 100 bp            | Narrow peak upstream of TSS       |
| TSS      | TSS upstream 100 bp to TSS downstream 1000 bp          | Narrow peak at TSS                |
| Broad    | TSS downstream 1001 bp to gene termination site        | Broad enrichment across gene body |

A 200 bp sliding window is used to smooth the coverage signal. Each gene is assigned to the region with the highest mean coverage, subject to the following criteria:

1. **Promoter**: Peak in the promoter region must be at least **25 % higher** than the maximum peak in any other region.
2. **TSS**: Peak in the TSS region must be at least **25 % higher** than the maximum peak in any other region.
3. **Broad**: At least **35 %** of the base‑pair positions in the gene body must have a signal above the regional mean coverage. This excludes genes with a single sharp internal peak that inflates the average.

Genes shorter than 5000 bp are excluded; they are too small to reliably differentiate Broad from TSS.

## Workflow: Two‑stage strategy (following Young et al. 2011)

The script follows the original design: **first identify marked genes, then classify them**.

### Stage 1 – Identify “marked” genes (union of two methods)

*   **Peak‑based (MACS / SICER)**: uses the SICER2 island file (bed format) to find genes that **overlap a called peak**.
*   **Poisson exact test (edgeR)**: counts reads over each gene body and tests whether the ChIP signal is significantly higher than the Input. This captures **broad enrichment** that may be missed by peak callers.

The final set of genes to be classified is the **union** of these two methods.

### Stage 2 – Classify coverage patterns (BAM only)

Each selected gene is then passed to `classify_gene()`, which inspects the **raw read coverage profile from the BAM file** to decide among Promoter, TSS, or Broad.

> **Key point**: The classification step **does not use the bed file**. It relies exclusively on the coverage distribution obtained from the BAM file.

## Directory contents

| File | Description |
|------|-------------|
| `functions.R` | Complete library of helper functions for coverage extraction and visualisation |
| `ClassifyGenes.R` | Main classification script: reads BAM and bed files, calls `classify_gene`, outputs results |
| `classify_gene_core.R` | Stand‑alone core classification function (same as Thesis Appendix C) |
| `dispersion_test.R` | Test the impact of dispersion parameter (0.05, 0.1, 0.2) on classification results |
| `create_gene_bed.py` | Extract gene‑level coordinates from GTF annotation and save as BED format |
| `create_transcript_bed.py` | Extract transcript‑level coordinates from GTF annotation and save as BED format |
| `split_gene_class.sh` | Extract gene‑level BED coordinates for Broad/Promoter/TSS genes from classification results (used for motif analysis) |
| `split_gene_class_trans.sh` | Extract transcript‑level BED coordinates for Broad/Promoter/TSS genes from classification results (used for deepTools TSS/TES plots) |

### BED file generation workflow

Two types of BED files are generated from the classification results for different downstream analyses:

1. **Gene‑level BED** (for motif analysis):
   - `create_gene_bed.py` generates a gene‑level reference BED (`gencode.vM23.annotation.genes.bed`) from the GTF annotation.
   - `split_gene_class.sh` extracts Broad, Promoter, and TSS gene coordinates from this reference file.

2. **Transcript‑level BED** (for TSS/TES signal plots):
   - `create_transcript_bed.py` generates a transcript‑level reference BED (`geneWithDiff_trans.bed`) from the GTF annotation.
   - `split_gene_class_trans.sh` extracts Broad, Promoter, and TSS transcript coordinates from this reference file.

**Key distinction**:
- Motif analysis (`08_Motif_analysis/`) uses **gene‑level BED** because only gene IDs are needed to match peaks in the `.annotate` file.
- TSS/TES visualisation uses **transcript‑level BED** because precise transcription start/end coordinates are required for accurate signal profiling.
## Usage

### Quick start (core function only)

```r
library(GenomicRanges)
library(caTools)
source("functions.R")            # provides getRegion() etc.
source("classify_gene_core.R")   # loads classify_gene()

result <- classify_gene("ENSMUSG00000051951", annot, cover)
print(result$classification)     # "Promoter", "TSS", "Broad", or "Unclassified"
```

### Full pipeline

1.  Edit the **USER PARAMETERS** block in `ClassifyGenes.R` to point to your files.
2.  Run:

```bash
Rscript ClassifyGenes.R
```

The script will:

*   Load annotation from a TxDb object
*   Read ChIP and Input BAMs, remove PCR duplicates
*   Read the SICER2 island file (bed) and run Poisson exact test (edgeR)
*   Classify every gene that passes either the peak‑overlap filter or the Poisson test
*   Output `Gene_classification.txt`

### Input files

| File | Format | Parameter | Role |
|------|--------|-----------|------|
| ChIP BAM | BAM | `mapped_file_names[1]` | **Core** – provides coverage for classification |
| Input BAM | BAM | `mapped_file_names[2]` | Background for edgeR exactTest (optional) |
| Peak file | bed | `experiment_MACS_file` | **Auxiliary** – marks genes overlapping a called peak |

*   **On the peak file**: The parameter name `experiment_MACS_file` accepts **any bed‑formatted peak file** (SICER2 islands, MACS2 broadPeaks, etc.). In this study, SICER2 islands were converted to bed format.
*   **When Input is unavailable**: The Poisson test is skipped. Genes are selected solely by the bed file. The classification step remains unchanged because it only needs the ChIP BAM.

## Output

`Gene_classification.txt` is a tab‑separated table with columns:

| Column | Description |
|--------|-------------|
| `gene_id` | Ensembl gene ID |
| `Poisson` | `TRUE` if gene is significant by edgeR exactTest |
| `MACS` | `TRUE` if gene overlaps a called peak (bed input) |
| `Classification` | `Promoter`, `TSS`, `Broad`, `Unclassified`, or `NA` |

Intermediate R objects (`*_reads.RData`, `annotation.RData`, etc.) are also saved.

## Positive control validation

The pipeline was validated with public H3K27me3 data from mouse ES and G1ME cells (GEO: GSE27970). It reproduced the three distribution patterns reported by Young et al. (2011). See `06_Positive_control/`.

## Key parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `peak_window` | 200 bp | Smoothing window |
| `min_gene_length` | 5000 bp | Genes shorter than this are excluded |
| `Promoter_region` | c(-3000, -100) | Relative to TSS |
| `TSS_region` | c(-99, 1000) | Relative to TSS |
| `Broad_region` | c(1001, NA) | NA = end of gene |
| `Promoter_ratio` / `TSS_ratio` | 1.25 | Peak‑height dominance threshold |
| `Broad_ratio` | 0.35 | Fraction of gene body above mean coverage |
| `dispersion` | 0.05 | edgeR dispersion (results robust 0.05–0.2) |

## Software requirements

*   R ≥ 4.0.2
*   Packages: `GenomicRanges`, `GenomicAlignments`, `GenomicFeatures`, `Rsamtools`, `edgeR`, `caTools`
*   A TxDb object (e.g., `TxDb.Mmusculus.UCSC.mm10.ensGene`) or custom SQLite database from GTF

```r
install.packages("caTools")
BiocManager::install(c("GenomicRanges","GenomicAlignments","GenomicFeatures",
                       "Rsamtools","edgeR"))
```

## Key findings

*   **Broad genes**: most abundant, low expression, species‑specific functions
*   **Promoter/TSS genes**: narrow peaks, some remain highly expressed
*   **Cross‑species**: all three patterns detected in mouse, human, rat, bovine, porcine
*   **Dispersion robustness**: classification insensitive to dispersion choice

## Reference

Young MD, Willson TA, Wakefield MJ, et al. ChIP‑seq analysis reveals distinct H3K27me3 profiles that correlate with transcriptional activity. *Nucleic Acids Research*, 2011, 39(17): 7415–7427.
