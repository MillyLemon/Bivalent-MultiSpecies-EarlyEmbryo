# Bivalent domains and H3K27me3 in early embryonic development across multiple species

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Introduction

This repository contains the complete analysis pipeline and code for the master's thesis:

**"Distribution patterns of Bivalent domains and H3K27me3 during early embryonic development across multiple species"**

*Li Xuemei, Shandong Medical and Pharmaceutical University, 2026*

The project systematically investigates whether bivalent domains (H3K4me3/H3K27me3 co‑marked regions) are universally present across mammalian species during early embryogenesis, how their modification patterns change dynamically, and whether H3K27me3 exhibits non‑canonical narrow distribution patterns *in vivo*. The analysis integrates publicly available ChIP‑seq and RNA‑seq data from five mammalian species: mouse, human, rat, bovine, and porcine.

The code of main figures in this paper can be found in the `Figures` directory, and some necessary annotation files can be found in the `Files` directory.

## Raw Public Sequencing Data

All raw data (fastq files) were downloaded from the accession numbers listed below.

### Mouse

1. Liu X, Wang C, Liu W, Li J et al. ***Distinct features of H3K4me3 and H3K27me3 chromatin domains in pre‑implantation embryos***. Nature 2016 Sep.
   - ULI‑NChIP–seq: [GSE73952](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE73952), SRP064741 (67 samples; we selected 36 samples containing H3K4me3 and H3K27me3)
   - RNA‑seq: [GSE70605](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70605), SRP060478 (49 WT samples)

**Table 1 – Samples selected in mouse1**

| Type | 2‑cell | 4‑cell | 8‑cell | morula | ICM | TE |
|:----:|:------:|:------:|:------:|:------:|:---:|:--:|
| H3K4me3 | 3 | 3 | 2 | 3 | 2 | 2 |
| H3K27me3 | 3 | 3 | 3 | 2 | 2 | 2 |
| Input | 1 | 1 | 1 | 1 | 1 | 1 |
| RNA‑seq | 5 | 6 | 6 | 6 | 5 | 6 |

2. Zhang B, Zheng H, Huang B, Li W et al. ***Allelic reprogramming of the histone modification H3K4me3 in early mammalian development***. Nature 2016.
   - STAR ChIP‑seq: [GSE71434](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE71434), SRP062106 (H3K4me3)
   - RNA‑seq: [GSE71434](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE71434)

3. Zheng H, Huang B, Zhang B, Xiang Y et al. ***Resetting Epigenetic Memory by Reprogramming of Histone Modifications in Mammals***. Mol Cell 2016.
   - STAR ChIP‑seq: [GSE76687](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76687), SRP068549 (H3K27me3)
   - RNA‑seq: [GSE71434](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE71434)

**Table 2 – Samples selected in mouse2**

| Type | zygote | 2‑cell early | 2‑cell late | 8‑cell | ICM |
|:----:|:------:|:------------:|:-----------:|:------:|:---:|
| H3K4me3 | 2 | 2 | 2 | 2 | 2 |
| H3K27me3 | 2 | 1 | 2 | 1 | 2 |
| RNA‑seq | 2 | 2 | 2 | 2 | 2 |

### Human

1. Xia W, Xu J, Yu G et al. ***Resetting histone modifications during human parental‑to‑zygotic transition***. Science 2019.
   - CUT&RUN: [GSE124718](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE124718), SRP176396

2. Wu J, Xu J, Liu B et al. ***Chromatin analysis in human early development reveals epigenetic transition during ZGA***. Nature 2018.
   - RNA‑seq: [GSE101571](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE101571), SRP112718

As paper 1 and paper 3 are from the same laboratory, we selected these two papers for the final analysis.

**Table 3 – Samples selected in human**

| Type | GV | 4‑cell 3PN | 8‑cell 2PN | 8‑cell 3PN | ICM 2PN |
|:----:|:--:|:---------:|:----------:|:----------:|:-------:|
| H3K4me3 | 2 | 1 | 2 | 2 | 2 |
| H3K27me3 | 2 | 1 | 1 | 2 | 2 |
| RNA‑seq | 2 | 2 | 2 | 2 | 2 |

### Rat, Bovine, Porcine

Lu X, Zhang Y, Wang L et al. ***Evolutionary epigenomic analyses in mammalian early embryos reveal species‑specific innovations and conserved principles of imprinting***. Sci Adv 2021.
[GSE163620](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163620), SRP298725

**Table 4 – Samples selected in rat**

| Type | FGO | MII | 1‑cell | 2‑cell | 4‑cell | 8‑cell | blastula |
|:----:|:---:|:---:|:------:|:------:|:------:|:------:|:--------:|
| H3K4me3 | 3 | 4 | 2 | 2 | 2 | 2 | 2 |
| H3K27me3 | 3 | 3 | 2 | 2 | 2 | 2 | 2 |
| Smart‑seq | 4 | 4 | 4 | 2 | 2 | 2 | 2 |

**Table 5 – Samples selected in bovine**

| Type | FGO | MII | 4‑cell | 8‑cell | 16‑cell | blastula |
|:----:|:---:|:---:|:------:|:------:|:-------:|:--------:|
| H3K4me3 | 2 | 2 | 2 | 2 | 2 | 2 |
| H3K27me3 | 2 | 2 | 2 | 2 | 2 | 2 |
| Smart‑seq | 2 | 2 | 2 | 2 | 2 | 3 |

**Table 6 – Samples selected in porcine**

| Type | FGO | MII | 2‑cell | 4‑cell | 8‑cell | morula | blastula |
|:----:|:---:|:---:|:------:|:------:|:------:|:------:|:--------:|
| H3K4me3 | 2 | 2 | 2 | 1 | 1 | 2 | 2 |
| H3K27me3 | 2 | 2 | 2 | 2 | 2 | 2 | 2 |
| Smart‑seq | 2 | 2 | 2 | 2 | 2 | 4 | 3 |

## Requirements

The software packages needed for the analysis are:
```
R v4.2.1
perl v5.16.3
python3
fastp v0.23.2
hisat2 v2.2.1
samtools v1.9
picard v2.25.5
macs2 v2.2.7
sicer2 v1.0.3
bedtools v2.30.0
homer v4.11
salmon v1.8.0
deepTools v3.5.5
```

All of them can be found in the `Packages` directory in this repository.

## 1. ChIP‑seq pipeline (using 2‑cell H3K27me3 rep1 as an example)

### 1.1 Data quality control
```bash
fastp -i SRR3208756_1.fastq.gz \
      -I SRR3208756_2.fastq.gz \
      -o 2_cell_H3K27me3_rep1_1_trim.fq.gz \
      -O 2_cell_H3K27me3_rep1_2_trim.fq.gz \
      -l 25 --detect_adapter_for_pe \
      -h 2_cell_H3K27me3_rep1_fastp.html \
      -j 2_cell_H3K27me3_rep1_fastp.json
```

### 1.2 Genome mapping
```bash
hisat2 -X 500 -p 4 \
       -x mm10_index \
       --no-temp-splicesite \
       --no-spliced-alignment \
       -3 1 --summary-file 2_cell_H3K27me3_rep1_aln_sum.txt \
       -1 2_cell_H3K27me3_rep1_1_trim.fq.gz \
       -2 2_cell_H3K27me3_rep1_2_trim.fq.gz | \
       samtools view -ShuF 4 -f 2 -q 30 - | \
       samtools sort -o 2_cell_H3K27me3_rep1_f2q30.sorted.bam
```

### 1.3 Blacklist removal
Blacklist: The ENCODE blacklist regions to exclude from analysis. See this publication and this repository for details.
```bash
intersectBed -a 2_cell_H3K27me3_rep1_f2q30.sorted.bam \
             -b mm10.blacklist.bed -v \
             > 2_cell_H3K27me3_rep1_final.bam
```

### 1.4 Duplication removal
```bash
java -Xmx12g \
     -jar picard.jar MarkDuplicates \
     INPUT=2_cell_H3K27me3_rep1_final.bam \
     OUTPUT=2_cell_H3K27me3_rep1_f2q30_pmd.bam \
     REMOVE_DUPLICATES=true ASSUME_SORTED=true \
     METRICS_FILE=2_cell_H3K27me3_rep1_f2q30_pmd.out \
     2>2_cell_H3K27me3_rep1.log
```

### 1.5 Peak calling
We compared MACS2 and SICER2 for peak calling and ultimately selected SICER2 for its higher sensitivity and stability in detecting broad domains such as H3K27me3.

Option A – SICER2 (used in final analysis):
```bash
sicer -t 2_cell_H3K27me3_rep1_f2q30_pmd.bam \
      -c 2_cell_input_f2q30_pmd.bam \
      -s mm10
```

Option B – MACS2 (used for comparison):
```bash
# PE
macs2 callpeak -t 2_cell_H3K27me3_rep1_f2q30_pmd.bam \
               -c 2_cell_input_f2q30_pmd.bam \
               -f BAMPE -g mm --keep-dup all \
               --broad --broad-cutoff 0.01 \
               --nolambda -n 2_cell_H3K27me3_rep1 \
               --outdir macs2_pk -B --SPMR --shift 0
```

### 1.6 Prepare bigwig files for track visualization
Bigwig files were generated using deepTools with CPM normalisation for visualisation:
```bash
bamCoverage -b 2_cell_H3K27me3_rep1_f2q30_pmd.bam \
            -o 2_cell_H3K27me3_rep1_uniform_25bp_CPM.bw \
            --binSize 25 \
            --normalizeUsing CPM \
            --smoothLength 75 \
            --extendReads 200 \
            --ignoreForNormalization chrX chrY chrM \
            --numberOfProcessors 4
```
### 1.7 PCA quality control and sample filtering
After peak calling, we used RPKM‑based PCA to assess replicate reproducibility. The PCA matrix was constructed by merging all broadPeak files into a union set, then calculating RPKM for each sample at each peak region.

Region RPKM was calculated as:

$$
\text{Region RPKM} = \frac{\text{ReadsWithinRegion}}{\frac{\text{RegionLength}}{1000} \times \text{SampleReadsDepth}}
$$
where
$$
\text{SampleReadsDepth} = \frac{\text{SampleReads}}{1,000,000}
$$



