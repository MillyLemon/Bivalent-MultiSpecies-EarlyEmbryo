# Bivalent domains and H3K27me3 in early embryonic development across multiple species

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Introduction

This repository contains the complete analysis pipeline and code for the master's thesis:

**"Distribution patterns of Bivalent domains and H3K27me3 during early embryonic development across multiple species"**

*Li Xuemei, Binzhou Medical University, 2026*

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
