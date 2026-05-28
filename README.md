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

where `SampleReadsDepth = SampleReads / 1,000,000`.
```bash
# Count reads per sample
samtools view -c 2_cell_H3K27me3_rep1_f2q30_pmd.bam
```

Outlier samples identified by PCA (e.g., 4‑cell H3K4me3 rep2) were removed. Low‑abundance peaks with RPKM < 0.5 were filtered, and peaks present in at least two independent samples were retained.

PCA plotting code and demo input files can be found in the Figures directory and the Files directory.

### 1.8 Merge biological replicates and call peaks
After quality filtering, replicates were merged for each stage:
```bash
samtools merge -@ 16 2_cell_H3K27me3.merged.bam \
               2_cell_H3K27me3_rep1_f2q30_pmd.bam \
               2_cell_H3K27me3_rep3_f2q30_pmd.bam
```
SICER2 peak calling was performed on merged BAM files as described in 1.5.

### 1.9 Identify bivalent / K27‑only / K4‑only peaks (HOMER: mergePeaks)
```bash
mergePeaks -d given \
           -venn 2_cell_given.venn \
           -prefix 2_cell_given \
           2_cell_H3K27me3_peaks.broadPeak \
           2_cell_H3K4me3_peaks.broadPeak

# The overlapping peaks are bivalent domains.
# The remaining peaks after removing overlaps are K27‑only or K4‑only peaks.
```

### 1.10 Annotate peaks (HOMER: annotatePeaks.pl)
```bash
annotatePeaks.pl 2_cell_bivalent mm10 \
                 -annStats 2_cell_bivalent.annStats \
                 > 2_cell_bivalent.annotate

# Remove intergenic peaks to retain gene‑associated peaks
perl filter.peak.intergenic.pl 2_cell_bivalent.annotate \
     > 2_cell_bivalent.annotate.filter
```

## 2. RNA‑seq pipeline
### 2.1 Merge raw data
For single‑cell RNA‑seq data, reads from multiple cells at the same developmental stage were merged (using WT 2‑cell as an example):
```bash
cat SRR2089603_1.fastq.gz SRR2089604_1.fastq.gz SRR2089605_1.fastq.gz \
    SRR2089606_1.fastq.gz SRR2089607_1.fastq.gz > WT_2-cell_1.fastq.gz

cat SRR2089603_2.fastq.gz SRR2089604_2.fastq.gz SRR2089605_2.fastq.gz \
    SRR2089606_2.fastq.gz SRR2089607_2.fastq.gz > WT_2-cell_2.fastq.gz
```

### 2.2 Transcript quantification by Salmon
Salmon is a tool for rapid transcript quantification from RNA‑seq data. See the [Salmon documentation](https://salmon.readthedocs.io/en/latest/salmon.html) for details.

#### 2.2.1 Indexing
We used the decoy‑aware index approach following [this tutorial](https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/).
```bash
# Mouse (mm10)
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/gencode.vM23.transcripts.fa.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/GRCm38.primary_assembly.genome.fa.gz
grep "^>" <(gunzip -c GRCm38.primary_assembly.genome.fa.gz) | cut -d " " -f 1 > decoys.txt
sed -i.bak -e 's/>//g' decoys.txt
cat gencode.vM23.transcripts.fa.gz GRCm38.primary_assembly.genome.fa.gz > gentrome.fa.gz
salmon index -t gentrome.fa.gz -d decoys.txt -p 12 -i salmon_mm_index --gencode
```
Reference genomes and annotations for other species:

| Species | Genome source | Annotation source |
|---------|---------------|-------------------|
| Human (hg38) | GENCODE release 43 | gencode.v43.annotation.gtf.gz |
| Rat (rn6) | Ensembl release 104 | Rattus_norvegicus.Rnor_6.0.104.gtf.gz |
| Bovine (bosTau9) | Ensembl release 104 | Bos_taurus.ARS-UCD1.2.104.gtf.gz |
| Porcine (susScr11) | Ensembl release 104 | Sus_scrofa.Sscrofa11.1.104.gtf.gz |

#### 2.2.2 Quantification
```bash
salmon quant -i salmon_mm_index \
             -g gencode.vM23.annotation.gtf.gz \
             --gcBias -l A \
             -1 WT_2-cell_1.fastq.gz \
             -2 WT_2-cell_2.fastq.gz \
             -p 8 \
             -o ./quants/WT_2-cell
```
### 2.3 Convert transcript quantification to gene‑level
```r
library(GenomicFeatures)
library(tximport)

gtf_file <- 'gencode.vM23.annotation.gtf.gz'
txdb <- makeTxDbFromGFF(gtf_file)
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")

count_dir <- './quants'
name <- c('WT_2-cell','WT_4-cell','WT_8-cell','WT_Blastocyst',
          'WT_ICM','WT_MII','WT_Morula','WT_TE','WT_Zygote')
files <- file.path(count_dir, name, "quant.sf")
names(files) <- paste0(name)

txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene)
tpm <- txi.salmon$abundance
write.table(tpm, 'mouse.GSE73952.gene.tpm', sep='\t', quote=F)
```
### 2.4 Convert Ensembl ID to gene name
```bash
# Download GTF file from Ensembl
wget ftp://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/Mus_musculus.GRCm38.102.gtf.gz

# Extract transcript-to-gene mapping
perl transcript2gene.pl Mus_musculus.GRCm38.102.gtf.gz \
     > Mus_musculus.GRCm38.102.gene2trans.info

# Merge with TPM matrix
perl merge.geneinfo.pl Mus_musculus.GRCm38.102.gene2trans.info \
     mouse.GSE73952.gene.tpm > mouse.GSE73952.gene.symbol.tpm
```
The file `Mus_musculus.GRCm38.102.gene2trans.info` contains three tab‑separated columns:
`ENSMUST00000082908 ENSMUSG00000064842 Gm26206`

### 2.5 Merge ChIP‑seq peak annotation with gene expression
```bash
perl merge.chip.rna.pl 2_cell_H3K27me3.merged.Peak.annotate.filter \
     mouse.GSE73952.gene.symbol.tpm \
     > 2_cell_H3K27me3.merged.Peak.annotate.filter.exp
```
## 3. Gene classification of H3K27me3‑modified genes
This method is based on the classification framework described in Young et al. (2011) Nucleic Acids Research, which categorises H3K27me3‑modified genes into three enrichment patterns: Promoter, TSS, and Broad.

### 3.1 Classification logic
Each gene is divided into three regions:

- Promoter region: TSS upstream 3000 bp to TSS upstream 100 bp

- TSS region: TSS upstream 100 bp to TSS downstream 1000 bp

- Broad region: TSS downstream 1001 bp to gene termination site

A 200 bp sliding window is used to smooth the coverage signal within each region, and the maximum peak height is determined for each region. Each gene is assigned to the region with the highest mean coverage, subject to the following criteria:

1. Promoter: The peak in the promoter region must be at least 25% higher than the maximum peak in any other region.

2. TSS: The peak in the TSS region must be at least 25% higher than the maximum peak in any other region.

3. Broad: At least 35% of the base‑pair positions in the gene body must have a signal higher than the mean signal across the region (to exclude genes with a single sharp internal peak driving a high average).

