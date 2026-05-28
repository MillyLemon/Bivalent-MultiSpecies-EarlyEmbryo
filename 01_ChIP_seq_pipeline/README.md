# 01_ChIP_seq_pipeline

Standard ChIP‑seq data processing pipeline from raw FASTQ files to cleaned BAM files ready for peak calling.

## Pipeline overview
FASTQ → fastp → hisat2 → intersectBed → picard → BAM

## Scripts

| Script | Description |
|--------|-------------|
| `chip_seq_pipeline_PE.sh` | Complete pipeline for paired‑end reads |
| `chip_seq_pipeline_SE.sh` | Complete pipeline for single‑end reads |
| `merge_replicates.sh` | Merge biological replicate BAM files before peak calling |

## Software requirements

| Software | Version | Purpose |
|----------|---------|---------|
| fastp | v0.23.2 | Read quality control and adapter trimming |
| HISAT2 | v2.2.1 | Genome alignment |
| SAMtools | v1.9 | BAM file sorting and filtering |
| Picard | v2.25.5 | PCR duplicate removal |
| BEDTools | v2.30.0 | Blacklist region removal |

## Usage

### Paired‑end sample

```bash
bash chip_seq_pipeline_PE.sh \
     sample_R1.fastq.gz \
     sample_R2.fastq.gz \
     sample_name \
     /path/to/hisat2_index/genome \
     /path/to/mm10.blacklist.bed
```
### Single‑end sample
```bash
bash chip_seq_pipeline_SE.sh \
     sample.fastq.gz \
     sample_name \
     /path/to/hisat2_index/genome \
     /path/to/mm10.blacklist.bed
```
### Merging replicates
After processing individual replicates, merge them for each stage:
```bash
bash merge_replicates.sh 2_cell_H3K27me3 \
     2_cell_H3K27me3_rep1.pmd.bam \
     2_cell_H3K27me3_rep2.pmd.bam \
     2_cell_H3K27me3_rep3.pmd.bam
```
### Peak calling
Peak calling was performed on merged BAM files using SICER2 (v1.0.3):
```bash
sicer -t sample.merged.bam -c input.pmd.bam -s mm10
```
SICER2 was selected over MACS2 after systematic comparison because of its higher sensitivity and better replicate consistency for broad histone marks such as H3K27me3. See `Figures/` for the comparison results.

### key parameters
```text
Step     Parameter               Value  Rationale
fastp    -l                      25     Minimum read length after trimming
hisat2   -X                      500    Maximum fragment length for paired‑end
hisat2   --no-spliced-alignment  on     ChIP‑seq uses genomic DNA, no splicing
samtools -q                      30     Mapping quality filter
samtools -f 2                    on     Keep only properly paired reads (PE only)
picard   REMOVE_DUPLICATES       true   Remove PCR duplicates
```
### Output directory structure
```
ChIP_seq/
├── fastp_fq/           # Trimmed FASTQ files
├── fastp_report/       # Quality control reports (HTML/JSON)
├── hisat2_log/         # Alignment summary files
├── hisat2_bam/         # Sorted BAM files (pre‑filter)
├── picard_bam/         # Final deduplicated BAM files
└── picard_log/         # Picard metrics files
```
### Notes
Blacklist files for mouse (mm10) and human (hg38) can be downloaded from the [Boyle Lab repository](https://github.com/Boyle-Lab/Blacklist). Other species were analysed without blacklist filtering as ENCODE blacklists are not available.

All output BAM files in `picard_bam/` are ready for downstream peak calling (SICER2) and signal visualisation (deepTools).

For samples without input controls, SICER2 can be run without the `-c` option with a manual E‑value threshold.


