# 02_RNA_seq_pipeline
RNA‑seq quantification pipeline using **Snakemake** + **Salmon** for processing merged single‑cell RNA‑seq data from early embryos across five mammalian species.

## Overview

This pipeline was used to quantify gene expression levels from scRNA‑seq data for all five species in the thesis. The core workflow is:

```
Raw FASTQ (single cells) → merge by stage → fastp → Salmon quant → gene‑level TPM
```
## Directory contents

| File | Description |
|------|-------------|
| `Snakefile` | Main Snakemake workflow: fastp trimming + Salmon quantification |
| `config.json` | Configuration file (Salmon index, gene mapping file, library type) |
| `run_local.sh` | Local execution script with automatic CPU detection |
| `merge_fastq.sh` | Generic script for merging fastq files from multiple single cells |
| `mergefq.sh` | Example merging commands for mouse dataset (GSE70605) |
| `../Files/mouse.tran2gene.txt` | Example gene mapping file (two‑column, transcript‑to‑gene) |

## Software requirements

| Software | Version | Purpose |
|----------|---------|---------|
| Snakemake | ≥5.0 | Workflow management and parallel processing |
| fastp | v0.23.2 | Read quality control and adapter trimming |
| Salmon | v1.8.0 | Transcript quantification (alignment‑free) |

## Quick start

### 1. Prepare input fastq files

Since each developmental stage contains multiple single cells, reads from all cells at the same stage are **merged** to obtain a representative expression profile.

**Using the generic merge script** (recommended for new datasets):
```bash
# Usage: bash merge_fastq.sh <output_stage_name> <comma‑separated_R1_list> <comma‑separated_R2_list>
bash merge_fastq.sh WT_2-cell "SRR1_R1.fq.gz,SRR2_R1.fq.gz" "SRR1_R2.fq.gz,SRR2_R2.fq.gz"
```

**Using the pre‑written example** for the mouse dataset (GSE70605):
```bash
bash mergefq.sh
```

After merging, place the resulting `*_1.fastq.gz` and `*_2.fastq.gz` files into a `fastq/` subdirectory. The Snakemake pipeline will automatically scan this directory for all paired‑end files.

### 2. Build Salmon index

Before quantification, build a **decoy‑aware index** for each species to improve accuracy:

```bash
# 1. Download transcriptome and genome FASTA files
wget <transcriptome_fasta_url>
wget <genome_fasta_url>
#wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/gencode.vM23.transcripts.fa.gz
#wget //ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/GRCm38.primary_assembly.genome.fa.gz

# 2. Create decoys file
grep "^>" <(gunzip -c genome.fa.gz) | cut -d " " -f 1 > decoys.txt
sed -i.bak -e 's/>//g' decoys.txt

# 3. Concatenate
cat transcripts.fa.gz genome.fa.gz > gentrome.fa.gz

# 4. Build index
salmon index -t gentrome.fa.gz -d decoys.txt -p 12 -i salmon_index --gencode
```

### 3. Prepare gene mapping file

The gene mapping file (e.g. `mouse.tran2gene.txt`) is a two‑column, tab‑separated file mapping transcript IDs to gene IDs:

```
ENSMUST00000193812.1	ENSMUSG00000102693.1
ENSMUST00000082908.1	ENSMUSG00000064842.1
```

This file can be generated from a GTF annotation using `transcript2gene.pl` (available in the `Files/` directory at the repository root):

```bash
perl ../Files/transcript2gene.pl annotation.gtf.gz | cut -f1,2 > mouse.tran2gene.txt
```
An example gene mapping file (`mouse.tran2gene.txt`) generated from the mouse GENCODE M23 annotation is provided in the `Files/` directory for reference.

### 4. Edit `config.json`

```json
{
    "idx": "/path/to/salmon_index",
    "gm": "mouse.tran2gene.txt",
    "libtype": "A"
}
```

- `idx`: Absolute path to the decoy‑aware Salmon index
- `gm`: Gene mapping file (two columns, tab‑separated, e.g. transcript_id → gene_id)
- `libtype`: `"A"` to let Salmon automatically infer library type

### 5. Run the pipeline

```bash
bash run_local.sh
```

The script will automatically create output folders, detect CPU cores (up to 16) and execute Snakemake with `--rerun-incomplete --keep-going`. All execution logs are saved in `logs/snakemake/`.

### 6. Output

For each sample, the following outputs are generated:

| Output | Path | Description |
|--------|------|-------------|
| Trimmed FASTQ | `fastp_fq/<sample>_r1_trimmed.fq.gz` | Quality‑controlled reads |
| QC Report | `fastp_report/<sample>_fastp.html` | fastp quality report |
| Transcript‑level | `salmon_output/<sample>/quant.sf` | Transcript abundance (TPM, NumReads) |
| Gene‑level | `salmon_output/<sample>/quant.genes.sf` | Gene abundance aggregated via `-g` |

Gene‑level TPM values are used for all downstream expression analyses.

## Gene‑level aggregation

Salmon's `-g` option aggregates transcript‑level estimates to the gene level:

- **Gene TPM** = Σ(transcript TPM)
- **Gene Length / EffectiveLength** = weighted average of transcript lengths, using transcript TPM as weights

This ensures that highly expressed transcripts contribute proportionally more to the gene‑level estimate.

## Reference genomes

| Species | Genome build | Transcriptome source |
|---------|-------------|---------------------|
| Mouse | mm10 (GRCm38) | GENCODE release M23 |
| Human | hg38 (GRCh38) | GENCODE release 43 |
| Rat | rn6 (Rnor_6.0) | Ensembl release 104 |
| Bovine | bosTau9 (ARS‑UCD1.2) | Ensembl release 104 |
| Porcine | susScr11 (Sscrofa11.1) | Ensembl release 104 |

## Notes

- The Snakemake pipeline automatically discovers samples by scanning `fastq/` for all `*_1.fastq.gz` files. It does **not** require a separate sample table.
- `libtype: "A"` tells Salmon to auto‑detect strandedness and read orientation.
- The pipeline uses `--validateMappings` and `--seqBias --gcBias` to improve quantification accuracy.
- All supplementary scripts (e.g., `transcript2gene.pl`) are located in the `Files/` directory at the repository root.