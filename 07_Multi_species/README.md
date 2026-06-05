# 07_Multi_species

Species‑specific configurations and adaptations required to extend the analysis pipeline from mouse to human, rat, bovine, and porcine.

## Overview

While the core analysis logic is species‑independent, several tools require explicit species configuration:

- **SICER2**: Uses an internal `GenomeData.py` file to look up chromosome names and lengths. Only a limited set of genomes is pre‑configured.
- **HOMER**: Genome databases must be installed separately for each species via `configureHomer.pl` or `loadGenome.pl`.
- **Salmon**: Requires decoy‑aware transcriptome indices built from species‑specific reference genomes and annotations.

## Directory contents

| File | Description |
|------|-------------|
| `add_genome_sicer2.py` | Chromosome lists and lengths added to SICER2 `GenomeData.py` for bosTau9 and rn6 |
| `install_homer_genomes.sh` | Commands to install HOMER genome databases for all five species |
| `bosTau9_fix_chr_prefix.sh` | Fix Bos taurus chromosome naming inconsistency (add `chr` prefix) and rebuild HOMER database |

## Species reference genomes

| Species | Genome build | Assembly | Annotation source |
|---------|-------------|----------|-------------------|
| Mouse | mm10 | GRCm38 | GENCODE release M23 |
| Human | hg38 | GRCh38 | GENCODE release 43 |
| Rat | rn6 | Rnor_6.0 | Ensembl release 104 |
| Bovine | bosTau9 | ARS‑UCD1.2 | Ensembl release 104 |
| Porcine | susScr11 | Sscrofa11.1 | Ensembl release 104 |

## Key adaptations by species

### Rat (rn6)
- SICER2: Added chromosome list and lengths to `GenomeData.py`
- HOMER: Standard installation via `configureHomer.pl -install rn6`
- Salmon: Decoy‑aware index built from Ensembl release 104

### Bovine (bosTau9)
- SICER2: Added chromosome list and lengths to `GenomeData.py`. **Note**: The original `GenomeData.py` only contained bosTau8 (UMD3.1.1). All analysis was unified to bosTau9 (ARS‑UCD1.2).
- HOMER: Custom installation via `loadGenome.pl`. Required chromosome prefix correction (`1` → `chr1`) because SICER2 output used UCSC‑style chromosome names while the Ensembl FASTA/GTF used numeric chromosome names.
- Salmon: Decoy‑aware index built from Ensembl release 104

### Porcine (susScr11)
- SICER2: Standard genome already supported
- HOMER: Standard installation via `configureHomer.pl -install susScr11`
- Salmon: Decoy‑aware index built from Ensembl release 104

### Human (hg38)
- SICER2: Standard genome already supported
- HOMER: Standard installation via `configureHomer.pl -install hg38`
- Salmon: Decoy‑aware index built from GENCODE release 43

## Lessons learned

1. **Always verify genome version consistency** across all tools (aligner, peak caller, annotation). The bovine dataset initially had a mismatch between SICER2 (bosTau8) and other tools (bosTau9).
2. **HOMER custom genomes**: The `loadGenome.pl` approach works for species not available via `configureHomer.pl`, but requires careful handling of chromosome naming conventions.
3. **Chromosome naming**: UCSC‑style (`chr1`) vs Ensembl‑style (`1`) is a common source of silent failures. Always check the first few lines of BAM, bed, and annotation files to confirm consistency.

## Software requirements

- SICER2 v1.0.3 (with modified `GenomeData.py`)
- HOMER v4.11
- Salmon v1.8.0
- HISAT2 v2.2.1