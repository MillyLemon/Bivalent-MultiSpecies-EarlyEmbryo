#!/bin/bash
# Fix Bos taurus chromosome prefix and rebuild HOMER genome database
# Problem: SICER2 output used "chr1" format, but HOMER bosTau9 used "1" format.
# Solution: Add "chr" prefix to FASTA and GTF, then rebuild HOMER database.
# Usage: bash bosTau9_fix_chr_prefix.sh

GENOME_FA="/path/to/Bos_taurus.primary_assembly.fa"
GTF_FILE="/path/to/Bos_taurus.ARS-UCD1.2.104.gtf"

echo "=== Adding chr prefix to FASTA ==="
sed '/^>/s/^>/>chr/' ${GENOME_FA} | sed '/^>/s/ .*//' > Bos_taurus.primary_assembly.chr.fa

echo "=== Adding chr prefix to GTF ==="
sed '1,5s/^/#/' ${GTF_FILE} | sed '6,$s/^/chr/' > Bos_taurus.ARS-UCD1.2.104.chr.gtf

echo "=== Rebuilding HOMER bosTau9 database ==="
loadGenome.pl -name bosTau9 \
    -org cow \
    -fasta Bos_taurus.primary_assembly.chr.fa \
    -gtf Bos_taurus.ARS-UCD1.2.104.chr.gtf \
    -force

echo "=== Done ==="