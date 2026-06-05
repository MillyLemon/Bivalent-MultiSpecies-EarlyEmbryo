#!/bin/bash
# Install HOMER genome databases for all five species
# Usage: bash install_homer_genomes.sh

echo "=== Installing HOMER genome databases ==="

# Standard genomes (available via configureHomer.pl)
perl /path/to/homer/configureHomer.pl -install hg38
perl /path/to/homer/configureHomer.pl -install mm10
perl /path/to/homer/configureHomer.pl -install rn6
perl /path/to/homer/configureHomer.pl -install susScr11

# bosTau9 requires manual installation (see bosTau9_fix_chr_prefix.sh)
perl /path/to/homer/configureHomer.pl -install cow

echo "=== Done ==="
echo "Note: bosTau9 was installed manually via loadGenome.pl (see bosTau9_fix_chr_prefix.sh)"