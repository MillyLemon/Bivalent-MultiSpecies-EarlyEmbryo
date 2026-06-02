#!/bin/bash
# Count unique genes for each peak category from filtered annotate files
# Usage: bash stat_peak.gene.sh <output_prefix>
# Example: bash stat_peak.gene.sh mouse1

PREFIX=${1:-output}

echo -e "sample\tGeneCounts\tGroup" > ${PREFIX}.gene.stat

# K27-only
for f in *K27only*genic.annotate; do
    [ -f "$f" ] || continue
    s=$(echo "$f" | sed 's/_K27only.*//')
    c=$(awk -F'\t' 'NR>1 && $16!="."{print $16}' "$f" | sort -u | wc -l)
    echo -e "${s}_K27Me3\t$c\tk27_gene"
done >> ${PREFIX}.gene.stat

# K4-only
for f in *K4only*genic.annotate; do
    [ -f "$f" ] || continue
    s=$(echo "$f" | sed 's/_K4only.*//')
    c=$(awk -F'\t' 'NR>1 && $16!="."{print $16}' "$f" | sort -u | wc -l)
    echo -e "${s}_K4Me3\t$c\tk4_gene"
done >> ${PREFIX}.gene.stat

# Bivalent
for f in *given*genic.annotate; do
    [ -f "$f" ] || continue
    s=$(echo "$f" | cut -d'_' -f1)
    c=$(awk -F'\t' 'NR>1 && $16!="."{print $16}' "$f" | sort -u | wc -l)
    echo -e "${s}_biv\t$c\tbiv_gene"
done >> ${PREFIX}.gene.stat

echo "Done! Output: ${PREFIX}.gene.stat"
cat ${PREFIX}.gene.stat