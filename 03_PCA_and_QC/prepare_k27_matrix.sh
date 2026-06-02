#!/bin/bash
# Preprocessing for H3K27me3: scoreisland → sorted BED → union BED → coverage matrix
# Usage: bash prepare_k27_matrix.sh
# Run this script in a directory containing all H3K27me3 scoreisland files

echo "=== Step 1: Converting scoreisland to sorted BED ==="
for file in *_f2q30_pmd-W200-G600.scoreisland; do
    base_name="${file%_f2q30_pmd-W200-G600.scoreisland}"
    output_name="${base_name}_sorted.bed"
    awk '{print $1"\t"$2"\t"$3}' "$file" | sort -k1,1 -k2,2n > "$output_name"
    echo "Processed: $file -> $output_name"
done

echo ""
echo "=== Step 2: Merging all K27 sorted BED files into union set ==="
cat *H3K27me3*sorted.bed | sort -k1,1 -k2,2n | bedtools merge -i - > k27.merged_peak_union.bed
echo "Union peaks: $(wc -l < k27.merged_peak_union.bed) regions"

echo ""
echo "=== Step 3: Building coverage matrix ==="
echo "peak_id" > header.txt
awk '{print $1 "_" $2 "_" $3}' k27.merged_peak_union.bed > peaks.txt
cat header.txt peaks.txt > k27.coverage_matrix.txt

for file in *H3K27me3*sorted.bed; do
    sample=$(basename "$file" .sorted.bed)
    echo "Processing $sample"

    bedtools coverage -a k27.merged_peak_union.bed -b "$file" -counts | \
        awk -v sample="$sample" 'BEGIN{print sample} {print $4}' > temp_col.txt

    paste k27.coverage_matrix.txt temp_col.txt > temp_matrix.txt
    mv temp_matrix.txt k27.coverage_matrix.txt
    rm temp_col.txt
done

rm header.txt peaks.txt
echo ""
echo "=== Done! ==="
echo "Output: k27.coverage_matrix.txt"
echo "Columns: $(head -1 k27.coverage_matrix.txt | awk '{print NF}')"