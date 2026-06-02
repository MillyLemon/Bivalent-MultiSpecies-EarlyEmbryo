#!/bin/bash
# Filter intergenic peaks and generate annotation statistics for all peak types
# Usage: bash filter_all_intergenic.sh

echo -e "File\tTotal\tIntergenic\tPromoter\tExon\tIntron\t3UTR\t5UTR\tTTS\tNonCoding" > stats.all.txt

for file in *.annotate; do
    [ -f "$file" ] || continue
    echo "Processing: $file"

    awk 'NR>1 {
        split($8, a, "(")
        anno = a[1]
        gsub(/ /, "", anno)

        if (anno ~ /Intergenic/) inter++
        else if (anno ~ /promoter/) pro++
        else if (anno ~ /exon/) exon++
        else if (anno ~ /intron/) intron++
        else if (anno ~ /3.UTR/) utr3++
        else if (anno ~ /5.UTR/) utr5++
        else if (anno ~ /TTS/) tts++
        else if (anno ~ /non-coding/) nc++
        else other++
    }
    END {
        print FILENAME, NR-1, inter+0, pro+0, exon+0, intron+0, utr3+0, utr5+0, tts+0, nc+0
    }' "$file" >> stats.all.txt

    # Save genic-only version
    awk 'NR==1 || $8 !~ /Intergenic/' "$file" > "${file%.annotate}_genic.annotate"
done

echo "Done!"
column -t stats.all.txt