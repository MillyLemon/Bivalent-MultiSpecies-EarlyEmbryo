#!/bin/bash
# Annotate bivalent, K4-only, and K27-only peaks using HOMER annotatePeaks.pl
# Usage: bash annotate_all_peaks.sh <genome> <stage1> <stage2> ...

GENOME=${1:-mm10}
shift  # 移除第一个参数，剩余的都是阶段名
STAGES="$@"

for stage in $STAGES; do
    echo "=== Processing ${stage} ==="

    # 1. Bivalent peaks
    biv_file="${stage}_given_${stage}_H3K27me3_peaks.broadPeak_${stage}_H3K4me3_peaks.broadPeak"
    if [[ -f "$biv_file" ]]; then
        annotatePeaks.pl "$biv_file" "$GENOME" -annStats "${stage}_biv.annStats" > "${stage}_biv.annotate"
        echo "  Bivalent: ${stage}_biv.annotate"
    else
        echo "  Bivalent: file not found, skipping"
    fi

    # 2. K27-only peaks
    k27_file="${stage}_given_${stage}_H3K27me3_peaks.broadPeak"
    if [[ -f "$k27_file" ]]; then
        annotatePeaks.pl "$k27_file" "$GENOME" > "${stage}_K27only.annotate"
        echo "  K27-only: ${stage}_K27only.annotate"
    else
        echo "  K27-only: file not found, skipping"
    fi

    # 3. K4-only peaks
    k4_file="${stage}_given_${stage}_H3K4me3_peaks.broadPeak"
    if [[ -f "$k4_file" ]]; then
        annotatePeaks.pl "$k4_file" "$GENOME" > "${stage}_K4only.annotate"
        echo "  K4-only: ${stage}_K4only.annotate"
    else
        echo "  K4-only: file not found, skipping"
    fi

    echo ""
done

echo "All done!"