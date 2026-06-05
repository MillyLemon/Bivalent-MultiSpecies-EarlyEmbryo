#!/bin/bash
# Generate TSS and TES profile plots for ES and G1ME validation data
# Usage: bash plot_validation_profiles.sh
# Prerequisites:
#   - BigWig files generated from processed BAM files (see 01_ChIP_seq_pipeline/)
#   - Gene classification results (broad_genes.trans.bed, promoter_genes.trans.bed, tss_genes.trans.bed)
#     placed in subdirectories named after each sample

BW_DIR="../01_ChIP_seq_pipeline/picard_bam"

for sample in GSM691591_ES.H3K27 GSM691594_G1ME.H3K27; do
    echo "=== Processing ${sample} ==="

    BW="${BW_DIR}/${sample}.uniform_25bp_CPM.bw"
    BROAD="${sample}/broad_genes.trans.bed"
    PROMOTER="${sample}/promoter_genes.trans.bed"
    TSS_BED="${sample}/tss_genes.trans.bed"

    # Check input files
    if [[ ! -f "${BW}" ]]; then echo "Missing: ${BW}"; continue; fi
    if [[ ! -f "${BROAD}" ]]; then echo "Missing: ${BROAD}"; continue; fi
    if [[ ! -f "${PROMOTER}" ]]; then echo "Missing: ${PROMOTER}"; continue; fi
    if [[ ! -f "${TSS_BED}" ]]; then echo "Missing: ${TSS_BED}"; continue; fi

    # TSS analysis
    computeMatrix reference-point \
        --referencePoint TSS \
        -S "${BW}" \
        -R "${BROAD}" "${PROMOTER}" "${TSS_BED}" \
        -b 30000 -a 30000 \
        --skipZeros \
        -o ${sample}/${sample}_matrix_trans_TSS.gz

    # TES analysis
    computeMatrix scale-regions \
        -S "${BW}" \
        -R "${BROAD}" "${PROMOTER}" "${TSS_BED}" \
        -b 30000 -a 30000 \
        --binSize 25 \
        --regionBodyLength 30000 \
        --skipZeros \
        -o ${sample}/${sample}_matrix_trans_TES.gz

    # TSS profile
    plotProfile -m ${sample}/${sample}_matrix_trans_TSS.gz \
        -out ${sample}/${sample}_trans_TSS_profile.pdf \
        --regionsLabel "Broad" "Promoter" "TSS" \
        --legendLocation best \
        --samplesLabel "" \
        --plotTitle "TSS Signal - ${sample}"

    # TES profile
    plotProfile -m ${sample}/${sample}_matrix_trans_TES.gz \
        -out ${sample}/${sample}_trans_TES_profile.pdf \
        --regionsLabel "Broad" "Promoter" "TSS" \
        --legendLocation best \
        --samplesLabel "" \
        --plotTitle "TES Signal - ${sample}"

    # TSS heatmap
    plotHeatmap -m ${sample}/${sample}_matrix_trans_TSS.gz \
        -out "${sample}/heatmap_${sample}.trans.TSS.pdf" \
        --colorMap RdBu \
        --whatToShow 'heatmap and colorbar' \
        --regionsLabel "Broad" "Promoter" "TSS" \
        --legendLocation none

    # TES heatmap
    plotHeatmap -m ${sample}/${sample}_matrix_trans_TES.gz \
        -out "${sample}/heatmap_${sample}.trans.TES.pdf" \
        --colorMap RdBu \
        --whatToShow 'heatmap and colorbar' \
        --regionsLabel "Broad" "Promoter" "TSS" \
        --legendLocation none

done

echo "=== All done ==="