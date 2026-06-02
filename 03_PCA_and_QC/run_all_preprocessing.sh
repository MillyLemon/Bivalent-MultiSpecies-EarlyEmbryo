#!/bin/bash
# Master script: run K4 and K27 preprocessing sequentially
# Usage: bash run_all_preprocessing.sh
# Prerequisites:
#   - SICER2 scoreisland files must be present in subdirectories
#   - K4 and K27 files should be processed separately

echo "=========================================="
echo "PCA Preprocessing Pipeline"
echo "Started at: $(date)"
echo "=========================================="

# ---------- H3K4me3 ----------
echo ""
echo "### Processing H3K4me3 ###"
# If K4 files are in a separate directory, cd there first
# cd /path/to/k4_sicer2_output

if ls *H3K4me3*f2q30_pmd-W200-G600.scoreisland 1> /dev/null 2>&1; then
    bash prepare_k4_matrix.sh
    echo "K4 coverage matrix ready."
else
    echo "Warning: No H3K4me3 scoreisland files found in current directory."
    echo "Please run prepare_k4_matrix.sh in the directory containing K4 SICER2 output."
fi

# ---------- H3K27me3 ----------
echo ""
echo "### Processing H3K27me3 ###"
# If K27 files are in a separate directory, cd there first
# cd /path/to/k27_sicer2_output

if ls *H3K27me3*f2q30_pmd-W200-G600.scoreisland 1> /dev/null 2>&1; then
    bash prepare_k27_matrix.sh
    echo "K27 coverage matrix ready."
else
    echo "Warning: No H3K27me3 scoreisland files found in current directory."
    echo "Please run prepare_k27_matrix.sh in the directory containing K27 SICER2 output."
fi

echo ""
echo "=========================================="
echo "Preprocessing finished at: $(date)"
echo "Output files: k4.coverage_matrix.txt, k27.coverage_matrix.txt"
echo "These matrices are ready for PCA analysis (see pca_rpkm.R)."
echo "=========================================="