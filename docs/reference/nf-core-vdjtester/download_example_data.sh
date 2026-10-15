#!/bin/bash

# Script to download example data for nf-core-vdjtester tutorial
# Based on: https://www.10xgenomics.com/support/software/cell-ranger/8.0/tutorials/cr-tutorial-vdj

set -e  # Exit on error

echo "=== nf-core-vdjtester Example Data Download Script ==="
echo ""

WORKING_DIR="$(pwd)"
PROJECT_NAME="dataset-vdj-practice"

# Create project directory structure
echo "1. Creating project directory..."
mkdir -p "${WORKING_DIR}/${PROJECT_NAME}"
cd "${WORKING_DIR}/${PROJECT_NAME}"

# Download FASTQ files
echo "2. Downloading input FASTQ files..."
curl --fail -L \
    https://cf.10xgenomics.com/samples/cell-vdj/6.0.0/sc5p_v2_hs_B_1k_multi_5gex_b_Multiplex/sc5p_v2_hs_B_1k_multi_5gex_b_Multiplex_fastqs.tar \
    -o sc5p_v2_hs_B_1k_multi_5gex_b_Multiplex_fastqs.tar

# Decompress FASTQs
echo "3. Extracting FASTQ tarball..."
tar -xf sc5p_v2_hs_B_1k_multi_5gex_b_Multiplex_fastqs.tar

# Clean up tarball after extraction
rm sc5p_v2_hs_B_1k_multi_5gex_b_Multiplex_fastqs.tar

echo "   ✓ FASTQs extracted:"
ls -la sc5p_v2_hs_B_1k_multi_5gex_b_fastqs/

# Change back to working directory
cd "${WORKING_DIR}"

# Download V(D)J reference
echo ""
echo "4. Downloading V(D)J reference..."
curl --fail -L \
    https://cf.10xgenomics.com/supp/cell-vdj/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0.tar.gz \
    -o refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0.tar.gz

# Extract reference
echo "5. Extracting V(D)J reference..."
tar -xf refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0.tar.gz
rm refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0.tar.gz

echo ""
echo "=== Download Complete! ==="
echo ""
echo "Directory structure created:"
echo "${WORKING_DIR}/"
echo "├── ${PROJECT_NAME}/"
echo "│   └── sc5p_v2_hs_B_1k_multi_5gex_b_fastqs/"
echo "└── refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0/"
echo ""