process DOWNLOAD_CELLRANGER_DATA {
    tag "CellRanger Data Download"
    storeDir("${params.outdir}/cellranger_data")
    publishDir("${params.outdir}/cellranger_downloads", mode: 'copy')

    output:
        path "*", emit: fastqs
        path "refdata-gex-GRCh38-*", emit: gex_reference
        path "refdata-cellranger-vdj-*", emit: vdj_reference
        path "cellranger_config.csv", emit: config_csv

    script:
    """
    set -eu

    # Create output directories
    mkdir -p fastq_samples
    mkdir -p references

    echo "Starting CellRanger example data download..."

    # =====================================================
    # STEP 1: Download Input FASTQ Files
    # =====================================================
    echo "Downloading FASTQ files..."
    
    local FASTQ_URL="https://cf.10xgenomics.com/samples/cell-vdj/6.0.0/sc5p_v2_hs_B_1k_multi_5gex_b_Multiplex/sc5p_v2_hs_B_1k_multi_5gex_b_Multiplex_fastqs.tar"
    local FASTQ_ARCHIVE="sc5p_v2_hs_B_1k_multi_5gex_b_Multiplex_fastqs.tar"

    curl -LO "${FASTQ_URL}"
    tar -xf "${FASTQ_ARCHIVE}"
    
    # Copy extracted FASTQ folders to standard location
    mkdir -p fastq_samples
    cp -r sc5p_v2_hs_B_1k_multi_5gex_b_fastqs fastq_samples/
    
    # Create symlinks for easier access
    ln -sf sc5p_v2_hs_B_1k_multi_5gex_b_fastqs/sc5p_v2_hs_B_1k_5gex_fastqs fastq_samples/
    ln -sf sc5p_v2_hs_B_1k_multi_5gex_b_fastqs/sc5p_v2_hs_B_1k_b_fastqs fastq_samples/
    
    # Clean up temporary archive
    rm -f "${FASTQ_ARCHIVE}"
    rm -rf sc5p_v2_hs_B_1k_multi_5gex_b_fastqs
    
    echo "FASTQ files downloaded successfully!"
    echo "  - Gene Expression: fastq_samples/sc5p_v2_hs_B_1k_5gex_fastqs"
    echo "  - VDJ/BCR: fastq_samples/sc5p_v2_hs_B_1k_b_fastqs"

    # =====================================================
    # STEP 2: Download Reference Transcriptome (Gene Expression)
    # =====================================================
    echo "Downloading Gene Expression reference..."
    
    local GEX_REF_URL="https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz"
    local GEX_REF_ARCHIVE="refdata-gex-GRCh38-2020-A.tar.gz"

    curl -LO "${GEX_REF_URL}"
    tar -xf "${GEX_REF_ARCHIVE}"
    
    # Copy to references folder
    cp -r refdata-gex-GRCh38-2020-A references/
    
    # Clean up temporary archive
    rm -f "${GEX_REF_ARCHIVE}"
    
    echo "Gene Expression reference downloaded successfully!"
    echo "  - Reference: references/refdata-gex-GRCh38-2020-A"

    # =====================================================
    # STEP 3: Download V(D)J Reference
    # =====================================================
    echo "Downloading V(D)J reference..."
    
    local VDJ_REF_URL="https://cf.10xgenomics.com/supp/cell-vdj/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0.tar.gz"
    local VDJ_REF_ARCHIVE="refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0.tar.gz"

    curl -LO "${VDJ_REF_URL}"
    tar -xf "${VDJ_REF_ARCHIVE}"
    
    # Copy to references folder
    cp -r refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0 references/
    
    # Clean up temporary archive
    rm -f "${VDJ_REF_ARCHIVE}"
    
    echo "V(D)J reference downloaded successfully!"
    echo "  - Reference: references/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0"

    # =====================================================
    # STEP 4: Create CellRanger Multi Config CSV
    # =====================================================
    echo "Creating CellRanger multi configuration CSV..."
    
    cat > cellranger_config.csv << 'CSVEOF'
[gene-expression]
reference,/references/refdata-gex-GRCh38-2020-A
expect-cells,1000
create-bam,true

[vdj]
reference,/references/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0

[libraries]
fastq_id,fastqs,lanes,feature_types,subsample_rate
sc5p_v2_hs_B_1k_5gex,/fastq_samples/sc5p_v2_hs_B_1k_5gex_fastqs,1|2,gene expression,
sc5p_v2_hs_B_1k_b,/fastq_samples/sc5p_v2_hs_B_1k_b_fastqs,1|2,vdj,
CSVEOF

    echo "CellRanger configuration CSV created: cellranger_config.csv"

    # =====================================================
    # Summary of downloaded data
    # =====================================================
    echo ""
    echo "================================================"
    echo "CellRanger Example Data Download Complete"
    echo "================================================"
    echo ""
    echo "FastQ Samples:"
    ls -la fastq_samples/
    echo ""
    echo "References:"
    ls -la references/
    echo ""
    echo "Configuration:"
    cat cellranger_config.csv
    echo ""
    
    # Copy reference config from project reference_data
    if [ -d "${projectDir}/reference_data" ]; then
        cp "${projectDir}/reference_data/cellranger_config.csv" . 2>/dev/null || true
        echo "Reference configuration copied from reference_data folder"
    fi
    """
}

workflow run_cellranger_data {
    take:

    channel FASTQ
    channel GEX_REFERENCE
    channel VDJ_REFERENCE
    channel CONFIG_CSV

    emit:
        ch_fastqs = FASTQ
        ch_gex_reference = GEX_REFERENCE
        ch_vdj_reference = VDJ_REFERENCE
        ch_config = CONFIG_CSV
}
