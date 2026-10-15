process RUN_CELLRANGER {
    tag "CellRanger Multi"
    storeDir("${params.outdir}/cellranger_store")
    publishDir("${params.outdir}/cellranger_output", mode: 'copy')

    input:
        path fastqs
        path gex_reference
        path vdj_reference
        path config_csv

    output:
        path "*/outs/*", emit: all_outputs
        path "*/outs/multi/count/*", emit: gex_outputs
        path "*/outs/multi/vdj_*/*", emit: vdj_outputs
        path "*/outs/per_sample_outs/*/web_summary.html", emit: web_summary
        path "*/outs/per_sample_outs/*/metrics_summary.csv", emit: metrics
        path "_cmdline", emit: cmdline
        path "cr_version.txt", emit: version

    script:
    """
    set -eu

    echo "========================================"
    echo "Running CellRanger Multi Pipeline"
    echo "========================================"
    echo ""
    echo "Inputs:"
    echo "  FASTQs: ${fastqs}"
    echo "  GEX Reference: ${gex_reference}"
    echo "  VDJ Reference: ${vdj_reference}"
    echo "  Config: ${config_csv}"
    echo ""

    # Create output directory structure
    RUN_ID="$params.cellranger_run_id"
    if [ -z "$RUN_ID" ]; then
        RUN_ID="HumanB_Cell_multi"
    fi
    mkdir -p "$RUN_ID"
    cd "$RUN_ID"

    # Copy configuration file and update paths if needed
    cp ${config_csv} .

    # Check if Cell Ranger is available
    if command -v cellranger &> /dev/null; then
        CELLRANGER_CMD="cellranger"
    elif [ -x /opt/cellranger-10.0.0/cellranger ]; then
        CELLRANGER_CMD="/opt/cellranger-10.0.0/cellranger"
    else
        echo "ERROR: cellranger binary not found in PATH or /opt" >&2
        exit 1
    fi

    # Get and save Cell Ranger version
    ${CELLRANGER_CMD} --version > ../cr_version.txt
    echo "Cell Ranger version saved"

    # Run Cell Ranger Multi pipeline
    echo "Starting Cell Ranger Multi analysis..."
    
    ${CELLRANGER_CMD} multi \\
        --id="${RUN_ID}" \\
        --csv="./cellranger_config.csv" \\
        --jobmode=local \\
        --localcores=${task.cpus} \\
        --localmem=${task.memory.toHuman().replace('GB', 'g').replace('MB', 'm')}

    echo "Cell Ranger Multi pipeline completed successfully!"
    echo "Output directory: ${RUN_ID}/"

    # List output structure
    echo ""
    echo "Output structure:"
    tree -L 2 outs/ 2>/dev/null || find outs -maxdepth 2 -type d | head -30

    # Save command line for reference
    echo "${CELLRANGER_CMD} multi --id=${RUN_ID} --csv=./cellranger_config.csv" > _cmdline
    """

    stub:
    """
    echo "CellRanger Multi stub execution"
    mkdir -p ${params.cellranger_run_id:-HumanB_Cell_multi}/outs
    mkdir -p ${params.cellranger_run_id:-HumanB_Cell_multi}/outs/multi/count
    mkdir -p ${params.cellranger_run_id:-HumanB_Cell_multi}/outs/multi/vdj_b
    mkdir -p ${params.cellranger_run_id:-HumanB_Cell_multi}/outs/per_sample_outs
    
    # Create stub outputs
    echo "Stub FASTQ received: ${fastqs}" > cr_version.txt
    echo "Stub GEX Reference: ${gex_reference}" >> cr_version.txt
    echo "Stub VDJ Reference: ${vdj_reference}" >> cr_version.txt
    echo "Stub Config: ${config_csv}" >> cr_version.txt
    
    # Create stub output files
    mkdir -p ${params.cellranger_run_id:-HumanB_Cell_multi}/outs/multi/count/raw_feature_bc_matrix
    touch ${params.cellranger_run_id:-HumanB_Cell_multi}/outs/multi/count/raw_feature_bc_matrix/matrix.h5
    touch ${params.cellranger_run_id:-HumanB_Cell_multi}/outs/multi/count/web_summary.html
    
    mkdir -p ${params.cellranger_run_id:-HumanB_Cell_multi}/outs/multi/vdj_b
    touch ${params.cellranger_run_id:-HumanB_Cell_multi}/outs/multi/vdj_b/all_contig_annotations.csv
    touch ${params.cellranger_run_id:-HumanB_Cell_multi}/outs/multi/vdj_b/all_contig.fasta
    touch ${params.cellranger_run_id:-HumanB_Cell_multi}/outs/multi/vdj_b/all_contig.fastq
    
    mkdir -p ${params.cellranger_run_id:-HumanB_Cell_multi}/outs/per_sample_outs
    touch ${params.cellranger_run_id:-HumanB_Cell_multi}/outs/per_sample_outs/web_summary.html
    touch ${params.cellranger_run_id:-HumanB_Cell_multi}/outs/per_sample_outs/metrics_summary.csv
    
    echo "${params.cellranger_run_id:-HumanB_Cell_multi} multi --id=${params.cellranger_run_id:-HumanB_Cell_multi} --csv=./cellranger_config.csv" > _cmdline
    """
}

workflow run_cellranger {
    input:
        ch_fastqs
        ch_gex_reference
        ch_vdj_reference
        ch_config

    output:
        ch_all_outputs
        ch_gex_outputs
        ch_vdj_outputs
        ch_web_summary
        ch_metrics
        ch_cmdline
        ch_version
}
