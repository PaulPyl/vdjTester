process RUN_MIXCR {
    tag "MixCR"
    storeDir("${params.outdir}/mixcr_store")
    publishDir("${params.outdir}/mixcr", mode: 'copy')

    input:
        path fastq_sample

    output:
        path "mixcr.tsv", emit: tsv
        path "mixcr_per_cell.tsv", emit: per_cell
        path "*.clns", emit: clns
        path "*.vdjca", emit: vdjca
        path "report.txt", emit: report

    script:
    """
    set -eu
    
    # Set environment variables to use local writable paths
    export TMPDIR=$(pwd)
    export HOME=$(pwd)
    export XDG_CACHE_HOME=$(pwd)/.cache
    
    # Register MixCR license
    export MI_LICENSE="${params.mixcr_license}"
    
    echo "Running MixCR pipeline..."
    echo "  Input FASTQ: $fastq_sample"
    
    # Extract sample name from path
    local SAMPLE_NAME=$(basename $fastq_sample)
    
    # Define FASTQ file paths
    local R1_FASTQ="${fastq_sample}/${SAMPLE_NAME}_S1_L001_R1_001.fastq.gz"
    local R2_FASTQ="${fastq_sample}/${SAMPLE_NAME}_S1_L001_R2_001.fastq.gz"
    
    # MixCR align
    echo "Step 1: Running MixCR align..."
    mixcr align -p 10x-sc-xcr-vdj -s human \\
        --report report.txt \\
        --force-overwrite \\
        "${R1_FASTQ}" "${R2_FASTQ}" \\
        align.vdjca
    
    # MixCR assemble
    echo "Step 2: Running MixCR assemble..."
    mixcr assemble --report report.txt \\
        --force-overwrite \\
        align.vdjca \\
        clones.clns
    
    # MixCR export - clonotype table
    echo "Step 3: Exporting clonotype table..."
    mixcr exportClones --force-overwrite \\
        clones.clns \\
        mixcr.tsv
    
    # MixCR export - per-cell clonotypes (one row per cell)
    echo "Step 4: Exporting per-cell clonotypes..."
    mixcr exportClones \\
        --drop-default-fields \\
        --split-by-tags Cell \\
        -tag cell \\
        -cellGroup \\
        -uniqueTagCount Molecule \\
        -count \\
        -vFamily -jFamily \\
        -aaFeature CDR3 \\
        -nFeatureImputed VDJRegion \\
        --force-overwrite \\
        clones.clns \\
        mixcr_per_cell_raw.tsv
    
    # Rename tagValueCELL column to sequence_id for compatibility with comparison script
    awk 'NR==1 {gsub(/tagValueCELL/, "sequence_id"); print; next} {print}' \\
        mixcr_per_cell_raw.tsv \\
        > mixcr_per_cell.tsv
    
    echo "MixCR pipeline completed successfully!"
    """
}

workflow mixcr_workflow {
    input
    channel ch_fastq_input
    
    ch_tsv_out = RUN_MIXCR(ch_fastq_input).tsv
    ch_per_cell_out = RUN_MIXCR(ch_fastq_input).per_cell
    ch_clns_out = RUN_MIXCR(ch_fastq_input).clns
    ch_vdjca_out = RUN_MIXCR(ch_fastq_input).vdjca
    ch_report_out = RUN_MIXCR(ch_fastq_input).report

    output
    channel ch_tsv = ch_tsv_out
    channel ch_per_cell = ch_per_cell_out
    channel ch_clns = ch_clns_out
    channel ch_vdjca = ch_vdjca_out
    channel ch_report = ch_report_out
}
