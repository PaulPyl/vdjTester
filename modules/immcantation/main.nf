process RUN_IMMCANTATION {
    tag "Immcantation 10x"
    storeDir("${params.outdir}/immcantation_store")
    publishDir("${params.outdir}/immcantation", mode: 'copy')

    input:
        path vdj_data
        path annotations

    output:
        path "immcantation.tsv", emit: tsv
        path "*_db-pass.tsv", emit: db_pass

    script:
    """
    set -eu
    
    echo "Running Immcantation 10x pipeline..."
    echo "  VDJ Input: $vdj_data"
    echo "  Annotations: $annotations"
    
    # Run changeo-10x for immune repertoire analysis
    changeo-10x -s $vdj_data -a $annotations -o . -g human -t ig -x 0.16
    
    # Rename output to standardized pipeline name
    mv *_germ-pass.tsv immcantation.tsv
    
    echo "Immcantation pipeline completed successfully!"
    """
}

workflow immcantation_workflow {
    input
    channel ch_vdj_input, ch_annotations
    
    ch_tsv_out = RUN_IMMCANTATION(ch_vdj_input, ch_annotations).tsv
    ch_db_pass_out = RUN_IMMCANTATION(ch_vdj_input, ch_annotations).db_pass

    output
    channel ch_tsv = ch_tsv_out
    channel ch_db_pass = ch_db_pass_out
}
