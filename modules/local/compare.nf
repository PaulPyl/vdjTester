// Harmonise every tool's calls and compare them cell by cell.
process COMPARE_VDJ {
    tag "${sample_id}"
    label 'process_low'
    publishDir "${params.outdir}/compare", mode: 'copy'

    input:
    val sample_id
    path calls

    output:
    path 'vdj_cell_classes.tsv', emit: cell_classes
    path 'vdj_calls_long.tsv'
    path 'vdj_tool_summary.tsv'
    path 'vdj_pairwise_agreement.tsv'

    script:
    """
    vdj_compare.py --outdir . ${calls}
    """

    stub:
    """
    touch vdj_cell_classes.tsv vdj_calls_long.tsv vdj_tool_summary.tsv vdj_pairwise_agreement.tsv
    """
}

// B/T identity from expression, doublet calls, and the cross-tab against VDJ calls.
process EXPRESSION_CHECK {
    tag "${sample_id}"
    label 'process_low'
    publishDir "${params.outdir}/expression", mode: 'copy'

    input:
    tuple val(sample_id), path(filtered_h5)
    path cell_classes
    path markers

    output:
    path 'expression_cells.tsv'
    path 'expression_vs_vdj.tsv'
    path 'expression_summary.tsv'
    path 'marker_scores.png'

    script:
    """
    export NUMBA_CACHE_DIR=\$PWD/.numba MPLCONFIGDIR=\$PWD/.mpl
    expression_check.py \\
        --h5 ${filtered_h5} \\
        --vdj ${cell_classes} \\
        --markers ${markers} \\
        --score-min ${params.score_min} \\
        --min-markers ${params.min_markers} \\
        --outdir .
    """

    stub:
    """
    touch expression_cells.tsv expression_vs_vdj.tsv expression_summary.tsv marker_scores.png
    """
}
