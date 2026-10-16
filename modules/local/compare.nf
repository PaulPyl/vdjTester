// Harmonise every tool's calls and compare them cell by cell.
process COMPARE_VDJ {
    tag "${sample_id}"
    label 'process_low'
    publishDir path: { "${params.outdir}/compare/${sample_id}" }, mode: 'copy'

    input:
    tuple val(sample_id), path(calls), path(filtered_h5), val(gex_libraries)
    path vdj_ref
    path germline_ref

    output:
    tuple val(sample_id), path('vdj_cell_classes.tsv'), emit: cell_classes
    tuple val(sample_id), path('vdj_calls_long.tsv'), emit: calls_long
    path 'vdj_tool_summary.tsv'
    path 'vdj_pairwise_agreement.tsv'
    path 'vdj_v_disagreements.tsv'

    script:
    """
    vdj_compare.py --gex-h5 ${filtered_h5} --vdj-ref ${vdj_ref} \\
        --germline-fasta ${germline_ref}/human_IMGT+C.fa \\
        --gex-libraries '${gex_libraries.join(',')}' \\
        --trust4-min-umis ${params.trust4_min_umis} --trust4-min-frac ${params.trust4_min_frac} \\
        --trust4-gex-min-umis ${params.trust4_gex_min_umis} \\
        --outdir . ${calls}
    """

    stub:
    """
    touch vdj_cell_classes.tsv vdj_calls_long.tsv vdj_tool_summary.tsv vdj_pairwise_agreement.tsv vdj_v_disagreements.tsv
    """
}

// B/T identity from expression, doublet calls, and the cross-tab against VDJ calls.
process EXPRESSION_CHECK {
    tag "${sample_id}"
    label 'process_low'
    publishDir path: { "${params.outdir}/expression/${sample_id}" }, mode: 'copy'

    input:
    tuple val(sample_id), path(filtered_h5), path(cell_classes)
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
        --score-method ${params.score_method} \\
        ${params.max_rank ? "--max-rank ${params.max_rank}" : ''} \\
        ${params.score_min != null ? "--score-min ${params.score_min}" : ''} \\
        --min-markers ${params.min_markers} \\
        --outdir .
    """

    stub:
    """
    touch expression_cells.tsv expression_vs_vdj.tsv expression_summary.tsv marker_scores.png
    """
}

// Clone sizes per timepoint and how they change over time, per tool, for one patient.
process CLONE_TRACKING {
    tag "${patient}"
    label 'process_low'
    publishDir path: { "${params.outdir}/clones/${patient}" }, mode: 'copy'

    input:
    tuple val(patient), val(sample_ids), val(timepoints), path(calls, stageAs: 'calls_?/*')

    output:
    path 'clonality.tsv'
    path 'top_clones.tsv'
    path 'clone_tracking.tsv'
    path 'timepoint_overlap.tsv'
    path 'dominant_clone_agreement.tsv'

    script:
    def samples = [sample_ids, timepoints, calls].transpose()
        .collect { sid, tp, f -> "--sample ${sid} ${tp} ${f}" }
        .join(' ')
    """
    clone_tracking.py ${samples} --outdir .
    """

    stub:
    """
    touch clonality.tsv top_clones.tsv clone_tracking.tsv timepoint_overlap.tsv dominant_clone_agreement.tsv
    """
}
