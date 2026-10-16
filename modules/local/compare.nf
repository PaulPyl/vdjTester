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
    tuple val(sample_id), path("${sample_id}.vdj_cell_classes.tsv"), emit: cell_classes
    tuple val(sample_id), path("${sample_id}.vdj_calls_long.tsv"), emit: calls_long
    path "${sample_id}.vdj_{tool_summary,pairwise_agreement,v_disagreements,cell_classes}.tsv", emit: report

    script:
    """
    vdj_compare.py --gex-h5 ${filtered_h5} --vdj-ref ${vdj_ref} \\
        --germline-fasta ${germline_ref}/human_IMGT+C.fa \\
        --gex-libraries '${gex_libraries.join(',')}' \\
        --trust4-min-umis ${params.trust4_min_umis} --trust4-min-frac ${params.trust4_min_frac} \\
        --trust4-gex-min-umis ${params.trust4_gex_min_umis} \\
        --outdir . ${calls}
    for f in vdj_*.tsv; do mv "\$f" "${sample_id}.\$f"; done
    """

    stub:
    """
    for f in vdj_cell_classes vdj_calls_long vdj_tool_summary vdj_pairwise_agreement vdj_v_disagreements; do
        touch ${sample_id}.\$f.tsv
    done
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
    path "${sample_id}.{expression_cells,expression_vs_vdj,expression_summary}.tsv", emit: report
    path "${sample_id}.marker_scores.png"

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
    for f in expression_*.tsv marker_scores.png; do mv "\$f" "${sample_id}.\$f"; done
    """

    stub:
    """
    touch ${sample_id}.expression_cells.tsv ${sample_id}.expression_vs_vdj.tsv ${sample_id}.expression_summary.tsv ${sample_id}.marker_scores.png
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
    path "${patient}.{clonality,top_clones,clone_tracking,timepoint_overlap,dominant_clone_agreement}.tsv", emit: report

    script:
    def samples = [sample_ids, timepoints, calls].transpose()
        .collect { sid, tp, f -> "--sample ${sid} ${tp} ${f}" }
        .join(' ')
    """
    clone_tracking.py ${samples} --outdir .
    for f in clonality top_clones clone_tracking timepoint_overlap dominant_clone_agreement; do
        mv \$f.tsv ${patient}.\$f.tsv
    done
    """

    stub:
    """
    for f in clonality top_clones clone_tracking timepoint_overlap dominant_clone_agreement; do
        touch ${patient}.\$f.tsv
    done
    """
}

// One self-contained HTML report for the whole run, from every step's tables.
process REPORT {
    label 'process_low'
    publishDir params.outdir, mode: 'copy'

    input:
    path data, stageAs: 'data/*'
    path qmd

    output:
    path 'vdjtester_report.html'

    script:
    """
    export HOME=\$PWD XDG_CACHE_HOME=\$PWD/.cache MPLCONFIGDIR=\$PWD/.mpl NUMBA_CACHE_DIR=\$PWD/.numba
    # Quarto renders next to the source, so work on a copy in the task directory.
    cp ${qmd} vdjtester_report.qmd
    quarto render vdjtester_report.qmd --to html -P data_dir:\$PWD/data
    """

    stub:
    """
    touch vdjtester_report.html
    """
}
