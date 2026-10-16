// Immcantation's 10x workflow (changeo-10x): IgBLAST re-annotation of Cell Ranger's
// assembled contigs, then Change-O MakeDb into AIRR format. It does not assemble
// from raw reads, so it is a second annotation of Cell Ranger's contigs, not an
// independent caller. Annotation is Python (IgBLAST + Change-O); clone definition
// (scoper) and the threshold plot run R inside the container.
process IMMCANTATION {
    tag "${sample_id}:${chain}"
    label 'process_low'
    publishDir path: { "${params.outdir}/immcantation/${sample_id}" }, mode: 'copy'

    input:
    tuple val(sample_id), val(chain), path(vdj_dir)

    output:
    tuple val(sample_id), path("immcantation.${sample_id}.${chain}.tsv"), emit: calls
    path "immcantation.${sample_id}.${chain}.log"

    script:
    // changeo-10x takes `ig` or `tr` for the receptor type.
    def receptor = chain == 'vdj_b' ? 'ig' : 'tr'
    """
    changeo-10x \\
        -s ${vdj_dir}/filtered_contig.fasta \\
        -a ${vdj_dir}/filtered_contig_annotations.csv \\
        -o out \\
        -n ${sample_id}_${chain} \\
        -g human \\
        -t ${receptor} \\
        -x 0.16 \\
        -p ${task.cpus} \\
        > immcantation.${sample_id}.${chain}.log 2>&1

    # MakeDb output: every annotated contig, before clonal grouping.
    cp out/${sample_id}_${chain}_db-pass.tsv immcantation.${sample_id}.${chain}.tsv
    """

    stub:
    """
    touch immcantation.${sample_id}.${chain}.tsv immcantation.${sample_id}.${chain}.log
    """
}
