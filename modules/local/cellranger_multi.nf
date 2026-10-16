// Cell Ranger multi on one sample. Runs from a host install (params.cellranger_dir),
// not a container: 10x does not allow redistributing Cell Ranger.
//
// Each library arrives as its own FASTQ directory (staged under fastqs/), and the
// multi config is written here from the library list so its paths point into the
// task directory; a hand-written CSV with absolute paths breaks on every host.
process CELLRANGER_MULTI {
    tag "${sample_id}"
    label 'process_cellranger'
    publishDir path: { "${params.outdir}/cellranger/${sample_id}" }, mode: 'copy', pattern: '{multi_config.csv,versions.txt,per_sample_outs/**}'

    input:
    tuple val(sample_id), val(libraries), path(lib_dirs, stageAs: 'fastqs/*'), path(gex_ref), path(vdj_ref)

    output:
    tuple val(sample_id), path('per_sample_outs/vdj_b'), emit: vdj_b, optional: true
    tuple val(sample_id), path('per_sample_outs/vdj_t'), emit: vdj_t, optional: true
    tuple val(sample_id), path('per_sample_outs/sample_filtered_feature_bc_matrix.h5'), emit: filtered_h5
    tuple val(sample_id), path("cellranger.${sample_id}.*.csv"), emit: calls
    path 'per_sample_outs/{web_summary.html,metrics_summary.csv}', emit: summary
    tuple val(sample_id), path('per_sample_outs/metrics_summary.csv'), emit: metrics
    path 'multi_config.csv'
    path 'versions.txt', emit: versions

    script:
    def cr_path = params.cellranger_dir ? "export PATH=${params.cellranger_dir}:\$PATH" : ''
    def lib_lines = libraries.collect { lib ->
        "printf '%s\\n' \"${lib.fastq_id},\$(readlink -f fastqs/${lib.dir}),${lib.feature_types}\""
    }.join('\n        ')
    """
    ${cr_path}
    # Cell Ranger writes telemetry under \$HOME; containers may give it an unwritable one.
    [ -w "\${HOME:-/nonexistent}" ] || export HOME=\$PWD
    cellranger --version > versions.txt

    {
        printf '%s\\n' '[gene-expression]' "reference,\$(readlink -f ${gex_ref})" 'create-bam,false' ''
        printf '%s\\n' '[vdj]' "reference,\$(readlink -f ${vdj_ref})" ''
        printf '%s\\n' '[libraries]' 'fastq_id,fastqs,feature_types'
        ${lib_lines}
    } > multi_config.csv

    cellranger multi \\
        --id=${sample_id} \\
        --csv=multi_config.csv \\
        --localcores=${task.cpus} \\
        --localmem=${task.memory.toGiga()} \\
        --disable-ui

    # Without a [samples] section there is exactly one per-sample directory.
    mv ${sample_id}/outs/per_sample_outs/* per_sample_outs_tmp
    mv per_sample_outs_tmp per_sample_outs
    rm -rf ${sample_id}/SC_MULTI_CS

    for chain in vdj_b vdj_t; do
        if [ -f per_sample_outs/\$chain/filtered_contig_annotations.csv ]; then
            cp per_sample_outs/\$chain/filtered_contig_annotations.csv cellranger.${sample_id}.\$chain.csv
        fi
    done
    """

    stub:
    def chains = libraries.findAll { lib -> lib.feature_types.startsWith('VDJ') }
        .collect { lib -> lib.feature_types == 'VDJ-B' ? 'vdj_b' : 'vdj_t' }
        .join(' ')
    """
    echo "cellranger stub" > versions.txt
    touch multi_config.csv
    mkdir -p per_sample_outs
    touch per_sample_outs/sample_filtered_feature_bc_matrix.h5 per_sample_outs/web_summary.html per_sample_outs/metrics_summary.csv
    for chain in ${chains}; do
        mkdir -p per_sample_outs/\$chain
        touch per_sample_outs/\$chain/filtered_contig.fasta per_sample_outs/\$chain/filtered_contig_annotations.csv
        touch cellranger.${sample_id}.\$chain.csv
    done
    """
}
