process RUN_SRA_DOWNLOAD {
    tag "SRA-Download"
    publishDir "${params.outdir}/fastq_samples", mode: 'copy'

    output:
    // emit each sample subdirectory separately so downstream can run per-sample
    path "fastq_samples/*"

    script:
    """
    sh ${projectDir}/bin/run_sra_download.sh \
        "${projectDir}/bin/srr_list.txt" \
        ${params.subsample_reads} \
        ${params.subsample_seed}
    """}
