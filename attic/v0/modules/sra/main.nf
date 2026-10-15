process RUN_SRA_DOWNLOAD {
    tag "SRA-Download"
    storeDir "${params.outdir}/fastq_samples"

    output:
    // emit each sample subdirectory separately so downstream can run per-sample
    path "fastq_samples/*"

    stub:
    """
    # Copy existing fastq samples from reference_data instead of generating new ones
    mkdir -p fastq_samples
    cp -r ${projectDir}/reference_data/fastq_samples/* fastq_samples/
    """

    script:
    """
    # Run the actual download script if it exists
    if [ -f "${projectDir}/bin/run_sra_download.sh" ]; then
        sh ${projectDir}/bin/run_sra_download.sh \
            "${projectDir}/bin/srr_list.txt" \
            ${params.subsample_reads} \
            ${params.subsample_seed}
    fi
    """}