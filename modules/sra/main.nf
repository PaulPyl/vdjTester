process RUN_SRA_DOWNLOAD {
    tag "SRA-Download"
    storeDir "${params.outdir}/fastq_samples"

    output:
    // emit each sample subdirectory separately so downstream can run per-sample
    path "fastq_samples/*"

    stub:
    """
    mkdir -p fastq_samples
    
    # Create a small example fastq file for each sample
    for i in {1..3}; do
        mkdir -p fastq_samples/sample_${i}_sample
        echo "@read1
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
+
III
@read2
TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA
+
III
@read3
CGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT
+
III" > fastq_samples/sample_${i}_sample/sample_${i}_S1_L001_R1_001.fastq.gz
        echo "@read1
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
+
III
@read2
TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA
+
III
@read3
CGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT
+
III" > fastq_samples/sample_${i}_sample/sample_${i}_S1_L001_R2_001.fastq.gz
    done
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
