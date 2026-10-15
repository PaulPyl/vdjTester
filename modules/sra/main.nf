process RUN_SRA_DOWNLOAD {
    tag "SRA-Download"

    // Store outputs in dedicated directory for efficient cache management
    storeDir "${params.outdir}"

    input: path sraSamples

    output: path "fastq_samples/*"

    script:
    """
    # Run the actual download script if it exists
    if [ -f "run_sra_download.sh" ]; then
        run_sra_download.sh \
            "bin/srr_list.txt" \
            \${params.subsample_reads} \
            \${params.subsample_seed}
    fi
    """

    stub:
    """
    # Check for test data directory or create stub outputs
    TEST_DATA_DIR="${task.configDir}/test_data/fastq_samples/SRR11233622"
    
    if [ -d "\$TEST_DATA_DIR" ]; then
        echo "Using test data from \$TEST_DATA_DIR"
        mkdir -p fastq_samples
        cp -r "\$TEST_DATA_DIR"/* fastq_samples/
        echo "SRR11233622 stub outputs copied to fastq_samples/"
    else
        echo "Test data not found at \$TEST_DATA_DIR"
        # Create fallback stub outputs
        mkdir -p fastq_samples/SRR11233622
        echo "Created stub outputs for SRR11233622"
    fi
    
    # Copy supporting reference data files
    cp -r ${projectDir}/reference_data/*.tsv reference_data/
    cp -r ${projectDir}/reference_data/*.csv reference_data/
    cp -r ${projectDir}/reference_data/*.rds reference_data/
    cp -r ${projectDir}/reference_data/filtered_contig.fasta reference_data/
    
    echo "SRA Download stub completed successfully"
    """
}
