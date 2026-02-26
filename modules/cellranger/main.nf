process RUN_CELLRANGER {
    tag "CellRanger"
    publishDir "${params.outdir}/cellranger", mode: 'copy'

    input:
    path fastq_dirs

    output:
    path "cr_version.txt"

    script:
    """
    set -eu
    echo "Received fastq dir: $fastq_dirs"

    sample_id=\${fastq_dirs##*/}

    # delegate to helper script
    sh ${projectDir}/bin/run_cellranger.sh "$fastq_dirs" "$sample_id"
    """
}
