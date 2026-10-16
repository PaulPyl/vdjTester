// MiXCR on the raw reads of one VDJ library: assembles independently of Cell Ranger.
// Exports one row per clone and cell from the preset's final, cell-assembled
// clonotypes (<prefix>.assembledCells.clns).
process MIXCR {
    tag "${sample_id}:${fastq_id}"
    label 'process_medium'
    publishDir path: { "${params.outdir}/mixcr/${sample_id}" }, mode: 'copy'
    secret 'MI_LICENSE'

    input:
    tuple val(sample_id), val(fastq_id), path(lib_dir)

    output:
    tuple val(sample_id), path("mixcr.${sample_id}.${fastq_id}.tsv"), emit: calls
    path "${fastq_id}.assembledCells.clns"
    path "${fastq_id}.*report*"

    script:
    """
    export HOME=\$PWD TMPDIR=\$PWD XDG_CACHE_HOME=\$PWD/.cache

    # Concatenate lanes; gzip members concatenate cleanly.
    cat ${lib_dir}/${fastq_id}_S*_L*_R1_001.fastq.gz > R1.fastq.gz
    cat ${lib_dir}/${fastq_id}_S*_L*_R2_001.fastq.gz > R2.fastq.gz

    mixcr analyze 10x-sc-xcr-vdj \\
        --species hs \\
        --threads ${task.cpus} \\
        R1.fastq.gz R2.fastq.gz \\
        ${fastq_id}

    mixcr exportClones \\
        --split-by-tags Cell \\
        --drop-default-fields \\
        -tag CELL -chains -readCount -uniqueTagCount Molecule \\
        -vGene -dGene -jGene -cGene \\
        -nFeature CDR3 -aaFeature CDR3 \\
        ${fastq_id}.assembledCells.clns \\
        mixcr.${sample_id}.${fastq_id}.tsv

    rm R1.fastq.gz R2.fastq.gz
    """

    stub:
    """
    touch mixcr.${sample_id}.${fastq_id}.tsv ${fastq_id}.assembledCells.clns ${fastq_id}.report.txt
    """
}
