// TRUST4 on the raw reads of one VDJ library: a second de novo assembler, besides
// MiXCR, that does not depend on Cell Ranger's assembly.
//
// 10x 5': the cDNA read (R2) is assembled; barcode (0-15) and UMI (16-25) come from
// R1. TRUST4 has no cell calling, so its barcodes are restricted to Cell Ranger's
// GEX cells later, in vdj_compare.py.
process TRUST4 {
    tag "${sample_id}:${fastq_id}"
    label 'process_medium'
    publishDir path: { "${params.outdir}/trust4/${sample_id}" }, mode: 'copy'

    input:
    tuple val(sample_id), val(fastq_id), path(lib_dir), path(trust4_ref)

    output:
    tuple val(sample_id), path("trust4.${sample_id}.${fastq_id}.tsv"), emit: calls
    path "${fastq_id}_report.tsv"

    script:
    """
    cat ${lib_dir}/${fastq_id}_S*_L*_R1_001.fastq.gz > R1.fastq.gz
    cat ${lib_dir}/${fastq_id}_S*_L*_R2_001.fastq.gz > R2.fastq.gz

    run-trust4 \\
        -t ${task.cpus} \\
        -f ${trust4_ref}/hg38_bcrtcr.fa \\
        --ref ${trust4_ref}/human_IMGT+C.fa \\
        -u R2.fastq.gz \\
        --barcode R1.fastq.gz \\
        --UMI R1.fastq.gz \\
        --readFormat bc:0:15,um:16:25 \\
        -o ${fastq_id} \\
        --od out

    cp out/${fastq_id}_barcode_airr.tsv trust4.${sample_id}.${fastq_id}.tsv
    cp out/${fastq_id}_report.tsv .
    rm R1.fastq.gz R2.fastq.gz
    """

    stub:
    """
    touch trust4.${sample_id}.${fastq_id}.tsv ${fastq_id}_report.tsv
    """
}
