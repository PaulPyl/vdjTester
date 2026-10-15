process RUN_IMMCANTATION {
    tag "Immcantation 10x"
    publishDir "${params.outdir}/immcantation", mode: 'copy'

    input:
    path fasta
    path annotations

    output:
    path "immcantation.tsv", emit: tsv
    path "*_db-pass.tsv", emit: db_pass

    script:
    """
    changeo-10x -s $fasta -a $annotations -o . -g human -t ig -x 0.16
    
    # Rename output to standardized pipeline name
    mv *_germ-pass.tsv immcantation.tsv
    """
}
