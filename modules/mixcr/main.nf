process RUN_MIXCR {
    tag "MixCR"
    publishDir "${params.outdir}/mixcr", mode: 'copy'

    input:
    path fasta

    output:
    path "mixcr.tsv", emit: tsv
    path "mixcr_per_cell.tsv", emit: per_cell
    path "*.clns", emit: clns
    path "*.vdjca", emit: vdjca
    path "report.txt", emit: report

    script:
    """
    # Set environment variables to use local writable paths
    export TMPDIR=\$(pwd)
    export HOME=\$(pwd)
    export XDG_CACHE_HOME=\$(pwd)/.cache
    
    # Register MixCR license
    # mixcr activate-license "${params.mixcr_license}"
    # using an environment variable here since activate-license is for interactive use 
    export MI_LICENSE="${params.mixcr_license}"

    # MixCR align
    mixcr align -p 10x-sc-xcr-vdj -s human --report report.txt --force-overwrite $fasta align.vdjca

    # MixCR assemble
    mixcr assemble --report report.txt --force-overwrite align.vdjca clones.clns

    # MixCR export - default clonotype table
    mixcr exportClones --force-overwrite clones.clns mixcr.tsv
    
    # MixCR export - per-cell clonotypes (one row per cell)
    mixcr exportClones \
        --drop-default-fields \
        --split-by-tags Cell \
        -tag cell \
        -cellGroup \
        -uniqueTagCount Molecule \
        -count \
        -vFamily -jFamily \
        -aaFeature CDR3 \
        -nFeatureImputed VDJRegion \
        --force-overwrite \
        clones.clns \
        mixcr_per_cell_raw.tsv
    
    # Rename tagValueCELL column to sequence_id for compatibility with comparison script
    awk 'NR==1 {gsub(/tagValueCELL/, "sequence_id"); print; next} {print}' mixcr_per_cell_raw.tsv > mixcr_per_cell.tsv
    """
}
