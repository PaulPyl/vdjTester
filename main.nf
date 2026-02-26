#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import modules
include { RUN_IMMCANTATION } from './modules/immcantation/main'
include { RUN_MIXCR        } from './modules/mixcr/main'
include { RUN_CELLRANGER   } from './modules/cellranger/main'
include { RUN_SRA_DOWNLOAD } from './modules/sra/main'
include { COMPARE_RESULTS  } from './modules/utils/comparison'

workflow test_cellranger {
    RUN_CELLRANGER()
}

workflow test_immcantation {
    ch_input_fasta = Channel.fromPath(params.input_fasta)
    ch_input_annotations = Channel.fromPath(params.input_annotations)
    RUN_IMMCANTATION(ch_input_fasta, ch_input_annotations)
}

workflow test_mixcr {
    ch_input_fasta = Channel.fromPath(params.input_fasta)
    RUN_MIXCR(ch_input_fasta)
}

workflow test_sra_download {
    // create tiny subsample to avoid heavy download
    RUN_SRA_DOWNLOAD()
}

workflow {
    log.info """
    ================================================================
    V D J   T O O L   E V A L U A T I O N   P I P E L I N E
    ================================================================
    Input Fasta       : ${params.input_fasta}
    Input Annotations : ${params.input_annotations}
    Reference Data    : ${params.reference_data}
    Output Directory  : ${params.outdir}
    ================================================================
    """

    // Channel for input files
    ch_input_fasta = Channel.fromPath(params.input_fasta)
    ch_input_annotations = Channel.fromPath(params.input_annotations)
    ch_reference = Channel.fromPath(params.reference_data)

    // Download / subsample SRA fastqs and then run CellRanger (version check)
    ch_fastqs = RUN_SRA_DOWNLOAD()
    RUN_CELLRANGER(ch_fastqs)

    // Run Immcantation 10x pipeline
    RUN_IMMCANTATION(ch_input_fasta, ch_input_annotations)

    // Run MixCR pipeline
    RUN_MIXCR(ch_input_fasta)

    // Collect all pipeline results into a single channel
    // Easy to extend: just add another .mix(NEW_PIPELINE.out.tsv) below
    ch_pipeline_tsvs = Channel.empty()
        .mix(
            RUN_IMMCANTATION.out.tsv,
            RUN_MIXCR.out.per_cell
        )
        .collect()

    // Compare results
    COMPARE_RESULTS(ch_pipeline_tsvs, ch_reference)
}
