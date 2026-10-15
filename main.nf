#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import module workflows
include { immcantation_workflow  } from './modules/immcantation/main'
include { mixcr_workflow         } from './modules/mixcr/main'
include { run_cellranger    } from './modules/cellranger/main'
include { comparison_workflow    } from './modules/utils/comparison'
include { TEST_LIMA         } from './modules/test_lima/main'
include { run_cellranger_data } from './modules/cellranger_data/main'

// ===========================================================================
// Main Pipeline Workflow - CellRanger Example Data Integration
// ===========================================================================

workflow cellranger_evaluation {
    
    log.info """
    ╔═══════════════════════════════════════════════════════════╗
    ║   V D J   T O O L   E V A L U A T I O N   P I P E L I N E ║
    ╠═══════════════════════════════════════════════════════════╣
    ║   Configuration: CellRanger Example Data Pipeline         ║
    ║   CellRanger Run ID  : ${params.cellranger_run_id}        ║
    ║   Output Directory   : ${params.outdir}                   ║
    ╚═══════════════════════════════════════════════════════════╝
    """

    // ===========================================================================
    // STEP 1: Download CellRanger Example Data
    // ===========================================================================
    log.info "\n[1/5] Downloading CellRanger example data..."
    
    ch_fastqs = DOWNLOAD_CELLRANGER_DATA.ch_fastqs
    ch_gex_reference = DOWNLOAD_CELLRANGER_DATA.ch_gex_reference
    ch_vdj_reference = DOWNLOAD_CELLRANGER_DATA.ch_vdj_reference
    ch_config = DOWNLOAD_CELLRANGER_DATA.ch_config_csv

    // ===========================================================================
    // STEP 2: Run CellRanger Multi Pipeline
    // ===========================================================================
    log.info "[2/5] Running CellRanger Multi analysis..."
    
    // Run CellRanger with all inputs
    cellranger_outputs = RUN_CELLRANGER(
        ch_fastqs,
        ch_gex_reference,
        ch_vdj_reference,
        ch_config
    )
    
    ch_cellranger_gex = cellranger_outputs.ch_gex_outputs
    ch_cellranger_vdj = cellranger_outputs.ch_vdj_outputs
    ch_cellranger_web = cellranger_outputs.ch_web_summary
    ch_cellranger_metrics = cellranger_outputs.ch_metrics

    // ===========================================================================
    // STEP 3: Run Immcantation Pipeline
    // ===========================================================================
    log.info "[3/5] Running Immcantation 10x pipeline..."
    
    ch_input_fasta = channel.fromPath(params.input_fasta)
    ch_input_annotations = channel.fromPath(params.input_annotations)
    
    // Run Immcantation with CellRanger VDJ output
    ch_immcantation_tsv = immcantation_workflow(ch_cellranger_vdj, ch_input_annotations).ch_tsv
    ch_immcantation_db_pass = immcantation_workflow(ch_cellranger_vdj, ch_input_annotations).ch_db_pass

    // ===========================================================================
    // STEP 4: Run MixCR Pipeline
    // ===========================================================================
    log.info "[4/5] Running MixCR pipeline..."
    
    // Run MixCR with VDJ FASTQ data
    ch_mixcr_tsv = mixcr_workflow(ch_fastqs).ch_tsv
    ch_mixcr_per_cell = mixcr_workflow(ch_fastqs).ch_per_cell
    ch_mixcr_clns = mixcr_workflow(ch_fastqs).ch_clns
    ch_mixcr_vdjca = mixcr_workflow(ch_fastqs).ch_vdjca

    // ===========================================================================
    // STEP 5: Compare and Generate Report
    // ===========================================================================
    log.info "[5/5] Comparing results and generating report..."
    
    ch_reference = channel.fromPath(params.reference_data)
    
    // Collect all pipeline TSV outputs for comparison
    ch_pipeline_tsvs = channel.empty()
        .mix(ch_immcantation_tsv, ch_mixcr_tsv)
        .collect()
    
    // Run comparison analysis
    ch_comparison_summary = comparison_workflow(ch_pipeline_tsvs, ch_reference).ch_summary_csv
    ch_comparison_plots = comparison_workflow(ch_pipeline_tsvs, ch_reference).ch_plots
    ch_comparison_report = comparison_workflow(ch_pipeline_tsvs, ch_reference).ch_report

    // ===========================================================================
    // Pipeline Complete - Final Summary
    // ===========================================================================
    log.info """
    ╔═══════════════════════════════════════════════════════════╗
    ║           PIPELINE EXECUTION COMPLETE                      ║
    ╚═══════════════════════════════════════════════════════════╝
    
    Output Directory: ${params.outdir}
    
    ┌──────────────────────────────────────────────────────────────┐
    │ Results Summary                                                  │
    ├──────────────────────────────────────────────────────────────┤
    │ CellRanger   : ${params.outdir}/cellranger_output/         │
    │ Immcantation : ${params.outdir}/immcantation/            │
    │ MixCR        : ${params.outdir}/mixcr/                    │
    │ Comparison   : ${params.outdir}/comparison/               │
    └──────────────────────────────────────────────────────────────┘
    
    Key Deliverables:
    • comparison_summary.csv  - Pipeline accuracy comparison
    • *.png files            - Visualization outputs
    • comparison_report.html - Interactive HTML report
    """
}

// ===========================================================================
// Test Workflows
// ===========================================================================

workflow test_cellranger_data {
    // Test CellRanger data download module only
    log.info "Testing CellRanger data download..."
    DOWNLOAD_CELLRANGER_DATA()
}

workflow test_cellranger {
    // Test CellRanger multi pipeline with downloaded data
    log.info "Testing CellRanger multi pipeline..."
    
    ch_fastqs = DOWNLOAD_CELLRANGER_DATA.ch_fastqs
    ch_gex_ref = DOWNLOAD_CELLRANGER_DATA.ch_gex_reference
    ch_vdj_ref = DOWNLOAD_CELLRANGER_DATA.ch_vdj_reference
    ch_config = DOWNLOAD_CELLRANGER_DATA.ch_config_csv
    
    RUN_CELLRANGER(ch_fastqs, ch_gex_ref, ch_vdj_ref, ch_config)
}

workflow test_immcantation {
    // Test Immcantation with CellRanger VDJ output
    log.info "Testing Immcantation pipeline..."
    
    ch_fastqs = DOWNLOAD_CELLRANGER_DATA.ch_fastqs
    ch_gex_ref = DOWNLOAD_CELLRANGER_DATA.ch_gex_reference
    ch_vdj_ref = DOWNLOAD_CELLRANGER_DATA.ch_vdj_reference
    ch_config = DOWNLOAD_CELLRANGER_DATA.ch_config_csv
    
    ch_vdj_data = RUN_CELLRANGER(ch_fastqs, ch_gex_ref, ch_vdj_ref, ch_config).ch_vdj_outputs
    
    ch_input_fasta = channel.fromPath(params.input_fasta)
    ch_input_annotations = channel.fromPath(params.input_annotations)
    
    RUN_IMMCANTATION(ch_vdj_data, ch_input_annotations)
}

workflow test_mixcr {
    // Test MixCR with CellRanger VDJ FASTQ data
    log.info "Testing MixCR pipeline..."
    
    ch_fastqs = DOWNLOAD_CELLRANGER_DATA.ch_fastqs
    RUN_MIXCR(ch_fastqs)
}

workflow test_comparison {
    // Test comparison module with all pipeline outputs
    log.info "Testing comparison analysis..."
    
    ch_fastqs = DOWNLOAD_CELLRANGER_DATA.ch_fastqs
    ch_gex_ref = DOWNLOAD_CELLRANGER_DATA.ch_gex_reference
    ch_vdj_ref = DOWNLOAD_CELLRANGER_DATA.ch_vdj_reference
    ch_config = DOWNLOAD_CELLRANGER_DATA.ch_config_csv
    
    ch_vdj_data = RUN_CELLRANGER(ch_fastqs, ch_gex_ref, ch_vdj_ref, ch_config).ch_vdj_outputs
    
    ch_input_fasta = channel.fromPath(params.input_fasta)
    ch_input_annotations = channel.fromPath(params.input_annotations)
    
    ch_immcantation_tsv = RUN_IMMCANTATION(ch_vdj_data, ch_input_annotations).ch_tsv
    ch_mixcr_tsv = RUN_MIXCR(ch_fastqs).ch_tsv
    ch_reference = channel.fromPath(params.reference_data)
    
    ch_pipeline_tsvs = channel.empty()
        .mix(ch_immcantation_tsv, ch_mixcr_tsv)
        .collect()
    
    RUN_COMPARISON(ch_pipeline_tsvs, ch_reference)
}

workflow test_lima {
    // Test LIMA infrastructure
    log.info "Testing LIMA infrastructure..."
    TEST_LIMA()
}
