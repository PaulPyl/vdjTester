process RUN_COMPARISON {
    tag "Comparison"
    storeDir("${params.outdir}/comparison_store")
    publishDir("${params.outdir}/comparison", mode: 'copy')

    input:
        path pipeline_tsvs
        path reference_tsv

    output:
        path "comparison_summary.csv", emit: summary_csv
        path "*.png", emit: plots
        path "comparison_report.html", emit: report

    script:
    """
    set -eu
    
    echo "Running comparison analysis..."
    echo "  Pipeline outputs: ${pipeline_tsvs.size} files"
    echo "  Reference: $reference_tsv"
    
    # Run the detailed comparison analysis
    ${projectDir}/bin/compare_results_detailed.R ${pipeline_tsvs.join(' ')} $reference_tsv .
    
    # Copy Quarto template and render HTML report
    cp ${projectDir}/bin/comparison_report.qmd .
    quarto render comparison_report.qmd --to html
    
    echo "Comparison analysis completed successfully!"
    echo "  Summary: comparison_summary.csv"
    echo "  Plots: *.png"
    echo "  Report: comparison_report.html"
    """
}

workflow comparison_workflow {
    input
    channel ch_pipeline_tsvs, ch_reference
    
    ch_summary_csv_out = RUN_COMPARISON(ch_pipeline_tsvs, ch_reference).summary_csv
    ch_plots_out = RUN_COMPARISON(ch_pipeline_tsvs, ch_reference).plots
    ch_report_out = RUN_COMPARISON(ch_pipeline_tsvs, ch_reference).report

    output
    channel ch_summary_csv = ch_summary_csv_out
    channel ch_plots = ch_plots_out
    channel ch_report = ch_report_out
}
