process COMPARE_RESULTS {
    tag "Comparison"
    publishDir "${params.outdir}/comparison", mode: 'copy'

    input:
    path pipeline_tsvs
    path reference_tsv

    output:
    path "comparison_summary.csv"
    path "*.png"
    path "comparison_report.html"

    script:
    """
    # Run the detailed comparison analysis
    compare_results_detailed.R ${pipeline_tsvs.join(' ')} $reference_tsv .
    
    # Copy Quarto template and render HTML report
    cp ${projectDir}/bin/comparison_report.qmd .
    quarto render comparison_report.qmd --to html
    """
}
