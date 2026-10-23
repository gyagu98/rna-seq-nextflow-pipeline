process DESEQ2 {
    label 'process_low'
    publishDir "${params.outdir}/deseq2", mode: 'copy'
    conda "${projectDir}/environment.yml"
    input:
    path quant_dirs
    path samplesheet
    path tx2gene
    val reference
    val comparison
    output:
    path "differential_expression_summary.csv", emit: dge_summary
    path "normalized_counts.csv"
    path "*.pdf"
    path "analysis_session.txt"
    script:
    """
    Rscript ${projectDir}/bin/deseq2_analysis.R \
        --samplesheet '${samplesheet}' --quants . \
        --tx2gene '${tx2gene}' --reference '${reference}' \
        --comparison '${comparison}' --out differential_expression_summary.csv
    """
    stub:
    """
    touch differential_expression_summary.csv normalized_counts.csv MA_plot.pdf PCA.pdf analysis_session.txt
    """
}
