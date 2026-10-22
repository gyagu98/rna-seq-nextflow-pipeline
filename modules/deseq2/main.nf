process DESEQ2 {
    label 'process_low'
    publishDir "${params.outdir}/deseq2", mode: 'copy'
    container 'rocker/tidyverse:4.3.2'

    input:
    path quant_dirs
    path samplesheet

    output:
    path "differential_expression_summary.csv", emit: dge_summary

    script:
    """
    Rscript ${projectDir}/bin/deseq2_analysis.R \
        --samplesheet ${samplesheet} \
        --quants . \
        --out differential_expression_summary.csv
    """

    stub:
    """
    touch differential_expression_summary.csv
    """
}
