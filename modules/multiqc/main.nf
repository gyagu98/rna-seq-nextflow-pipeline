process MULTIQC {
    label 'process_single'
    publishDir "${params.outdir}/multiqc", mode: 'copy'
    container 'quay.io/biocontainers/multiqc:1.21--pyhdfd78af_0'

    input:
    path ('*')

    output:
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data

    script:
    """
    multiqc .
    """
}
