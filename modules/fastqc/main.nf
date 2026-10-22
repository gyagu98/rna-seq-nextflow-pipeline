process FASTQC {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}/fastqc", mode: 'copy'
    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip") , emit: zip

    script:
    """
    fastqc -q ${reads}
    """

    stub:
    """
    touch ${meta.id}_1_fastqc.html
    touch ${meta.id}_1_fastqc.zip
    touch ${meta.id}_2_fastqc.html
    touch ${meta.id}_2_fastqc.zip
    """
}
