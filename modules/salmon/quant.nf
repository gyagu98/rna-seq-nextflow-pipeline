process SALMON_QUANT {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}/salmon_quant", mode: 'copy'
    container 'quay.io/biocontainers/salmon:1.10.2--hecfa306_0'

    input:
    tuple val(meta), path(reads)
    path index

    output:
    tuple val(meta), path("${meta.id}"), emit: results
    tuple val(meta), path("${meta.id}/quant.sf"), emit: quant_sf

    script:
    """
    salmon quant \
        -i ${index} \
        -l A \
        -1 ${reads[0]} \
        -2 ${reads[1]} \
        -o ${meta.id} \
        --threads ${task.cpus} \
        --validateMappings
    """

    stub:
    """
    mkdir -p ${meta.id}
    touch ${meta.id}/quant.sf
    """
}
