process SALMON_INDEX {
    tag "$fasta"
    label 'process_medium'
    publishDir "${params.outdir}/salmon_index", mode: 'copy'
    container 'quay.io/biocontainers/salmon:1.10.2--hecfa306_0'

    input:
    path fasta
    val kmer

    output:
    path "salmon_index", emit: index

    script:
    """
    salmon index -t ${fasta} -i salmon_index -k ${kmer} --threads ${task.cpus}
    """
}
