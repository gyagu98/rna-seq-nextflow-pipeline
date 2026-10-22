#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
========================================================================================
   Gretta Y. Rozenberg - rna-seq-nextflow-pipeline
   A reproducible RNA-seq pipeline using Nextflow DSL2, FastQC, Salmon, and MultiQC
========================================================================================
*/

def helpMessage() {
    log.info"""
    Usage:
      nextflow run gyagu98/rna-seq-nextflow-pipeline [options]

    Mandatory arguments:
      --input             Path to input samplesheet CSV (columns: sample,fastq_1,fastq_2,condition)
      --transcriptome     Path to reference transcriptome FASTA

    Optional arguments:
      --outdir            The output directory where the results will be saved (default: 'results')
      --salmon_kmer       K-mer length for Salmon indexing (default: 21)
      -profile            Configuration profile to use: docker, singularity, conda, test (default: 'standard')
    """.stripIndent()
}

if (params.help) {
    helpMessage()
    exit 0
}

if (!params.input && !params.test) {
    exit 1, "ERROR: Please specify an input samplesheet with '--input samplesheet.csv' or run with '-profile test'"
}
if (!params.transcriptome && !params.test) {
    exit 1, "ERROR: Please specify a reference transcriptome with '--transcriptome reference.fasta' or run with '-profile test'"
}

/*
========================================================================================
   IMPORT LOCAL MODULES
========================================================================================
*/
include { FASTQC }       from './modules/fastqc/main'
include { SALMON_INDEX } from './modules/salmon/index'
include { SALMON_QUANT } from './modules/salmon/quant'
include { MULTIQC }      from './modules/multiqc/main'
include { DESEQ2 }       from './modules/deseq2/main'

/*
========================================================================================
   WORKFLOW DEFINITION
========================================================================================
*/
workflow {
    ch_multiqc_files = Channel.empty()

    // 1. Channel from samplesheet CSV
    Channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map { row ->
            def meta = [id: row.sample, condition: row.condition]
            def reads = [file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true)]
            return [meta, reads]
        }
        .set { ch_reads }

    // 2. Reference Transcriptome Channel
    ch_transcriptome = file(params.transcriptome, checkIfExists: true)

    // 3. Quality Control (FastQC)
    FASTQC(ch_reads)
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.map { meta, zip -> zip })

    // 4. Salmon Indexing
    SALMON_INDEX(ch_transcriptome, params.salmon_kmer)

    // 5. Salmon Quantification
    SALMON_QUANT(ch_reads, SALMON_INDEX.out.index)
    ch_multiqc_files = ch_multiqc_files.mix(SALMON_QUANT.out.results.map { meta, dir -> dir })

    // 6. Differential Expression (DESeq2)
    ch_quants = SALMON_QUANT.out.quant_sf.map { meta, sf -> sf }.collect()
    ch_samplesheet = file(params.input)
    DESEQ2(ch_quants, ch_samplesheet)

    // 7. Aggregate Quality Control (MultiQC)
    MULTIQC(ch_multiqc_files.collect())
}

workflow.onComplete {
    log.info """
    Pipeline completed at : ${workflow.complete}
    Execution status     : ${workflow.success ? 'SUCCESS' : 'FAILED'}
    Duration             : ${workflow.duration}
    Output directory     : ${params.outdir}
    """.stripIndent()
}
