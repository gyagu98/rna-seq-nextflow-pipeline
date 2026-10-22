#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
========================================================================================
   Gretta Y. Rozenberg - rna-seq-nextflow-pipeline
   A reproducible RNA-seq pipeline using Nextflow DSL2, FastQC, Salmon, and MultiQC
========================================================================================
*/

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
            def r1 = file(row.fastq_1).exists() ? file(row.fastq_1) : file("${projectDir}/${row.fastq_1}")
            def r2 = file(row.fastq_2).exists() ? file(row.fastq_2) : file("${projectDir}/${row.fastq_2}")
            return [meta, [r1, r2]]
        }
        .set { ch_reads }

    // 2. Reference Transcriptome Channel
    def transcriptome_file = file(params.transcriptome).exists() ? file(params.transcriptome) : file("${projectDir}/${params.transcriptome}")
    ch_transcriptome = Channel.value(transcriptome_file)

    // 3. Quality Control (FastQC)
    FASTQC(ch_reads)
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.map { meta, zip -> zip })

    // 4. Salmon Indexing
    SALMON_INDEX(ch_transcriptome, params.salmon_kmer)

    // 5. Salmon Quantification
    SALMON_QUANT(ch_reads, SALMON_INDEX.out.index)
    ch_multiqc_files = ch_multiqc_files.mix(SALMON_QUANT.out.results.map { meta, dir -> dir })

    // 6. Differential Expression (DESeq2)
    ch_quants = SALMON_QUANT.out.results.map { meta, dir -> dir }.collect()
    ch_samplesheet = file(params.input).exists() ? file(params.input) : file("${projectDir}/${params.input}")
    if (params.run_de) {
        if (!workflow.profile.tokenize(',').contains('conda'))
            error 'DE currently requires -profile conda (or test,conda)'
        if (!params.tx2gene || !params.reference || !params.comparison)
            error 'DE requires --tx2gene, --reference and --comparison'
        if (!(params.reference ==~ /[A-Za-z0-9_.-]+/) || !(params.comparison ==~ /[A-Za-z0-9_.-]+/))
            error 'Condition names must contain only letters, digits, underscores, dots and hyphens'
        DESEQ2(ch_quants, ch_samplesheet, file(params.tx2gene, checkIfExists: true),
               params.reference, params.comparison)
    }

    // 7. Aggregate Quality Control (MultiQC)
    MULTIQC(ch_multiqc_files.collect())
}

