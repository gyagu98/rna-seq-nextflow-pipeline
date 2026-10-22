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
    ch_quants = SALMON_QUANT.out.quant_sf.map { meta, sf -> sf }.collect()
    ch_samplesheet = file(params.input).exists() ? file(params.input) : file("${projectDir}/${params.input}")
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
