# RNA-seq workflow: quantification and optional differential expression

A Nextflow DSL2 portfolio project connecting FastQC, Salmon and MultiQC.
An optional **two-group, unpaired** analysis imports Salmon estimates with
tximport and tests gene-level differences with DESeq2.

## Status and scope

This is a workflow demonstration under validation, not a validated production
or clinical pipeline. No real-cohort biological findings are claimed.
The bundled tiny FASTQ fixture has one sample per condition: it exercises
workflow wiring and quantification only, not differential expression.
Preview/stub checks do not validate numerical results.

The previous script labeled DESeq2 only averaged abundances and could emit
dummy results. It has been replaced; missing inputs now fail explicitly.

## Quantification quick start

Install a compatible Nextflow/Java runtime and Docker, or Conda.
From a clone:

```bash
nextflow run main.nf -profile test -stub-run
nextflow run main.nf -profile test,docker
```

For your own paired-end reads, supply a CSV with
`sample,fastq_1,fastq_2,condition`, an appropriate reference transcriptome, and
a destination:

```bash
nextflow run main.nf -profile docker --input samplesheet.csv \
  --transcriptome transcripts.fa --outdir results
```

Use unique, filesystem-safe sample IDs. Use absolute FASTQ paths for your own
data. Keep reference accessions, releases and checksums with the study metadata.
Use `-resume` to reuse completed tasks.

## Optional differential expression

DE currently requires the Conda profile; Docker/Singularity execution of the DE
stage is not configured. The R environment includes DESeq2 and tximport.

Supply at least two **independent biological replicates** per condition (three
or more is a better starting point for many studies). Technical replicates do
not establish biological independence. This implementation supports exactly
two conditions and the design `~ condition`; it does not handle paired,
batch-adjusted or other multifactor designs.

Create a tab-separated mapping with the exact header `transcript<TAB>gene`.
Every transcript in the Salmon tables must map to one gene, using the same
reference release and exact IDs, including version suffixes.

```bash
nextflow run main.nf -profile conda --input samplesheet.csv \
  --transcriptome transcripts.fa --run_de true --tx2gene tx2gene.tsv \
  --reference control --comparison treated --outdir results
```

Positive log2 fold changes mean higher expression in comparison versus
reference. Genes with total estimated count below 10 are removed. DESeq2 uses
an alpha of 0.05 and BH adjusted p-values; output fold changes are unshrunken.
NA adjusted p-values must not be counted as significant.

Outputs in `results/deseq2/`:
- `differential_expression_summary.csv`: gene IDs, effect sizes and statistics.
- `normalized_counts.csv`: normalized counts for inspection.
- `PCA.pdf`, `MA_plot.pdf`: exploratory diagnostics.
- `analysis_session.txt`: contrast, assumptions and package/session versions.

Inspect QC, confounding and sample independence before interpreting results.
A significant gene is not automatically a validated biomarker.

## Validation

The existing CI checks workflow preview and stubs. The analysis CI separately
runs real tximport/DESeq2 on synthetic Salmon tables with eight samples and
500 genes, checks an introduced effect and output artifacts, and verifies
that missing inputs and absent replication fail without creating results:

```bash
Rscript tests/test_deseq2.R
```

This tests the R analysis, not the full FASTQ-to-report workflow. A replicated
public-cohort case study, full real-tool workflow integration, pinned dependency
resolution, pathway analysis, and measured resource usage remain future work.
No HPC/cloud execution has been validated here. A container profile alone
does not configure SLURM or a cloud executor.

## Methods references

- [tximport documentation](https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html)
- [DESeq2 documentation](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)

Code: MIT license. Author: Gretta Y. Rozenberg.
